# Functions and rules for calling germline variants
from scripts.common import (
    abstract_location, 
    allocated,
    provided
)

rule deepvariant:
    """
    Data processing step to call germline variants using deep neural 
    network. DeepVariant is a deep learning-based variant caller that 
    takes aligned reads (in BAM or CRAM format), produces pileup image 
    tensors from them, classifies each tensor using a convolutional 
    neural network, and finally reports the results in a standard VCF 
    or gVCF file. This rule runs all three steps in the deepvariant
    pipeline as a single step, i.e.: make_examples, call_variants, and
    postprocess_variants. This is not optimal for large-scale projects
    as it will consume a lot of resources inefficently (only the 2nd
    step in the dv pipeline can make use of GPU-computing). As so, it
    is better to run the 1st/3rd step on a normal compute node and run
    the 2nd step on a GPU node. This rule is depreciated. Please see 
    the deepvariant_makeexamples, deepvariant_callvariants, and
    deepvariant_postprocessvariants rules for the optimal way to
    run this tool.
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Single-sample gVCF file with called variants
    """
    input: 
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai"),
    output:
        gvcf = join(workpath, "deepvariant", "gVCFs", "{name}.g.vcf.gz"),
        vcf  = join(workpath, "deepvariant", "VCFs", "{name}.vcf.gz"),
    params: 
        rname  = "deepvar",
        genome = config['references']['GENOME'],
        tmpdir = tmpdir,
        # Building option for glnexus config, where:
        #  @WES = --model_type=WES
        #  @WGS = --model_type=WGS
        dv_model_type = lambda _: "WES" if run_wes else "WGS",
    message: "Running DeepVariant on '{input.bam}' input file"
    threads: int(allocated("threads", "deepvariant", cluster))
    container: config['images']['deepvariant']
    envmodules: config['tools']['deepvariant']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    run_deepvariant \\
        --model_type={params.dv_model_type} \\
        --ref={params.genome} \\
        --reads={input.bam} \\
        --output_gvcf={output.gvcf} \\
        --output_vcf={output.vcf} \\
        --num_shards={threads} \\
        --intermediate_results_dir=${{tmp}}
    """


rule deepvariant_make_examples:
    """
    Data processing step to call germline variants using deep neural 
    network. The make_examples step prepares the input data for the
    deepvariant's CNN. DeepVariant is a deep learning-based variant 
    caller composed of multiple steps that takes aligned reads (in 
    BAM or CRAM format), produces pileup image tensors from them, 
    classifies each tensor using a convolutional neural network, 
    and finally reports the results in a standard VCF or gVCF file. 
    This rule is the first step in the deepvariant pipeline:
      1. make_examples        (CPU, parallelizable with gnu-parallel)
      2. call_variants        (GPU, use a GPU node)
      3. postprocess_variants (CPU)
    Running deepvariant in a single step using run_deepvariant is not 
    optimal for large-scale projects as it will consume resources very
    inefficently. As so, it is better to run the 1st/3rd step on a compute 
    node and run the 2nd step on a GPU node.
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Single-sample gVCF file with called variants
    """
    input: 
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai"),
    output:
        success = join(workpath, "deepvariant", "mk_examples", "{name}.make_examples.success"),
    params: 
        rname  = "dv_mkexamples",
        genome = config['references']['GENOME'],
        tmpdir = tmpdir,
        nshards = int(allocated("threads", "deepvariant", cluster))-1,
        example = lambda w: join(workpath, "deepvariant", "mk_examples", "{0}.make_examples.tfrecord@{1}.gz".format(
            w.name,
            int(allocated("threads", "deepvariant", cluster))
        )),
        gvcf = lambda w: join(workpath, "deepvariant", "mk_examples", "{0}.gvcf.tfrecord@{1}.gz".format(
            w.name,
            int(allocated("threads", "deepvariant", cluster))
        )),
    message: "Running DeepVariant make_examples on '{input.bam}' input file"
    threads: int(allocated("threads", "deepvariant", cluster))
    container: config['images']['deepvariant']
    envmodules: config['tools']['deepvariant']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    # Run DeepVariant make_examples and
    # parallelize it using gnu-parallel 
    time seq 0 {params.nshards} \\
        | parallel \\
            --eta \\
            -q \\
            --halt 2 \\
            --line-buffer \\
            make_examples \\
                --mode calling \\
                --ref {params.genome} \\
                --reads {input.bam} \\
                --examples {params.example} \\
                --channels "insert_size" \\
                --gvcf {params.gvcf} \\
                --task {{}} \\
    && touch {output.success}
    """


rule glnexus:
    """
    Data processing step to merge and joint call a set of gVCF files.
    GLnexus is a scalable gVCF merging and joint variant calling for 
    population-scale sequencing projects. For more information about
    GLnexus & deepvariant, please check out this comparison with GATK: 
    https://academic.oup.com/bioinformatics/article/36/24/5582/6064144
    @Input:
        Set of gVCF files (gather)
    @Output: 
        Multi-sample jointly called, normalized VCF file
    """
    input: 
        gvcf = expand(join(workpath,"deepvariant","gVCFs","{name}.g.vcf.gz"), name=samples),
        bed  = provided(join(workpath, "references", "wes_regions_50bp_padded.bed"), run_wes),
    output:
        gvcfs = join(workpath, "deepvariant", "VCFs", "gvcfs.list"),
        bcf   = join(workpath, "deepvariant", "VCFs", "joint.bcf"),
        norm  = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
        jvcf  = join(workpath, "deepvariant", "VCFs", "joint.glnexus.vcf.gz"),
    params: 
        rname  = "glnexus",
        tmpdir =  tmpdir,
        gvcfdir = join(workpath, "deepvariant", "gVCFs"),
        memory  = allocated("mem", "glnexus", cluster).rstrip('G'),
        genome = config['references']['GENOME'],
        # Building option for glnexus config, where:
        #  @WES = --config DeepVariantWES
        #  @WGS = --config DeepVariant_unfiltered
        gl_config = lambda _: "DeepVariantWES" if run_wes else "DeepVariant_unfiltered",
        # Building option for GLnexus bed file:
        #  @WES = --bed wes_regions_50bp_padded.bed
        #  @WGS = '' 
        wes_bed_option = lambda _: "--bed {0}".format(
            join(workpath, "references", "wes_regions_50bp_padded.bed"),
        ) if run_wes else "",
    message: "Running GLnexus on a set of gVCF files"
    threads: int(allocated("threads", "glnexus", cluster))
    container: config['images']['glnexus']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit, 
    # GLnexus tmpdir should NOT exist
    # prior to running it. If it does
    # exist prior to runnning, it will
    # immediately error out.
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp_parent=$(mktemp -d -p "{params.tmpdir}")
    tmp_dne=$(echo "${{tmp_parent}}"| sed 's@$@/GLnexus.DB@')
    trap 'rm -rf "${{tmp_parent}}"' EXIT

    # Avoids ARG_MAX issue which will
    # limit max length of a command
    find {params.gvcfdir} -iname '*.g.vcf.gz' \\
    > {output.gvcfs}

    glnexus_cli \\
        --dir ${{tmp_dne}} \\
        --config {params.gl_config} {params.wes_bed_option} \\
        --list {output.gvcfs} \\
        --threads {threads} \\
        --mem-gbytes {params.memory} \\
    > {output.bcf}

    bcftools norm \\
        -m - \\
        -Oz \\
        --threads {threads} \\
        -f {params.genome} \\
        -o {output.norm} \\
        {output.bcf}

    bcftools view \\
        -Oz \\
        --threads {threads} \\
        -o {output.jvcf} \\
        {output.bcf}

    bcftools index \\
        -f -t \\
        --threads {threads} \\
        {output.norm}
    
    bcftools index \\
        -f -t \\
        --threads {threads} \\
        {output.jvcf}
    """


rule gatk_selectvariants:
    """
    Make per-sample jointly called, normalized VCFs from GLnexus output.
    @Input:
        Multi-sample joint, normalized VCF file (indirect-gather-due-to-aggregation)
    @Output:
        Single-sample joint, normalized VCF file
    """
    input: 
        vcf = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
        bed = provided(wes_bed_file, run_wes),
    output: 
        vcf = join(workpath, "deepvariant", "VCFs", "{name}.germline.vcf.gz"),
    params: 
        rname  = "varselect",
        genome = config['references']['GENOME'], 
        sample = "{name}", 
        memory = allocated("mem", "gatk_selectvariants", cluster).rstrip('G'),
        # Building WES options for gatk selectvariants, where:
        #  @WES = --intervals gencode_v44_protein-coding_exons.bed -ip 100
        #  @WGS = ''
        wes_bed_option  = lambda _: "--intervals {0}".format(wes_bed_file) if run_wes else "",
        wes_bed_padding = lambda _: "-ip 100" if run_wes else "",
    message: "Running GATK4 SelectVariants on '{input.vcf}' input file"
    threads: int(allocated("threads", "gatk_selectvariants", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    gatk --java-options '-Xmx{params.memory}g -XX:ParallelGCThreads={threads}' SelectVariants \\
        -R {params.genome} {params.wes_bed_option} {params.wes_bed_padding} \\
        --variant {input.vcf} \\
        --sample-name {params.sample} \\
        --exclude-non-variants \\
        --output {output.vcf}
    """
