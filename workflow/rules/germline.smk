# Functions and rules for calling germline variants
from scripts.common import abstract_location, allocated


rule deepvariant:
    """
    Data processing step to call germline variants using deep neural 
    network. DeepVariant is a deep learning-based variant caller that 
    takes aligned reads (in BAM or CRAM format), produces pileup image 
    tensors from them, classifies each tensor using a convolutional 
    neural network, and finally reports the results in a standard VCF 
    or gVCF file.
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
        tmpdir = join(workpath, "deepvariant", "tmp", "{name}"),
    message: "Running DeepVariant on '{input.bam}' input file"
    threads: int(allocated("threads", "deepvariant", cluster))
    container: config['images']['deepvariant']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    run_deepvariant \\
        --model_type=WGS \\
        --ref={params.genome} \\
        --reads={input.bam} \\
        --output_gvcf={output.gvcf} \\
        --output_vcf={output.vcf} \\
        --num_shards={threads} \\
        --intermediate_results_dir=${{tmp}}
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
        Multi-sample joint called BCF file
    """
    input: 
        gvcf = expand(join(workpath,"deepvariant","gVCFs","{name}.g.vcf.gz"), name=samples),
    output:
        gvcfs = join(workpath, "deepvariant", "VCFs", "gvcfs.list"),
        bcf   = join(workpath, "deepvariant", "VCFs", "joint.bcf"),
        norm  = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
    params: 
        rname  = "glnexus",
        gvcfdir = join(workpath, "deepvariant", "gVCFs"),
        memory  = allocated("mem", "glnexus", cluster).rstrip('G'),
        genome = config['references']['GENOME'],
    message: "Running GLnexus on a set of gVCF files"
    threads: int(allocated("threads", "glnexus", cluster))
    container: config['images']['glnexus']
    shell: """
    # Avoids ARG_MAX issue which will
    # limit max length of a command
    find {params.gvcfdir} -iname '*.g.vcf.gz' \\
    > {output.gvcfs}

    glnexus_cli \\
        --config DeepVariant \\
        --list {output.gvcfs} \\
        --threads {threads} \\
        --mem-gbytes {params.memory} \\
    > {output.bcf}

    bcftools norm \\
        -Oz \\
        --threads {threads} \\
        -f {params.genome} \\
        -o {output.norm} \\
        {output.bcf}

    bcftools index \\
        -f -t \\
        --threads {threads} \\
        {output.norm}
    """
