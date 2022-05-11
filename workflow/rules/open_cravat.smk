# Functions and rules for running OpenCRAVAT
from scripts.common import abstract_location, allocated


rule chrom_selectvariants:
    """
    Makes per-chromosome jointly called, normalized VCFs from GLnexus output.
    This step is used to speed up the remaining OpenCRAVAT rules. At the end, 
    the resulting SQL-lite databases are merged back together (gather-step) 
    for each chromosome.
    @Input:
        Multi-sample joint, normalized VCF file (scatter)
    @Output:
        Chromosome chunked joint, normalized VCF file
    """
    input: 
        vcf = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
    output: 
        vcf = join(workpath, "OpenCRAVAT", "VCFs", "{chunk}.germline.vcf.gz"),
    params: 
        rname  = "chrselect",
        genome = config['references']['GENOME'], 
        chrom  = "{chunk}", 
        memory = allocated("mem", "chrom_selectvariants", cluster).rstrip('G')
    message: "Running Chromosome SelectVariants on '{input.vcf}' input file"
    threads: int(allocated("threads", "chrom_selectvariants", cluster))
    envmodules: 
        config['tools']['gatk4'], config['tools']['bcftools']
    shell: """
    gatk --java-options '-Xmx{params.memory}g -XX:ParallelGCThreads={threads}' SelectVariants \\
        -R {params.genome} \\
        --variant {input.vcf} \\
        -L {params.chrom} \\
        --exclude-non-variants \\
        --output {output.vcf}
    
    tabix --force \\
        -p vcf {output.vcf}
    """


rule open_cravat:
    """
    Performs genomic variant interpretation including variant impact, annotation,
    and scoring using OpenCRAVAT. Creates a SQL-lite database that can be used 
    with OpenCRAVAT's user interface. Here is more information about OpenCRAVAT: 
    https://open-cravat.readthedocs.io/en/latest/
    @Input:
        Chromosome chunked joint, normalized VCF file (scatter)
    @Output:
        Chromosome chunked joint, normalized VCF file
    """
    input: 
        vcf = join(workpath, "OpenCRAVAT", "VCFs", "{chunk}.germline.vcf.gz"),
    output: 
        db = join(workpath, "OpenCRAVAT", "cravat_{chunk}.sqlite"),
    params: 
        rname  = "cravat",
        prefix = "cravat_{chunk}",
        outdir = join(workpath, "OpenCRAVAT"),
        annot  = ' '.join(config['options']['oc_annotators']),
        genome = config['references']['OC_LIFTOVER'],
        module = config['references']['OC_MODULES'],
    message: "Running OpenCRAVAT run on '{input.vcf}' input file"
    threads: int(allocated("threads", "open_cravat", cluster))
    envmodules: config['tools']['open_cravat']
    shell: """
    oc run \\
        -t vcf \\
        -x \\
        --newlog \\
        --cleanrun \\
        --system-option "modules_dir={params.module}" \\
        -a {params.annot} \\
        -n {params.prefix} \\
        -l {params.genome} \\
        -d {params.outdir} \\
        --mp {threads} \\
        {input.vcf}
    """
