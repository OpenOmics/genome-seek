# Functions and rules for calling Copy Number Variation
from scripts.common import abstract_location, allocated


rule peddy:
    """
    Peddy compares familial-relationships and sexes as reported in 
    a PED/FAM file with those inferred from a VCF. It samples the VCF 
    at about 25000 sites (plus chrX) to accurately estimate relatedness, 
    IBS0, heterozygosity, sex and ancestry. Note that somalier is a more 
    scalable, faster, replacement for peddy that uses some of the same 
    methods as peddy along with some new ones; however, we are passing 
    information from peddy to CANVAS.
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        CSV file information about ped-reported & genotype-inferred sex
    """
    input:
        vcf = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
    output:
        ped  = temp(join(workpath, "deepvariant", "VCFs", "{batch}_batch.ped")),
        html = join(workpath, "deepvariant", "VCFs", "{batch}_peddy.html"),
        csv  = join(workpath, "deepvariant", "VCFs", "{batch}_peddy.sex_check.csv"),
    params:
        rname = "peddy",
        prefix = join(workpath, "deepvariant", "VCFs", "{batch}_peddy"),
        peddy_chr = config['references']['PEDDY_FILTER'],
        intermediate = join(workpath, "deepvariant", "VCFs", "intermediate")
    message: "Running peddy on '{input.vcf}' input file"
    threads: int(allocated("threads", "peddy", cluster))
    envmodules: 
        config['tools']['vcftools'],
        config['tools']['peddy']
    shell: """
    vcftools \\
        --gzvcf {input.vcf} \\
        --plink \\
        --out {params.intermediate} \\
        --chr {params.peddy_chr}
    cut -f1-6 {params.intermediate}.ped \\
    > {output.ped}
    peddy -p {threads} \\
        --prefix {params.prefix} \\
        {input.vcf} \\
        {output.ped}
    """
