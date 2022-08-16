# Functions and rules for calling Copy Number Variation
from scripts.common import abstract_location, allocated

# Helper functions for tumor, normal pairs 
def get_normal_recal_bam(wildcards):
    """
    Returns a tumor samples paired normal
    See config['pairs'] for tumor, normal pairs.
    """
    normal = tumor2normal[wildcards.name]
    if normal:
        # Runs in a tumor, normal mode
        return join(workpath, "BAM", "{0}.recal.bam".format(normal))
    else:
        # Runs in tumor-only mode
        return []

# Germline Copy Number Variation
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


rule canvas:
    """
    Canvas is a tool for calling copy number variants (CNVs). Germline-WGS 
    will call CNVs in germline samples from whole genome sequencing data.
    It is worth noting that Germline-WGS mode has been deprecated. The 
    authors recommend running SmallPedigree-WGS instead, even for a single 
    sample analysis. In the near future, we may move to this option.
    For more information about CANVAS, please read through the paper:
    https://academic.oup.com/bioinformatics/article/32/15/2375/1743834
    Predicted CNVs from CANVAS are then filtered, annotated  and ranked
    with AnnotSV. AnnotSV Github: https://github.com/lgmgeo/AnnotSV
    @Input:
        Per sample VCF file (scatter)
    @Output:
        Per sample VCF with annotated and ranked CNVs
    """
    input:
        vcf = join(workpath, "deepvariant", "VCFs", "{name}.vcf.gz"),
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai"),
        csv = join(workpath, "deepvariant", "VCFs", "{0}_peddy.sex_check.csv".format(batch_id)),
        joint = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
    output:
        vcf = join(workpath, "CANVAS", "{name}", "CNV.vcf.gz"),
        ploidy = join(workpath, "CANVAS", "{name}", "ploidy.vcf"),
        filtered  = join(workpath, "CANVAS", "{name}", "CNV.filtered.vcf"),
        annotated = join(workpath, "CANVAS", "{name}", "{name}.segments.annotations.tsv"),
    params:
        rname  = "canvas",
        sample = "{name}",
        outdir = join(workpath, "CANVAS", "{name}"),
        checkpoints   = join(workpath, "CANVAS", "{name}", "Checkpoints"),
        male_ploidy   = join(workpath, "resources", "male_ploidy.vcf"),
        female_ploidy = join(workpath, "resources", "female_ploidy.vcf"),
        canvas_filter = config['references']['CANVAS_FILTER'],
        canvas_kmer   = config['references']['CANVAS_KMER'],
        canvas_genome = config['references']['CANVAS_GENOME'],
        canvas_balleles = config['references']['CANVAS_BALLELES'],
        annotsv_build = config['references']['ANNOTSV_BUILD'],
    message: "Running canvas on '{input.vcf}' input file"
    threads: int(allocated("threads", "canvas", cluster))
    envmodules: 
        config['tools']['canvas'],
        config['tools']['bcftools'],
        config['tools']['annotsv'],
    shell: """
    # Get biological sex
    # predicted by peddy 
    predicted_sex=$(awk -F ',' \\
        '$8=="{params.sample}" \\
        {{print $7}}' \\
        {input.csv}
    )

    # Copy over base ploidy 
    # vcf for predicted sex 
    # and add name to header
    if [ "$predicted_sex" == "male" ]; then 
        cp {params.male_ploidy} {output.ploidy}
    else
        # prediction is female
        cp {params.female_ploidy} {output.ploidy}
    fi
    sed -i 's/SAMPLENAME/{params.sample}/g' \\
        {output.ploidy}

    # Delete Canvas checkpoints
    if [ -d {params.checkpoints} ]; then
        # Forces Canvas to start 
        # over from the beginning
        rm -rf '{params.checkpoints}'
    fi

    # CANVAS in Germline WGS mode
    export COMPlus_gcAllowVeryLargeObjects=1
    Canvas.dll Germline-WGS \\
        -b {input.bam} \\
        -n {params.sample} \\
        -o {params.outdir} \\
        -r {params.canvas_kmer} \\
        --ploidy-vcf={output.ploidy} \\
        -g {params.canvas_genome} \\
        -f {params.canvas_filter} \\
        --sample-b-allele-vcf={input.vcf}
    
    # Filter predicted CNVs
    bcftools filter \\
        --include 'FILTER="PASS" && INFO/SVTYPE="CNV"' \\
        {output.vcf} \\
    > {output.filtered}

    # Rank and annotate CNVs
    AnnotSV \\
        -genomeBuild {params.annotsv_build} \\
        -outputDir {params.outdir} \\
        -outputFile {output.annotated} \\
        -SVinputFile {output.filtered} \\
        -vcfFiles {input.joint} \\
        -vcfSamples {params.sample}
    
    # Check if AnnotSV silently failed
    if [ ! -f "{output.annotated}" ]; then
        # AnnotSV failed to process
        # provided SVinputFile file, 
        # usually due to passing an 
        # empty filtered SV file
        echo "WARNING: AnnotSV silently failed..." 1>&2
        touch {output.annotated}
    fi
    """


# Somatic Copy Number Variation
rule hmftools_amber:
    """Data-processing step to generate a tumor BAF file that is required
    for Purple's copy number fitting. HMF Tools is a suite of tools the
    Hartwig Medical Foundation developed to analyze genomic data. Amber 
    and cobalt must be run prior to running purple. For more information 
    about hmftools visit: https://github.com/hartwigmedical/hmftools
    @Input:
       Realigned, recalibrated BAM file (scatter-per-tumor-sample)
    @Output:
        BAF file for purple's copy number fit
    """
    input:
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        snp = join(workpath, "hmftools", "amber", "{name}.amber.snp.vcf.gz"),
        con = join(workpath, "hmftools", "amber", "{name}.amber.contamination.vcf.gz"),
        baf = join(workpath, "hmftools", "amber", "{name}.amber.baf.tsv"),
    params:
        rname     = 'hmfamber',
        tumor     = '{name}',
        outdir    = join(workpath, "hmftools", "amber"),
        amber_jar = config['references']['HMFTOOLS_AMBER_JAR'],
        loci_ref  = config['references']['HMFTOOLS_AMBER_LOCI'],
        memory    = allocated("mem", "hmftools_amber", cluster).lower().rstrip('g'),
        # Building optional argument for paired normal
        normal_name = lambda w: "-reference {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
        normal_bam = lambda w: "-reference_bam {0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
        # Building optional flag for tumor-only
        tumor_flag = lambda w: "" if tumor2normal[w.name] else "-tumor_only",
    threads: 
        int(allocated("threads", "hmftools_amber", cluster)),
    envmodules:
        config['tools']['rlang'],
    shell: """
    java -Xmx{params.memory}g -cp {params.amber_jar} \\
        com.hartwig.hmftools.amber.AmberApplication \\
            -tumor {params.tumor} {params.normal_name} \\
            -tumor_bam {input.tumor}  {params.normal_bam}\\
            -output_dir {params.outdir} \\ 
            -threads {threads} {params.tumor_flag} \\
            -loci {params.loci_ref}
    """
