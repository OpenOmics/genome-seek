# Functions and rules for quality-control
from scripts.common import (
    abstract_location, 
    allocated,
    ignore,
    provided 
)

# Pre alignment QC-related rules
rule fc_lane:
    """
    Quality-control step to get flowcell and lane information from FastQ file.
    FastQ files generated with older versions of Casava or downloaded from
    SRA have a different format than newer FastQ files generated with the
    current version of Casava. It is worth noting that FastQ files downloaded 
    from SRA or FastQ files generated with Casava version < 1.8 do not have 
    Flowcell IDs in its sequence indentifer. If a FastQ file does not have 
    Flowcell IDs, the Machine or Instrument ID is grabbed instead. This rule
    relies on python2; however, it could easily be modified to run in both
    majors versions of python (python2.X and python3.X). 
    @Input:
        Raw FastQ R1 file (scatter)
    @Output:
        Text file containing information about the FastQ file
    """
    input:
        r1 = join(workpath, "{name}.R1.fastq.gz"),
    output:
        txt = join(workpath,"rawQC","{name}.fastq.info.txt")
    params:
        rname = 'fclane',
        get_flowcell_lanes = join("workflow", "scripts", "flowcell_lane.py"),
    threads: int(allocated("threads", "fc_lane", cluster))
    envmodules: config['tools']['python2']
    shell: """
    python {params.get_flowcell_lanes} \\
        {input.r1} \\
        {wildcards.name} \\
    > {output.txt}
    """


rule fastqc_raw:
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
        r1 = join(workpath,"{name}.R1.fastq.gz"),
        r2 = join(workpath,"{name}.R2.fastq.gz"),
    output:
        join(workpath,"rawQC","{name}.R1_fastqc.zip"),
        join(workpath,"rawQC","{name}.R2_fastqc.zip"),
    params:
        rname  = 'rawfqc',
        outdir = join(workpath,"rawQC"),
    threads: int(allocated("threads", "fastqc_raw", cluster))
    envmodules: config['tools']['fastqc']
    shell: """
    fastqc \\
        {input.r1} \\
        {input.r2} \\
        -t {threads} \\
        -o {params.outdir}
    """


rule fastq_screen:
    """
    Quality-control step to screen for different sources of contamination.
    FastQ Screen compares your sequencing data to a set of different reference
    genomes to determine if there is contamination. It allows a user to see if
    the composition of your library matches what you expect.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        FastQ Screen report and logfiles
    """
    input:
        fq1 = join(workpath,"fastqs","{name}.R1.trimmed.fastq.gz"),
        fq2 = join(workpath,"fastqs","{name}.R2.trimmed.fastq.gz")
    output:
        txt1 = join(workpath,"FQscreen","{name}.R1.trimmed_screen.txt"),
        txt2 = join(workpath,"FQscreen","{name}.R2.trimmed_screen.txt"),
        png1 = join(workpath,"FQscreen","{name}.R1.trimmed_screen.png"),
        png2 = join(workpath,"FQscreen","{name}.R2.trimmed_screen.png")
    params:
        rname  = "fqscreen",
        outdir = join(workpath,"FQscreen"),
        # Exposed Parameters: modify resources/fastq_screen.conf to change 
        # default locations to bowtie2 indices
        fastq_screen_config = config['references']['FASTQ_SCREEN_CONFIG'],
    envmodules: 
        config['tools']['fastq_screen'],
        config['tools']['bowtie']
    threads: int(allocated("threads", "fastq_screen", cluster))
    shell: """
    fastq_screen --conf {params.fastq_screen_config} \\
        --outdir {params.outdir} \\
        --threads {threads} \\
        --subset 1000000 \\
        --aligner bowtie2 \\
        --force \\
        {input.fq1} \\
        {input.fq2}
    """

# Post-alignment QC-related rules
rule fastqc_bam:
    """
    Quality-control step to assess sequencing quality of each sample.
    FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        FastQC report and zip file containing sequencing quality information
    """
    input:
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
    output:
        zipfile =  join(workpath,"QC","{name}.sorted_fastqc.zip"),
        report  =  join(workpath,"QC","{name}.sorted_fastqc.html")
    params:
        rname  = "bamfqc",
        outdir =  join(workpath,"QC"),
        adapters = config['references']['FASTQC_ADAPTERS']
    message: "Running FastQC with {threads} threads on '{input.bam}' input file"
    threads: int(allocated("threads", "fastqc_bam", cluster))
    envmodules: config['tools']['fastqc']
    shell: """
    fastqc -t {threads} \\
        -f bam \\
        --contaminants {params.adapters} \\
        -o {params.outdir} \\
        {input.bam} 
    """


rule qualimap:
    """
    Quality-control step to assess various post-alignment metrics 
    and a secondary method to calculate insert size. Please see
    QualiMap's website for more information about BAM QC:
    http://qualimap.conesalab.org/
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Report containing post-aligment quality-control metrics
    """
    input:
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
    output: 
        txt  = join(workpath,"QC","{name}","genome_results.txt"),
        html = join(workpath,"QC","{name}","qualimapReport.html")
    params:
        outdir = join(workpath,"QC","{name}"),
        rname  = "qualibam"
    message: "Running QualiMap BAM QC with {threads} threads on '{input.bam}' input file"
    threads: int(allocated("threads", "qualimap", cluster))
    envmodules: config['tools']['qualimap']
    shell: """
    unset DISPLAY
    qualimap bamqc -bam {input.bam} \\
        --java-mem-size=92G \\
        -c -ip --gd HUMAN \\
        -outdir {params.outdir} \\
        -outformat HTML \\
        -nt {threads} \\
        --skip-duplicated \\
        -nw 500 \\
        -p NON-STRAND-SPECIFIC
    """


rule samtools_flagstats:
    """
    Quality-control step to assess alignment quality. Flagstat provides 
    counts for each of 13 categories based primarily on bit flags in the 
    FLAG field. Information on the meaning of the flags is given in the 
    SAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Text file containing alignment statistics
    """
    input:
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
    output:
        txt = join(workpath,"QC","{name}.samtools_flagstat.txt")
    params: 
        rname = "flagstat"
    message: "Running SAMtools flagstat on '{input.bam}' input file"
    threads: int(allocated("threads", "samtools_flagstats", cluster))
    envmodules: config['tools']['samtools']
    shell: """
    samtools flagstat --threads {threads} \\
        {input.bam} \\
    > {output.txt}
    """


# Post-variant calling QC-related rules
rule bcftools_stats:
    """
    Quality-control step to collect summary statistics from bcftools stats.
    When bcftools stats is run with one VCF file then stats by non-reference
    allele frequency, depth distribution, stats by quality and per-sample 
    counts, singleton statsistics are calculated. Please see bcftools' 
    documentation for more information: 
    http://samtools.github.io/bcftools/bcftools.html#stats
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Text file containing a collection of summary statistics
    """
    input: 
        vcf = join(workpath, "deepvariant", "VCFs", "{name}.germline.vcf.gz"),
    output: 
        txt = join(workpath, "QC", "BCFStats", "{name}.germline.bcftools_stats.txt"),
    params: 
        rname = "bcfstats",
    message: "Running BCFtools stats on '{input.vcf}' input file"
    threads: int(allocated("threads", "bcftools_stats", cluster))
    envmodules: config['tools']['bcftools']
    shell: """
    bcftools stats \\
        {input.vcf} \\
    > {output.txt}
    """


rule gatk_varianteval:
    """
    Quality-control step to calculate various quality control metrics from a 
    variant callset. These metrics include the number of raw or filtered SNP 
    counts; ratio of transition mutations to transversions; concordance of a
    particular sample's calls to a genotyping chip; number of s per sample.
    Please see GATK's documentation for more information: 
    https://gatk.broadinstitute.org/hc/en-us/articles/360040507171-VariantEval
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    """
    input: 
        vcf = join(workpath, "deepvariant", "VCFs", "{name}.germline.vcf.gz"),
    output: 
        grp = join(workpath, "QC", "VariantEval", "{name}.germline.eval.grp"),
    params:
        rname    = "vareval",
        genome   = config['references']['GENOME'],
        dbsnp    = config['references']['DBSNP'],
        ver_gatk = config['tools']['gatk4'],
        memory   = allocated("mem", "gatk_varianteval", cluster).lower().rstrip('g')
    message: "Running GATK4 VariantEval on '{input.vcf}' input file"
    threads: int(allocated("threads", "gatk_varianteval", cluster))
    envmodules: config['tools']['gatk4']
    shell: """
    gatk --java-options '-Xmx{params.memory}g -XX:ParallelGCThreads={threads}' VariantEval \\
        -R {params.genome} \\
        -O {output.grp} \\
        --dbsnp {params.dbsnp} \\
        --eval {input.vcf} 
    """


rule snpeff:
    """
    Data processing and quality-control step to annotate variants, predict its
    functional effects, and collect various summary statistics about variants and
    their annotations. Please see SnpEff's documentation for more information: 
    https://pcingola.github.io/SnpEff/
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    """
    input:  
        vcf = join(workpath, "deepvariant", "VCFs", "{name}.germline.vcf.gz"),
    output: 
        vcf  = join(workpath, "QC", "SNPeff", "{name}.germline.snpeff.ann.vcf"),
        csv  = join(workpath, "QC", "SNPeff", "{name}.germline.snpeff.ann.csv"),
        html = join(workpath, "QC", "SNPeff", "{name}.germline.snpeff.ann.html"),
    params: 
        rname  = "snpeff",
        genome = config['references']['SNPEFF_GENOME'],
        config = config['references']['SNPEFF_CONFIG'],
        bundle = config['references']['SNPEFF_BUNDLE'],
        memory = allocated("mem", "snpeff", cluster).lower().rstrip('g')
    threads: int(allocated("threads", "snpeff", cluster))
    envmodules: config['tools']['snpeff']
    shell: """
    java -Xmx{params.memory}g -jar ${{SNPEFF_JAR}} \\
        -v -canon -c {params.config} \\
        -csvstats {output.csv} \\
        -stats {output.html} \\
        {params.genome} \\
        {input.vcf} > {output.vcf}
    """


rule vcftools:
    """
    Quality-control step to calculates a measure of heterozygosity on 
    a per-individual basis. The inbreeding coefficient, F, is estimated
    for each individual using a method of moments. Please see VCFtools
    documentation for more information: 
    https://vcftools.github.io/man_latest.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a measure of heterozygosity
    """
    input: 
        vcf = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
    output: 
        het = join(workpath, "QC", "{batch}_variants.het"),
    params: 
        prefix = join(workpath, "QC", "{batch}_variants"),
        rname  = "vcftools",
    message: "Running VCFtools on '{input.vcf}' input file"
    threads: int(allocated("threads", "vcftools", cluster))
    envmodules: config['tools']['vcftools']
    shell: """
    vcftools \\
        --gzvcf {input.vcf} \\
        --het \\
        --out {params.prefix}
    """


rule collectvariantcallmetrics:
    """
    Quality-control step to collect summary metrics about snps and indels
    called in a multisample VCF file. Please see the Broad's documentation
    for more information about each field in the generated log file:
    https://broadinstitute.github.io/picard/picard-metric-definitions.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a collection of metrics relating to snps and indels 
    """
    input: 
        vcf = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
    output: 
        metrics = join(workpath, "QC", "{batch}_variants.variant_calling_detail_metrics"),
    params:
        rname  = "varcallmetrics",
        dbsnp  = config['references']['DBSNP'],
        prefix = join(workpath, "QC", "{batch}_variants"),
        memory = allocated("mem", "collectvariantcallmetrics", cluster).lower().rstrip('g')
    message: "Running Picard CollectVariantCallingMetrics on '{input.vcf}' input file"
    threads: int(allocated("threads", "collectvariantcallmetrics", cluster))
    envmodules: config['tools']['picard']
    shell: """
    java -Xmx{params.memory}g -jar ${{PICARDJARPATH}}/picard.jar \\
        CollectVariantCallingMetrics \\
        INPUT={input.vcf} \\
        OUTPUT={params.prefix} \\
        DBSNP={params.dbsnp} \\
        Validation_Stringency=SILENT
    """


rule somalier:
    """
    To estimate ancestry, Somalier first extracts known sites from mapped 
    reads. Somalier estimates relatedness using the extracted site information 
    to compare across all samples. This step also runs the ancestry estimation
    function in Somalier.
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Exracted sites in (binary) somalier format
    """
    input:
        vcf = join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
    output: 
        somalier = expand(join(workpath, "somalier", "{name}.somalier"), name=samples),
        related  =  join(workpath, "somalier", "relatedness.samples.tsv"),
        ancestry = join(workpath, "somalier", "ancestry.somalier-ancestry.tsv"),
    params:
        rname  = 'somalier',
        outdir = join(workpath, "somalier"),
        # Pedigree file from patient database,
        # pulled locally from API call
        ped    = join(workpath, "BSI_sex_ped.ped"),
        genome = config['references']['GENOME'],
        sites  = config['references']['SOMALIER']['SITES_VCF'],
        ancestry_db = config['references']['SOMALIER']['ANCESTRY_DB'],
    threads: int(allocated("threads", "somalier", cluster))
    container: config['images']['base']
    shell: """ 
    echo "Extracting sites to estimate ancestry."
    somalier extract \\
        -d {params.outdir} \\
        --sites {params.sites} \\
        -f {params.genome} \\
        {input.vcf}
    
    # Check if pedigree file exists,
    # pulled from patient database
    pedigree_option=""
    if [ -f "{params.ped}" ]; then
        # Use PED with relate command
        pedigree_option="-p {params.ped}"
    fi
    echo "Estimating relatedness with pedigree $pedigree_option"
    somalier relate "$pedigree_option" \\
        -i -o {params.outdir}/relatedness \\
        {output.somalier}
    
    echo "Estimating ancestry."
    somalier ancestry \\
        --n-pcs=10 \\
        -o {params.outdir}/ancestry \\
        --labels {params.ancestry_db}/ancestry-labels-1kg.tsv \\
        {params.ancestry_db}/*.somalier ++ \\
        {output.somalier}
    """


rule multiqc:
    """
    Reporting step to aggregate sample summary statistics and quality-control
    information across all samples. This will be one of the last steps of the 
    pipeline. The inputs listed here are to ensure that this step runs last. 
    During runtime, MultiQC will recurively crawl through the working directory
    and parse files that it supports.
    @Input:
        Lists of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report
    """
    input:
        # FastQ Information, flowcell and lanes
        expand(join(workpath,"rawQC","{name}.fastq.info.txt"), name=ignore(samples, skip_qc)),
        # FastQC (before and after trimming)
        expand(join(workpath,"rawQC","{name}.R1_fastqc.zip"), name=ignore(samples, skip_qc)),
        expand(join(workpath,"QC","{name}.sorted_fastqc.html"), name=ignore(samples, skip_qc)),
        # fastp, remove adapter sequences
        expand(join(workpath,"fastqs","{name}.R1.trimmed.fastq.gz"), name=samples),
        # FastQ Screen, contamination screening
        expand(join(workpath,"FQscreen","{name}.R1.trimmed_screen.txt"), name=ignore(samples, skip_qc)),
        # bwa-mem2, align to reference genome
        expand(join(workpath, "BAM", "{name}.sorted.bam.bai"), name=samples),
        # QualiMap2, bam quality control
        expand(join(workpath,"QC","{name}","qualimapReport.html"), name=ignore(samples, skip_qc)),
        # Samtools Flagstat, bam quality control
        expand(join(workpath,"QC","{name}.samtools_flagstat.txt"), name=ignore(samples, skip_qc)),
        # Deepvariant, call germline variants
        expand(join(workpath, "deepvariant", "VCFs", "{name}.vcf.gz"), name=samples),
        # GLnexus, jointly-called norm multi-sample VCF file
        join(workpath, "deepvariant", "VCFs", "joint.glnexus.norm.vcf.gz"),
        # GATK4 SelectVariants, jointly-called norm per-sample VCF file
        expand(join(workpath, "deepvariant", "VCFs", "{name}.germline.vcf.gz"), name=samples),
        # BCFtools Stats, variant quality control
        expand(join(workpath, "QC", "BCFStats", "{name}.germline.bcftools_stats.txt"), name=ignore(samples, skip_qc)),
        # GATK4 VariantEval, variant quality control
        expand(join(workpath, "QC", "VariantEval", "{name}.germline.eval.grp"), name=ignore(samples, skip_qc)),
        # SNPeff, variant quality control and annotation
        expand(join(workpath, "QC", "SNPeff", "{name}.germline.snpeff.ann.vcf"), name=samples),
        # VCFtools, variant quality control 
        expand(join(workpath, "QC", "{batch}_variants.het"), batch=ignore([batch_id], skip_qc)),
        # Picards CollectVariantCallingMetrics, variant quality control
        expand(join(workpath, "QC", "{batch}_variants.variant_calling_detail_metrics"), batch=ignore([batch_id], skip_qc))
    output: 
        report  = join(workpath, "QC", "MultiQC_Report.html"),
    params: 
        rname  = "multiqc",
        workdir = workpath,
    threads: int(allocated("threads", "multiqc", cluster))
    envmodules: config['tools']['multiqc']
    shell: """
    multiqc --ignore '*/.singularity/*' \\
        --ignore 'slurmfiles/' \\
        --exclude peddy \\
        -f --interactive \\
        -n {output.report} \\
        {params.workdir}
    """