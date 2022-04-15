# Functions and rules for quality-control
from scripts.common import abstract_location, allocated

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
        get_flowcell_lanes = os.path.join("workflow", "scripts", "flowcell_lane.py"),
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
    envmodules: config['tools']['fastq_screen']
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
        -c -ip --gd hg19 \\
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
