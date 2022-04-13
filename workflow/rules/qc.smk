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


rule rawfastqc:
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
    threads: int(allocated("threads", "rawfastqc", cluster))
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