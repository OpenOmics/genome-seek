# Functions and rules for processing data
from scripts.common import abstract_location, allocated

# Pre alignment QC-related rules
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
        join(workpath,"{name}.R1.fastq.gz"),
        join(workpath,"{name}.R2.fastq.gz"),
    output:
        join(workpath,"rawQC","{name}.R1_fastqc.zip"),
        join(workpath,"rawQC","{name}.R2_fastqc.zip"),
    params:
        rname='rawfqc',
        outdir=join(workpath,"rawQC"),
    threads: int(allocated("threads", "rawfastqc", cluster))
    envmodules: config['tools']['fastqc']
    shell: """
    fastqc {input} \\
        -t {threads} \\
        -o {params.outdir}
    """