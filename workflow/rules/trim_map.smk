# Functions and rules for processing data
from scripts.common import abstract_location, allocated

# Data processing rules
rule fastp:
    """
    Data-processing step to remove adapter sequences and perform quality trimming
    prior to alignment the reference genome.  Adapters are composed of synthetic
    sequences and should be removed prior to alignment. Fastp is much faster than
    trimmomatic and cutadapt, and it will also auto-detect adapter sequences.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Trimmed FastQ file
    """
    input:
        r1=join(workpath,"{name}.R1.fastq.gz"),
        r2=join(workpath,"{name}.R2.fastq.gz"),
    output:
        r1=join(workpath,"fastqs","{name}.R1.trimmed.fastq.gz"),
        r2=join(workpath,"fastqs","{name}.R2.trimmed.fastq.gz"),
        json=join(workpath,"fastqs","{name}.fastp.json"),
        html=join(workpath,"fastqs","{name}.fastp.html")
    params:
        rname='trim',
    threads: int(allocated("threads", "fastp", cluster))
    envmodules: config['tools']['fastp']
    shell: """
    fastp -w {threads} \\
        --detect_adapter_for_pe \\
        --in1 {input.r1} \\
        --in2 {input.r2} \\
        --out1 {output.r1} \\
        --out2 {output.r2} \\
        --json {output.json} \\
        --html {output.html}
    """


rule bwa_mem2:
    """
    Data processing rule to align trimmed reads to the reference genome using 
    bwa-mem2 aligner. The latest version of the aligner is faster than the previous
    version and the size of the resulting reference genome index is much smaller.
    The resulting bam file will have duplicates marked and will be sorted. A BAM
    index will also be created.
    @Input:
        Trimmed FastQ file (scatter)
    @Output:
        Aligned reads in BAM format
    """
    input:
        r1=join(workpath,"fastqs","{name}.R1.trimmed.fastq.gz"),
        r2=join(workpath,"fastqs","{name}.R2.trimmed.fastq.gz"),
    output: 
        bam = temp(join(workpath, "BAM", "{name}.sorted.bam")),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai")
    params:
        rname = "bwamem2",
        genome = config['references']['GENOME'],
        sample = "{name}"
    threads: 
        int(allocated("threads", "bwa_mem2", cluster))
    envmodules: 
        config['tools']['samtools'],
        config['tools']['bwa_mem2'],

    shell: """
    bwa-mem2 mem \\
        -t {threads} \\
        -K 100000000 \\
        -M \\
        -R \'@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:illumina\\tLB:{params.sample}\\tPU:{params.sample}\\tCN:ncbr\\tDS:wgs\' \\
        {params.genome} \\
        {input.r1} \\
        {input.r2} \\
    | samblaster -M \\
    | samtools sort -@{threads} \\
        --write-index \\
        -m 10G - \\
        -o {output.bam}
    """
