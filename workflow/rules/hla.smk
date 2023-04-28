# Functions and rules for HLA typing
from scripts.common import abstract_location, allocated


rule hla:
    """
    Data processing step to infer HLA type. HLA*LA provides HLA 
    typing based on a population reference graph and employs a new 
    linear projection method to align reads to the graph. HLA*LA is 
    faster and less resource-intensive than its predecessor HLA*PRG.
    For more information about HLA*LA, check out this blog post:
    https://genomeinformatics.github.io/HLA-PRG-LA/
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Text file containing HLA inference results
    """
    input:
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
    output:
        txt = join(workpath,"HLA", "{name}", "sample", "hla", "R1_bestguess_G.txt")
    params: 
        rname = "hlala",
        outdir = join(workpath, "HLA", "{name}"),
        graph  = config['references']['HLA_LA_GRAPH'],
    message: "Running HLA*LA on '{input.bam}' input file"
    threads: int(allocated("threads", "hla", cluster))
    container: config['images']['genome-seek_hla']
    envmodules: config['tools']['hla_la']
    shell: """
    HLA-LA.pl \\
        --BAM {input.bam} \\
        --graph {params.graph} \\
        --sampleID sample \\
        --maxThreads {threads} \\
        --workingDir {params.outdir}
    """
