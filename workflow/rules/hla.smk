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

    NOTE: HLA-LA graph reference files cannot be bundled within the 
    docekr image due to their size (~29GB) AND HLA-LA reference files
    must exist in a location relative to its installation. As so, we 
    must run singularity manually to bind the local HLA-LA graph dir
    to where HLA-LA expects to find it. For more info, see this issue:
    https://github.com/DiltheyLab/HLA-LA/issues/96

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
        rname  = "hlala",
        outdir = join(workpath, "HLA", "{name}"),
        graph  = config['references']['HLA_LA_GRAPH'],
        sif    = config['images']['genome-seek_hla'],
        bind   = "{0},{1}:{2}".format(
            ','.join(sorted(config['bindpaths'])),
            config['references']['HLA_LA_GRAPH'],
            '/opt2/hla-la/1.0.3/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT'
        ),
    message: "Running HLA*LA on '{input.bam}' input file"
    threads: int(allocated("threads", "hla", cluster))
    # Must run singularity manually,
    # cannot decouple graph reference 
    # files from installation path.
    # See issue for more information:
    # https://github.com/DiltheyLab/HLA-LA/issues/96
    # container: config['images']['genome-seek_hla']
    envmodules: config['tools']['hla_la']
    shell: """
    singularity exec -B {params.bind} {params.sif} \\
    HLA-LA.pl \\
        --BAM {input.bam} \\
        --graph PRG_MHC_GRCh38_withIMGT \\
        --sampleID sample \\
        --maxThreads {threads} \\
        --workingDir {params.outdir}
    """
