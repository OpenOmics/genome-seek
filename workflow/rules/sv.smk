# Functions and rules for calling structural variants
from scripts.common import abstract_location, allocated


rule manta:
    """
    Data processing step to call structural variants (SV) using manta. 
    Manta is optimized for optimized for analysis of germline variation 
    in small sets of individuals and somatic variation in tumor/normal 
    sample pairs. Manta  will discover, assemble, and score large-scale 
    SVs, medium-sized indels and large insertions.  Manta combines paired 
    and split-read evidence during SV discovery and scoring to improve 
    accuracy, but does not require split-reads or successful breakpoint 
    assemblies to report a variant in cases where there is strong evidence 
    otherwise. It provides scoring models for germline variants in small 
    sets of diploid samples and somatic variants in matched tumor/normal 
    sample pairs. Manta's Github Repo: https://github.com/Illumina/manta
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Single-sample VCF file with called structural variants
    """
    input: 
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai"),
    output:
        vcf  = join(workpath, "MANTA", "{name}", "results", "variants", "diploidSV.vcf.gz"),
    params: 
        rname  = "manta",
        outdir = join(workpath, "MANTA", "{name}"),
        workflow = join(workpath, "MANTA", "{name}", "runWorkflow.py"),
        regions  = config['references']['MANTA_CALLREGIONS'],
        genome   = config['references']['GENOME'],
        memory   = allocated("mem", "manta", cluster).rstrip('G'),
    message: "Running manta on '{input.bam}' input file"
    threads: int(allocated("threads", "manta", cluster))
    envmodules: config['tools']['manta']
    shell: """
    # Delete previous attempts output
    # directory to ensure hard restart
    if [ -d "{params.outdir}" ]; then
        rm -rf "{params.outdir}"
    fi

    # Configure Manta SV workflow 
    configManta.py \\
        --callRegions {params.regions} \\
        --bam {input.bam} \\
        --referenceFasta {params.genome} \\
        --runDir {params.outdir}
    
    # Call SV with Manta workflow
    echo "Starting Manta workflow..."
    {params.workflow} \\
        -j {threads} \\
        -g {params.memory} 
    """