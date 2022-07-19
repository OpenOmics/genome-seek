# Functions and rules for calling structural variants
from scripts.common import (
    abstract_location, 
    allocated
)

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

# Germline SV calling 
rule manta_germline:
    """
    Data processing step to call germline structural variants using manta. 
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
        Single-sample VCF file with called germline structural variants
    """
    input: 
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai"),
    output:
        vcf  = join(workpath, "MANTA", "germline", "{name}", "results", "variants", "diploidSV.vcf.gz"),
    params: 
        rname  = "manta_germline",
        outdir = join(workpath, "MANTA", "germline", "{name}"),
        workflow = join(workpath, "MANTA", "germline", "{name}", "runWorkflow.py"),
        regions  = config['references']['MANTA_CALLREGIONS'],
        genome   = config['references']['GENOME'],
        memory   = allocated("mem", "manta_germline", cluster).rstrip('G'),
    threads: int(allocated("threads", "manta_germline", cluster))
    envmodules: config['tools']['manta']
    shell: """
    # Delete previous attempts output
    # directory to ensure hard restart
    if [ -d "{params.outdir}" ]; then
        rm -rf "{params.outdir}"
    fi

    # Configure Manta germline SV workflow 
    configManta.py \\
        --callRegions {params.regions} \\
        --bam {input.bam} \\
        --referenceFasta {params.genome} \\
        --runDir {params.outdir}
    
    # Call germline SV with Manta workflow
    echo "Starting Manta workflow..."
    {params.workflow} \\
        -j {threads} \\
        -g {params.memory} 
    """


# Somatic SV calling 
rule manta_somatic:
    """
    Data processing step to call somatic structural variants using manta. 
    Manta is optimized for optimized for analysis of germline variation 
    in small sets of individuals and somatic variation in tumor/normal 
    sample pairs. Manta's Github Repo: https://github.com/Illumina/manta
    @Input:
        Realigned and recalibrated BAM file (scatter)
    @Output:
        Single-sample VCF file with called somatic structural variants
    """
    input: 
        tumor = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        vcf  = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.vcf.gz"),
    params: 
        rname  = "manta_somatic",
        outdir = join(workpath, "MANTA", "somatic", "{name}"),
        workflow = join(workpath, "MANTA", "somatic", "{name}", "runWorkflow.py"),
        regions  = config['references']['MANTA_CALLREGIONS'],
        genome   = config['references']['GENOME'],
        memory   = allocated("mem", "manta_somatic", cluster).rstrip('G'),
        # Building optional argument for paired normal
        normal_option = lambda w: "--normalBam {0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
    threads: int(allocated("threads", "manta_somatic", cluster))
    envmodules: config['tools']['manta']
    shell: """
    # Delete previous attempts output
    # directory to ensure hard restart
    if [ -d "{params.outdir}" ]; then
        rm -rf "{params.outdir}"
    fi

    # Configure Manta somatic SV workflow 
    configManta.py {params.normal_option} \\
        --tumorBam {input.tumor} \\
        --referenceFasta {params.genome} \\
        --runDir {params.outdir} \\
        --outputContig
    
    # Call somatic SV with Manta workflow
    echo "Starting Manta workflow..."
    {params.workflow} \\
        -j {threads} \\
        -g {params.memory} 
    """
