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
    container: config['images']['genome-seek_sv']
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
    container: config['images']['genome-seek_sv']
    envmodules: config['tools']['manta']
    shell: """
    # Delete previous attempts output
    # directory to ensure hard restart
    if [ -d "{params.outdir}" ]; then
        rm -rf "{params.outdir}"
    fi

    # Configure Manta somatic SV workflow 
    configManta.py {params.normal_option} \\
        --callRegions {params.regions} \\
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

# Filter Somatic SV calls
rule manta_somatic_filter:
    """
    Data processing step to call somatic structural variants using manta. 
    Manta is optimized for optimized for analysis of germline variation 
    in small sets of individuals and somatic variation in tumor/normal 
    sample pairs. Manta's Github Repo: https://github.com/Illumina/manta
    @Input:
        Single-sample VCF file with called somatic structural variants (scatter)
    @Output:
        Filtered single-sample VCF file with called somatic structural variants
    """
    input: 
        vcf  = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.vcf.gz"),
    output:
        samplevcf  = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.sample.vcf.gz"),
        flagged    = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.flagged.vcf.gz"),
        filtered   = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.filtered.vcf"),
    params: 
        rname       = "manta_somatic_filter",
        sample      = "{name}",
        outdir      = join(workpath, "MANTA", "somatic", "{name}"),
        filtermanta = join(workpath, "workflow", "scripts", "FilterManta.pl"),
        genome      = config['references']['GENOME'],
        filter_ref  = config['references']['MANTA_FILTER_CHROMOSEQ_TRANSLOCATION'],
        memory      = allocated("mem", "manta_somatic", cluster).rstrip('G'),
    threads: int(allocated("threads", "manta_somatic", cluster))
    container: config['images']['genome-seek_sv']
    envmodules: 
        config['tools']['bcftools'],
        config['tools']['samtools'],
        config['tools']['bedtools'],
        config['tools']['svtools'],
        config['tools']['minimap2'],
        config['tools']['perl'],
    shell: """
    # Flag and filter SVs based on the
    # following: read support, contig 
    # re-mapping, and allele fraction. 
    # Also removes non-PASS SVs from
    # somatic SV callset.
    bcftools view \\
        -s {params.sample} \\
        -O z \\
        -o {output.samplevcf} \\
        {input.vcf}
    # Script to filter manta results 
    # from chromoseq pipeline,
    # https://github.com/genome/docker-basespace_chromoseq
    # Needs paths to its depedencies:
    bedtools_path="$(which bedtools)"
    minimap2_path="$(which minimap2)"
    svtools_path="$(which svtools)"
    echo "Paths of all dependencies... ${{bedtools_path}}:${{minimap2_path}}:${{svtools_path}}"
    perl {params.filtermanta} \\
        -r {params.genome} \\
        -k {params.filter_ref} \\
        -b "${{bedtools_path}}" \\
        -p "${{minimap2_path}}" \\
        -s "${{svtools_path}}" \\
        -t {params.outdir} \\
        {output.samplevcf} \\
        {output.flagged}
    # Remove non-PASS-ing SVs
    bcftools view \\
        -O v \\
        -i 'FILTER=="PASS"' \\
        {output.flagged} \\
    > {output.filtered}
    """
