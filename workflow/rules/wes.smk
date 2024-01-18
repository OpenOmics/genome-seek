# Functions and rules related to whole-exome sequencing
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


rule build_exome_bed:
    """
    Data processing step reformat and build a new exome capture
    kit bed files. This step will create a padded WES BED file
    and a bgzip-ed and tabix-ed BED file.  
    @Input:
        WES capture kit BED file (singleton)
    @Output:
        Chromomsome to size TSV file
        Padded WES capture kit BED file
        bgzip padded WES capture kit BED file
        Tabix index of bgzip padded WES capture kit BED file
    """
    input: 
        bed = wes_bed_file,
    output:
        sizes  = join(workpath, "references", "genome_chrom_sizes.tsv"),
        padded = join(workpath, "references", "wes_regions_50bp_padded.bed"),
        bgzip  = join(workpath, "references", "wes_regions_50bp_padded.bed.gz"),
        tabix  = join(workpath, "references", "wes_regions_50bp_padded.bed.gz.tbi"),
    params: 
        rname   = "wes_bed",
        padding = '50',
        tmpdir  = tmpdir,
        genome  = config['references']['GENOME'],
    threads: int(allocated("threads", "build_exome_bed", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['vcftools']
    shell: """
    # Get sizes of each chrom
    samtools faidx {params.genome} -o - \\
        | cut -f1,2 \\
    > {output.sizes}

    # Padding features +/- N base pairs  
    # according to strand
    bedtools slop \\
        -s \\
        -l {params.padding} \\
        -r {params.padding} \\
        -g {output.sizes} \\
        -i {input.bed} \\
    > {output.padded}

    # Bgzip padded file and create
    # a tabix index 
    bgzip -c {output.padded} > {output.bgzip}
    tabix -f -p bed {output.bgzip}
    """


rule cnvkit_access:
    """
    CNVkit is a Python library and command-line toolkit to infer and 
    visualize Copy Number Variant (CNV) detection from targeted DNA 
    sequencing. This rule calculates the sequence-accessible coordinates
    in chromosomes from the given reference genome, output as a BED file.
    This BED file is used as an optional argument for CNVkit batch cmd.
    For more information, please visit: https://github.com/etal/cnvkit
    @Input:
        Genomic FASTA file of reference genome (singleton)
    @Output:
        BED file of sequence-accessible coordinates
    """
    input:
        fa  = config['references']['GENOME'],
    output:
        bed = join(workpath, "references", "genome.cnvkit_access.bed"),
    params:
        rname   = "cnvkit_access",
    threads: int(allocated("threads", "cnvkit_access", cluster)),
    container: config['images']['genome-seek_cnv'],
    envmodules: 
        config['tools']['cnvkit'],
    shell: """
    # Run CNVkit access to calculate
    # sequence-accessible coordinates
    cnvkit.py access \\
        -o {output.bed} \\
        {input.fa}
    """


# Rule(s) for calling Copy Number Variation (CNV) from WES data
rule cnvkit_batch:
    """
    CNVkit is a Python library and command-line toolkit to infer and 
    visualize Copy Number Variant (CNV) detection from targeted DNA 
    sequencing. It is designed for use with hybrid capture, including
    both whole-exome and custom target panels, and short-read sequencing.
    It can also be used with whole-genome sequencing, although we are 
    not using it for that purpose (better tools are available for WGS).
    CNVkit will only run if a tumor-normal pair is provided. For more 
    information, please visit: https://github.com/etal/cnvkit
    @Input:
        Realigned, recalibrated BAM file (scatter-per-tumor-sample),
        WES capture kit BED file
    """
    input:
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        bed    = join(workpath, "references", "genome.cnvkit_access.bed"),
        normal = get_normal_recal_bam,
    output:
        fit    = join(workpath, "cnvkit", "{name}", "{name}.call.cns"),
    params:
        rname   = "cnvkit_batch",
        outdir  = workpath,
        tumor   = "{name}",
        fasta   = config['references']['GENOME'],
        bed     = wes_bed_file,
        # Building optional argument for paired normal,
        # if normal not provided, cnvkit will NOT run
        normal_option = lambda w: "-n {0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
    threads: int(allocated("threads", "cnvkit_batch", cluster)),
    container: config['images']['genome-seek_cnv'],
    envmodules: 
        config['tools']['cnvkit'],
    shell: """
    # Run CNVkit with a TN pair
    cnvkit.py batch \\
        --scatter {params.normal_option} \\
        -d cnvkit/{params.tumor}/ \\
        --diagram \\
        -t {params.bed} \\
        -f {params.fasta} \\
        -p {threads} \\
        --drop-low-coverage \\
        -g {input.bed} \\
        -y \\
        --short-names {input.tumor}
    """