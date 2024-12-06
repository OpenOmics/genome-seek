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

    # Enforcing BED is in BED6 format,  
    # Padding features +/- N base pairs  
    # according to strand and merge any
    # overlapping features after padding.
    # Tools like deepvariant/deepsomatic
    # strictly if check the BED input BED
    # BED file is valid. The score must 
    # be a numeric value and strand info
    # must also be valid: "+", "-", or "."
    # Set strand as "." as its not important
    # to retain and only causes problems
    # later with bedtools merge distinct.
    # If originial strand information is
    # kept, it can lead to strand equaling 
    # "+,-" after the bedtools merge. This
    # causes problems with DV/DS.
    awk -F '\\t' -v OFS='\\t' \\
        'function trim(s) {{ sub(/^[ ]+/, "", s); sub(/[ ]+$/, "", s); return s }}; \\
        function isnum(x) {{ return(x==x+0) }} {{ if (isnum($1)) {{print $1}} }} ; \\
        NF>="3" {{print \\
            trim($1), \\
            trim($2), \\
            trim($3), \\
            ((trim($4)=="" || trim($4)==".") ? trim($1)"_"trim($2)"_"trim($3)"_"NR : trim($4)), \\
            (!(isnum(trim($5))) ? 0 : trim($5) ), \\
            "." \\
    }}' {input.bed} \\
    | bedtools slop \\
        -s \\
        -l {params.padding} \\
        -r {params.padding} \\
        -g {output.sizes} \\
        -i - \\
    | bedtools sort -i - \\
    | bedtools merge \\
        -i - \\
        -c 4,5,6 \\
        -o distinct,mean,distinct \\
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
    output:
        bed = join(workpath, "references", "genome.cnvkit_access.bed"),
    params:
        rname   = "cnvkit_access",
        fa      = config['references']['GENOME'],
    threads: int(allocated("threads", "cnvkit_access", cluster)),
    container: config['images']['genome-seek_cnv'],
    envmodules: 
        config['tools']['cnvkit'],
    shell: """
    # Run CNVkit access to calculate
    # sequence-accessible coordinates
    cnvkit.py access \\
        -o {output.bed} \\
        {params.fa}
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


rule sequenza_utils:
    """
    Sequenza can be used to estimate and quantify of purity/ploidy 
    and copy number alteration in sequencing experiments. It contains
    a set of tools to analyze genomic sequencing data from paired 
    normal-tumor samples, including cellularity and ploidy estimation; 
    mutation and copy number (allele-specific and total copy number) 
    detection, quantification and visualization. For more information,
    please visit: https://bitbucket.org/sequenzatools/sequenza/src/master/
    Sequenza will only run if a tumor-normal pair is provided.
    @Input:
        Realigned, recalibrated BAM file (scatter-per-tumor-sample),
        WES capture kit BED file
    @Output:
        1kb Binned Sequenza SEQZ file
    """
    input:
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        bins = join(workpath, "sequenza_out", "{name}", "{name}.1kb.seqz.gz"),
    params:
        rname   = "sequenza_utils",
        outdir  = workpath,
        tumor   = "{name}",
        fasta   = config['references']['GENOME'],
        species = config['references']['SEQUENZA_SPECIES'],
        gc      = config['references']['SEQUENZA_GC'],
        shell_script  = join("workflow", "scripts", "run_sequenza.sh"),
        # Building optional argument for paired normal,
        # if normal not provided, sequenza will NOT run
        normal = lambda w: "{0}".format(
           tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
        # Building optional flag for wes capture kit
        wes_flag = lambda _: "-b {0}".format(wes_bed_file) if run_wes else "",
    threads: int(allocated("threads", "sequenza_utils", cluster)),
    container: config['images']['genome-seek_cnv'],
    envmodules: 
        config['tools']['sequenza'],
    shell: """
    # Create output directory
    mkdir -p sequenza_out/{params.tumor}
    # Runs sequenza-utils bam2seqz to convert BAM 
    # to seqz file format AND runs sequenza-utils 
    # seqz_binning to bin seqz file into 1kb bins
    {params.shell_script} \\
        -s {params.tumor} \\
        -t {input.tumor} \\
        -n {input.normal} \\
        -r {params.fasta} \\
        -c {threads} \\
        -g {params.gc} \\
        -e {params.species} {params.wes_flag}
    """


rule sequenza:
    """
    Sequenza can be used to estimate and quantify of purity/ploidy 
    and copy number alteration in sequencing experiments. It contains
    a set of tools to analyze genomic sequencing data from paired 
    normal-tumor samples, including cellularity and ploidy estimation; 
    mutation and copy number (allele-specific and total copy number) 
    detection, quantification and visualization. For more information,
    please visit: https://bitbucket.org/sequenzatools/sequenza/src/master/
    Sequenza will only run if a tumor-normal pair is provided.
    @Input:
        Realigned, recalibrated BAM file (scatter-per-tumor-sample),
        WES capture kit BED file
    """
    input:
        bins      = join(workpath, "sequenza_out", "{name}", "{name}.1kb.seqz.gz"),
    output:
        solutions = join(workpath, "sequenza_out", "{name}_alternative_solutions.txt"),
    params:
        rname   = "sequenza",
        outdir  = workpath,
        tumor   = "{name}",
        rlang_script  = join("workflow", "scripts", "run_sequenza.R"),
        # Building optional argument for paired normal,
        # if normal not provided, sequenza will NOT run
        normal = lambda w: "{0}".format(
           tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: int(allocated("threads", "sequenza", cluster)),
    container: config['images']['sequenza'],
    envmodules: 
        config['tools']['sequenza'],
    shell: """
    # Runs sequenza to estimate purity/ploidy and CNVs
    Rscript {params.rlang_script} \\
        {input.bins} \\
        {params.outdir}/sequenza_out/{params.tumor} \\
        {threads} \\
        {params.normal}+{params.tumor}
    mv \\
        {params.outdir}/sequenza_out/{params.tumor}/{params.normal}+{params.tumor}_alternative_solutions.txt \\
        {output.solutions}
    """
