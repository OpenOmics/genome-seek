# Functions and rules related to whole-exome sequencing
from scripts.common import (
    abstract_location, 
    allocated
)


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



