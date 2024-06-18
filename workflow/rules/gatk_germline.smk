# Functions and rules for calling germline variants
# using GATK4's recommended set of best practices.
# These rules are a subset of the total set of rules
# for calling germline variants, and are intended to
# be conditionally used in conjunction with the full
# set of rules germline rules (see "germline.smk")
# via the --gatk-germline option. If this cli option
# is provided then these rules be executed.
from scripts.common import (
    abstract_location, 
    allocated,
    provided
)


rule gatk_germline_haplotypecaller:
    """
    Call germline SNPs and indels via local re-assembly of haplotypes.
    The HaplotypeCaller is capable of calling SNPs and indels simultaneously 
    via local de-novo assembly of haplotypes in an active region. In other 
    words, whenever the program encounters a region showing signs of variation,
    it discards the existing mapping information and completely reassembles the
    reads in that region. This allows the HaplotypeCaller to be more accurate
    when calling regions that are traditionally difficult to call, for example
    when they contain different types of variants close to each other. It also 
    makes the HaplotypeCaller much better at calling indels than position-based
    callers like UnifiedGenotyper.
    @Input:
        Realigned, recalibrated BAM file (scatter-per-sample-per-chrom)
    @Output:
        Single-sample per-chromosome GVCF file with raw, unfiltered SNP and indel calls. 
    """
    input: 
        bam = join(workpath, "BAM", "{name}.recal.bam"),
    output: 
        gvcf = temp(join(workpath, "haplotypecaller", "gVCFs", "chunks", "{name}", "{chrom}.g.vcf.gz")),
        idx  = temp(join(workpath, "haplotypecaller", "gVCFs", "chunks", "{name}", "{chrom}.g.vcf.gz.tbi")),
    params:
        rname  = "haplotype",
        genome = config['references']['GENOME'], 
        sample = "{name}",
        chrom  = "{chrom}",
        snpsites=config['references']['DBSNP'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "gatk_germline_haplotypecaller", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_haplotypecaller", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "gatk_germline_haplotypecaller", cluster).lower().rstrip('g'),
    message: "Running GATK4 HaplotypeCaller on '{input.bam}' input file"
    threads: int(allocated("threads", "gatk_germline_haplotypecaller", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Call germline variants (SNPs and indels)
    # using GATK4's HaplotypeCaller
    gatk --java-options '-Xmx{params.memory}g' HaplotypeCaller \\
        --reference {params.genome} \\
        --input {input.bam} \\
        --use-jdk-inflater \\
        --use-jdk-deflater \\
        --emit-ref-confidence GVCF \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation \\
        --dbsnp {params.snpsites} \\
        --output {output.gvcf} \\
        --max-alternate-alleles 3 \\
        --intervals {params.chrom}
    """


rule gatk_germline_merge_gvcfs_across_chromosomes:
    """
    Data-processing step to merge gVCFs for each chromosome into a single gVCF.
    This rule is intended to be used in conjunction with the HaplotypeCaller to
    merge the scattered germline calling to speed up the germline calling process.
    @Input:
        Single-sample per-chromosome GVCF file with raw, unfiltered SNP and indel calls.
        (gather-per-sample-per-chrom).
    @Output:
        Single-sample GVCF file with raw, unfiltered SNP and indel calls. 
    """
    input: 
        gvcfs = expand(
            join(workpath, "haplotypecaller", "gVCFs", "chunks", "{{name}}", "{chrom}.g.vcf.gz"),
            chrom=chunks
        ),
        idxs = expand(
            join(workpath, "haplotypecaller", "gVCFs", "chunks", "{{name}}", "{chrom}.g.vcf.gz.tbi"),
            chrom=chunks
        ),
    output: 
        gvcf = join(workpath, "haplotypecaller", "gVCFs", "{name}.g.vcf.gz"),
        idx  = join(workpath, "haplotypecaller", "gVCFs", "{name}.g.vcf.gz.tbi"),
        lsl  = join(workpath, "haplotypecaller", "gVCFs", "{name}.gvcfs.list"),
    params:
        rname  = "gatk_gl_merge_chroms",
        sample = "{name}",
        gvcfdir = join(workpath, "haplotypecaller", "gVCFs", "chunks", "{name}"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory  = lambda _: int(
            int(allocated("mem", "gatk_germline_merge_gvcfs_across_chromosomes", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_merge_gvcfs_across_chromosomes", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "gatk_germline_merge_gvcfs_across_chromosomes", cluster).lower().rstrip('g'),
    threads: int(allocated("threads", "gatk_germline_merge_gvcfs_across_chromosomes", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Get a list of per-chromosome gVCFs
    # to merge into a single gVCF. Avoids
    # ARG_MAX issue which will limit max
    # length of a command.
    find {params.gvcfdir}/ \\
        -type f \\
        -iname '*.g.vcf.gz' \\
    > {output.lsl}

    # Merge GATK Germline gVCFs for each
    # chromosome into a single gVCF file.
    gatk --java-options '-Xmx{params.memory}g' MergeVcfs \\
        --OUTPUT {output.gvcf} \\
        --INPUT {output.lsl}
    """


rule gatk_germline_list_gvcfs_across_samples:
    """
    Data-processing step to create a TSV file to map each sample to its gVCF file.
    This tab-delimited sample map is used by the GenomicsDBImport tool to import
    gVCFs into a GenomicsDB database before joint genotyping. The GATK4 Best Practice
    Workflow for SNP and Indel calling uses GenomicsDBImport to merge GVCFs from 
    multiple samples. GenomicsDBImport for the most part offers the same functionality
    as CombineGVCFs and comes from the Intel-Broad Center for Genomics. The datastore 
    transposes sample-centric variant information across genomic loci to make data 
    more accessible to downstream tools.
    @Input:
        Single-sample GVCF file with raw, unfiltered SNP and indel calls. 
        (gather-across-all-samples)
    @Output:
        GenomicsDBImport sample map file (TSV). 
    """
    input: 
        gvcfs = expand( 
            join(workpath, "haplotypecaller", "gVCFs", "{name}.g.vcf.gz"),
            name=samples
        ),
    output:
        lsl  = join(workpath, "haplotypecaller", "gVCFs", "genomicsdbimport_gvcfs.lsl"),
        tsv  = join(workpath, "haplotypecaller", "gVCFs", "genomicsdbimport_gvcfs.tsv"),
    params:
        rname  = "gatk_gl_list_samples",
        gvcfdir = join(workpath, "haplotypecaller", "gVCFs"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory  = lambda _: int(
            int(allocated("mem", "gatk_germline_list_gvcfs_across_samples", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_list_gvcfs_across_samples", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "gatk_germline_list_gvcfs_across_samples", cluster).lower().rstrip('g'),
    threads: int(allocated("threads", "gatk_germline_list_gvcfs_across_samples", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Get a list of gVCFs across samples
    # for joint genotyping. This avoids
    # ARG_MAX issue which will limit max
    # length of a command.
    find {params.gvcfdir}/ \\
        -maxdepth 1 \\
        -type f \\
        -iname '*.g.vcf.gz' \\
    > {output.lsl}

    # Create a TSV file to map each sample
    # to its gVCF file for GenomicsDBImport.
    paste \\
        <(awk -F '/' '{{print $NF}}' {output.lsl} | sed 's/\.g\.vcf\.gz$//g') \\
        {output.lsl} \\
    > {output.tsv}
    """


rule gatk_germline_genomicsdbimport:
    """
    Data-processing step to import single-sample GVCFs into a GenomicsDB before 
    joint genotyping. The GATK4 Best Practice workflow for SNP and Indel calling
    uses GenomicsDBImport to merge GVCFs from multiple samples. GenomicsDBImport
    for the most part offers the same functionality as CombineGVCFs and comes from
    the Intel-Broad Center for Genomics. The datastore transposes sample-centric
    variant information across genomic loci to make data more accessible to 
    downstream tools.
    @Input:
        GenomicsDBImport sample map file (TSV) 
        (scatter-across-all-samples-per-chrom-chunks)
    @Output:
        Chromosome chunk-ed (chr:start-stop) TileDB for joint genotyping. 
    """
    input: 
        tsv   = join(workpath, "haplotypecaller", "gVCFs", "genomicsdbimport_gvcfs.tsv"),
    output:
        flg   = join(workpath, "haplotypecaller", "gVCFs", "genomicsdb", "{region}", "gvcf_to_tiledb_import.done"),
    params:
        rname = "gatk_gl_gdb_import",
        gdb   = join(workpath, "haplotypecaller", "gVCFs", "genomicsdb", "{region}"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory.
        # NOTE: We subtract 4GB from the total memory
        # as the -Xmx value the tool is run with should
        # be less than the total amount of physical 
        # memory available by at least a few GB, as 
        # the native TileDB library requires additional 
        # memory on top of the Java memory. Failure to 
        # leave enough memory for the native code can
        # result in confusing error messages!
        memory  = lambda _: int(
            int(allocated("mem", "gatk_germline_genomicsdbimport", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_genomicsdbimport", cluster)) 
        )-4 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_genomicsdbimport", cluster).lower().rstrip('g')) - 4,
        chunk = "{region}",
    threads: int(allocated("threads", "gatk_germline_genomicsdbimport", cluster)) - 2,
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # The --genomicsdb-workspace-path must point 
    # to a non-existent or empty directory.
    if [ -d "{params.gdb}" ]; then
        rm -rf "{params.gdb}";
        mkdir -p "{params.gdb}";
    fi

    # Import gVCFs into a TileDB GenomicsDB 
    # datastore before joint genotyping.
    gatk --java-options '-Xmx{params.memory}g' GenomicsDBImport \\
        --use-jdk-inflater \\
        --use-jdk-deflater \\
        --batch-size {threads} \\
        --genomicsdb-workspace-path {params.gdb} \\
        --sample-name-map {input.tsv} \\
        --intervals {params.chunk}
    
    # Touch a flag file to indicate that the
    # GenomicsDBImport step has completed.
    touch {output.flg}
    """


rule gatk_germline_genotypegvcfs:
    """
    Data-processing step to perform joint genotyping on one or more samples 
    pre-called with HaplotypeCaller. A final VCF in which all samples have 
    been jointly genotyped for a given region is produced.
    @Input:
        Chromosome chunk-ed (chr:start-stop) TileDB for joint genotyping. 
        (scatter-across-all-samples-per-chrom-chunks)
    @Output:
        Joint genotyped VCF for a given region (chr:start-stop). 
    """
    input: 
        flg   = join(workpath, "haplotypecaller", "gVCFs", "genomicsdb", "{region}", "gvcf_to_tiledb_import.done"),
    output:
        vcf  = join(workpath, "haplotypecaller", "VCFs", "chunks", "raw_variants_{region}.vcf.gz"),
    params:
        rname = "gatk_gl_genotype",
        genome   = config['references']['GENOME'], 
        snpsites = config['references']['DBSNP'],
        dburi    = join("gendb://haplotypecaller", "gVCFs", "genomicsdb", "{region}"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory.
        memory  = lambda _: int(
            int(allocated("mem", "gatk_germline_genotypegvcfs", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_genotypegvcfs", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_genotypegvcfs", cluster).lower().rstrip('g')) - 2,
        chunk = "{region}",
    threads: int(allocated("threads", "gatk_germline_genotypegvcfs", cluster)),
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Perform joint genotyping on one or more samples
    # at a given region (chr:start-stop).
    gatk --java-options '-Xmx{params.memory}g' GenotypeGVCFs \\
        --reference {params.genome} \\
        --use-jdk-inflater \\
        --use-jdk-deflater \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation \\
        --dbsnp {params.snpsites} \\
        --output {output.vcf} \\
        --variant {params.dburi} \\
        --intervals {params.chunk}
    """


rule gatk_germline_create_cohort_vcf_across_regions:
    """
    Data-processing step to merge the jointly genotyped VCFs for each region 
    into a single VCF file containing variants across all samples.
    @Input:
        Joint genotyped VCF for a given region (chr:start-stop).
        (gather-per-cohort-across-regions) ~ singleton.
    @Output:
        Joint genotyped VCF of the entire cohort across all regions. 
    """
    input: 
        vcfs = expand(
            join(workpath, "haplotypecaller", "VCFs", "chunks", "raw_variants_{region}.vcf.gz"),
            region=regions
        ),
    output: 
        vcf = join(workpath, "haplotypecaller", "VCFs", "raw_variants.vcf.gz"),
        lsl = join(workpath, "haplotypecaller", "VCFs", "raw_variants.region.vcfs.list"),
    params:
        rname  = "gatk_gl_merge_cohort",
        vcfdir = join(workpath, "haplotypecaller", "VCFs", "chunks"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory  = lambda _: int(
            int(allocated("mem", "gatk_germline_create_cohort_vcf_across_regions", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_create_cohort_vcf_across_regions", cluster)) 
        )-2 if run_mode == "uge" \
        else allocated("mem", "gatk_germline_create_cohort_vcf_across_regions", cluster).lower().rstrip('g'),
    threads: int(allocated("threads", "gatk_germline_create_cohort_vcf_across_regions", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Get a list of region VCF files to
    # merge into a single VCF. Avoids an
    # ARG_MAX issue which will limit max
    # length of a command.
    find {params.vcfdir}/ \\
        -maxdepth 1 \\
        -type f \\
        -iname '*.vcf.gz' \\
    > {output.lsl}

    # Merge GATK Germline joint genotyped VCFs 
    # across regions into a single VCF file.
    gatk --java-options '-Xmx{params.memory}g' MergeVcfs \\
        --OUTPUT {output.vcf} \\
        --INPUT {output.lsl}
    """
