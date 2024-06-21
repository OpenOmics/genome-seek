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
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_haplotypecaller", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_haplotypecaller", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_haplotypecaller", cluster).lower().rstrip('g')) - 2,
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
        rname   = "gatk_gl_merge_chroms",
        sample  = "{name}",
        gvcfdir = join(workpath, "haplotypecaller", "gVCFs", "chunks", "{name}"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory  = lambda _: int(
            int(allocated("mem", "gatk_germline_merge_gvcfs_across_chromosomes", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_merge_gvcfs_across_chromosomes", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_merge_gvcfs_across_chromosomes", cluster).lower().rstrip('g')) - 2,
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
        rname   = "gatk_gl_list_samples",
        gvcfdir = join(workpath, "haplotypecaller", "gVCFs"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory  = lambda _: int(
            int(allocated("mem", "gatk_germline_list_gvcfs_across_samples", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_list_gvcfs_across_samples", cluster)) 
        )-1 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_list_gvcfs_across_samples", cluster).lower().rstrip('g')) - 1,
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
        memory = lambda _: int(
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
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_create_cohort_vcf_across_regions", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_create_cohort_vcf_across_regions", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_create_cohort_vcf_across_regions", cluster).lower().rstrip('g')) - 2,
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


rule gatk_germline_build_vqsr_snps:
    """
    Data-processing step to filter SNPs by Variant Quality Score 
    Recalibration (VQSR). GATK's variant calling tools are designed to be 
    very lenient in order to achieve a high degree of sensitivity. This is 
    good because it minimizes the chance of missing real variants, but it 
    does mean that we need to filter the raw callset they produce in order 
    to reduce the amount of false positives, which can be quite large. The 
    established way to filter the raw variant callset is to use variant 
    quality score recalibration (VQSR), which uses machine learning to 
    identify annotation profiles of variants that are likely to be real, 
    and assigns a VQSLOD score to each variant that is much more reliable 
    than the QUAL score calculated by the caller. This tool performs the 
    first pass in VQSR. Specifically, it builds the model that will be 
    used in the second step to actually filter variants.
    @Input:
        Joint genotyped VCF of the entire cohort across all regions. 
        (singleton)
    @Output:
        Recalibration table file for ApplyVQSR.
        Tranches file with various recalibration metrics.
        Rscript file to visualize the recalibration plots.
    """
    input: 
        vcf = join(workpath, "haplotypecaller", "VCFs", "raw_variants.vcf.gz"),
    output:
        recal   = join(workpath, "haplotypecaller", "VCFs", "SNP.output.AS.recal"),
        tranch  = join(workpath, "haplotypecaller", "VCFs", "SNP.output.AS.tranches"),
        rscript = join(workpath, "haplotypecaller", "VCFs", "SNP.output.AS.R"),
    params:
        rname  = "gatk_gl_bld_vqsr_snps",
        # Reference files for SNP VQSR
        genome = config['references']['GENOME'],
        dbsnp  = config['references']['DBSNP'],
        onekg  = config['references']['1000GSNP'],
        hapmap = config['references']['HAPMAP'],
        omni   = config['references']['OMNI'], 
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_build_vqsr_snps", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_build_vqsr_snps", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_build_vqsr_snps", cluster).lower().rstrip('g')) - 2,
    threads: int(allocated("threads", "gatk_germline_build_vqsr_snps", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Build a recalibration model to score 
    # variant quality for filtering purposes.
    # This is the first pass in a two-stage 
    # process called Variant Quality Score 
    # Recalibration (VQSR). 
    gatk --java-options '-Xmx{params.memory}g' VariantRecalibrator \\
        --trust-all-polymorphic \\
        -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 \\
        -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 \\
        -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \\
        --max-gaussians 6 \\
        --reference {params.genome} \\
        --use-jdk-inflater --use-jdk-deflater \\
        -V {input.vcf} \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \\
        {params.hapmap} \\
        --resource:omni,known=false,training=true,truth=false,prior=12.0 \\
        {params.omni} \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 \\
        {params.onekg} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \\
        {params.dbsnp} \\
        -an QD -an DP -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
        -mode SNP \\
        -O {output.recal} \\
        --tranches-file {output.tranch} \\
        --rscript-file {output.rscript}
    """


rule gatk_germline_apply_vqsr_snps:
    """
    Data-processing step to apply a score cutoff to filter variants based
    on a recalibration table. This tool performs the second pass in a two-
    stage process called Variant Quality Score Recalibration (VQSR).
    Specifically, it applies filtering to the input variants based on the
    recalibration table produced in the first step by VariantRecalibrator
    and a target sensitivity value, which the tool matches internally to
    a VQSLOD score cutoff based on the model's estimated sensitivity to
    a set of true variants.
    @Input:
        Recalibration table file for ApplyVQSR.
        Tranches file with various recalibration metrics.
        Rscript file to visualize the recalibration plots.
        (gather-per-cohort-scatter-across-regions)
    @Output:
        A recalibrated VCF file in which each variant of the requested
        type is annotated with its VQSLOD and marked as filtered if
        the score is below the desired quality level.
    """
    input:
        vcf    = join(workpath, "haplotypecaller", "VCFs", "raw_variants.vcf.gz"),
        recal  = join(workpath, "haplotypecaller", "VCFs", "SNP.output.AS.recal"),
        tranch = join(workpath, "haplotypecaller", "VCFs", "SNP.output.AS.tranches"),
    output:
        vcf = join(workpath, "haplotypecaller", "VCFs", "snp_recal_chunks", "snps_recal_variants_{region}.vcf.gz"),
    params:
        rname  = "gatk_gl_apl_vqsr_snps",
        # Reference files for SNP ApplyVQSR
        genome = config['references']['GENOME'],
        chunk  = "{region}",
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_apply_vqsr_snps", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_apply_vqsr_snps", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_apply_vqsr_snps", cluster).lower().rstrip('g')) - 2,
    threads: int(allocated("threads", "gatk_germline_apply_vqsr_snps", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Apply a score cutoff to filter variants
    # based on a recalibration table. This is
    # the last step in VQSR to filter out low
    # quality variants.
    gatk --java-options '-Xmx{params.memory}g' ApplyVQSR \\
        --intervals {params.chunk} \\
        --create-output-variant-index true \\
        --reference {params.genome} \\
        --use-jdk-inflater --use-jdk-deflater \\
        -V {input.vcf} \\
        -mode SNP \\
        --recal-file {input.recal} \\
        --tranches-file {input.tranch} \\
        --truth-sensitivity-filter-level 99.7 \\
        -O {output.vcf}
    """


rule gatk_germline_build_vqsr_indels:
    """
    Data-processing step to filter INDELs by Variant Quality Score 
    Recalibration (VQSR). GATK's variant calling tools are designed to be 
    very lenient in order to achieve a high degree of sensitivity. This is 
    good because it minimizes the chance of missing real variants, but it 
    does mean that we need to filter the raw callset they produce in order 
    to reduce the amount of false positives, which can be quite large. The 
    established way to filter the raw variant callset is to use variant 
    quality score recalibration (VQSR), which uses machine learning to 
    identify annotation profiles of variants that are likely to be real, 
    and assigns a VQSLOD score to each variant that is much more reliable 
    than the QUAL score calculated by the caller. This tool performs the 
    first pass in VQSR. Specifically, it builds the model that will be 
    used in the second step to actually filter variants.
    @Input:
        Joint genotyped VCF of the entire cohort across all regions. 
        (singleton)
    @Output:
        Recalibration table file for ApplyVQSR.
        Tranches file with various recalibration metrics.
        Rscript file to visualize the recalibration plots.
    """
    input: 
        vcf = join(workpath, "haplotypecaller", "VCFs", "raw_variants.vcf.gz"),
    output:
        recal   = join(workpath, "haplotypecaller", "VCFs", "INDEL.output.AS.recal"),
        tranch  = join(workpath, "haplotypecaller", "VCFs", "INDEL.output.AS.tranches"),
        rscript = join(workpath, "haplotypecaller", "VCFs", "INDEL.output.AS.R"),
    params:
        rname  = "gatk_gl_bld_vqsr_indels",
        # Reference files for INDEL VQSR
        genome = config['references']['GENOME'],
        dbsnp  = config['references']['DBSNP'],
        mills  = config['references']['MILLS'],
        axiom  = config['references']['AXIOM'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_build_vqsr_indels", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_build_vqsr_indels", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_build_vqsr_indels", cluster).lower().rstrip('g')) - 2,
    threads: int(allocated("threads", "gatk_germline_build_vqsr_indels", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Build a recalibration model to score 
    # variant quality for filtering purposes.
    # This is the first pass in a two-stage 
    # process called Variant Quality Score 
    # Recalibration (VQSR). 
    gatk --java-options '-Xmx{params.memory}g' VariantRecalibrator \\
        --trust-all-polymorphic \\
        -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 \\
        -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \\
        -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 \\
        -tranche 91.0 -tranche 90.0 \\
        --reference {params.genome} \\
        --use-jdk-inflater --use-jdk-deflater \\
        -V {input.vcf} \\
        --resource:mills,known=false,training=true,truth=true,prior=12.0 \\
        {params.mills} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \\
        {params.dbsnp} \\
        --resource:axiomPoly,known=false,training=true,truth=false,prior=10 \\
        {params.axiom} \\
        -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\
        -mode INDEL \\
        -O {output.recal} \\
        --tranches-file {output.tranch} \\
        --rscript-file {output.rscript} \\
        --max-gaussians 4
    """


rule gatk_germline_apply_vqsr_indels:
    """
    Data-processing step to apply a score cutoff to filter variants based
    on a recalibration table. This tool performs the second pass in a two-
    stage process called Variant Quality Score Recalibration (VQSR).
    Specifically, it applies filtering to the input variants based on the
    recalibration table produced in the first step by VariantRecalibrator
    and a target sensitivity value, which the tool matches internally to
    a VQSLOD score cutoff based on the model's estimated sensitivity to
    a set of true variants.
    @Input:
        Recalibrated, filtered VCF file from apply VQSR for SNPs.
        Recalibration table file for ApplyVQSR.
        Tranches file with various recalibration metrics.
        (gather-per-cohort-scatter-across-regions)
    @Output:
        A recalibrated VCF file in which each variant of the requested
        type is annotated with its VQSLOD and marked as filtered if
        the score is below the desired quality level containing both
        SNPs and INDELs for a given region.
    """
    input:
        vcf    = join(workpath, "haplotypecaller", "VCFs", "snp_recal_chunks", "snps_recal_variants_{region}.vcf.gz"),
        recal  = join(workpath, "haplotypecaller", "VCFs", "INDEL.output.AS.recal"),
        tranch = join(workpath, "haplotypecaller", "VCFs", "INDEL.output.AS.tranches"),
    output:
        vcf = join(workpath, "haplotypecaller", "VCFs", "snp_indel_recal_chunks", "snps_and_indels_recal_variants_{region}.vcf.gz"),
    params:
        rname  = "gatk_gl_apl_vqsr_indels",
        # Reference files for INDEL ApplyVQSR
        genome = config['references']['GENOME'],
        chunk  = "{region}",
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_apply_vqsr_indels", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_apply_vqsr_indels", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_apply_vqsr_indels", cluster).lower().rstrip('g')) - 2,
    threads: int(allocated("threads", "gatk_germline_apply_vqsr_indels", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Apply a score cutoff to filter variants
    # based on a recalibration table. This is
    # the last step in VQSR to filter out low
    # quality variants. The output VCF file
    # contains both SNPs and INDELs.
    gatk --java-options '-Xmx{params.memory}g' ApplyVQSR \\
        --intervals {params.chunk} \\
        --reference {params.genome} \\
        --use-jdk-inflater --use-jdk-deflater \\
        -V {input.vcf} \\
        -mode INDEL \\
        --recal-file {input.recal} \\
        --tranches-file {input.tranch} \\
        --truth-sensitivity-filter-level 99.7 \\
        -O {output.vcf}
    """


rule gatk_germline_scatter_genotype_refinement:
    """
    Data-processing step to calculate genotype posterior probabilities
    given known population genotypes. The tool calculates the posterior
    genotype probability for each sample genotype in a given VCF format
    callset. This tool uses use extra information like allele frequencies
    in relevant populations to further refine the genotype assignments.
    @Input:
        A recalibrated VCF file in which each variant of the requested
        type is annotated with its VQSLOD and marked as filtered if
        the score is below the desired quality level containing both
        SNPs and INDELs for a given region.
        (gather-per-cohort-scatter-across-regions)
    @Output:
        Genotype refined VCF with the following information at a given
        region (chr:start-stop): Genotype posteriors added to the FORMAT
        fields ("PP"), genotypes and GQ assigned according to these
        posteriors, per-site genotype priors added to the INFO field
        ("PG").
    """
    input:
        vcf = join(workpath, "haplotypecaller", "VCFs", "snp_indel_recal_chunks", "snps_and_indels_recal_variants_{region}.vcf.gz"),
    output:
        tmp = temp(
            join(workpath, "haplotypecaller", "VCFs", "gtype_temp_chunks", "snps_and_indels_recal_refinement_variants_{region}.vcf.gz")
        ),
        vcf = temp(
            join(workpath, "haplotypecaller", "VCFs", "gtype_fixed_chunks", "snps_and_indels_recal_refinement_variants_{region}.GTfix.vcf.gz")
        ),
    params:
        rname  = "gatk_gl_scatter_gtype_refine",
        # Reference files for GType Refinement
        genome = config['references']['GENOME'],
        onekg  = config['references']['1000G'],
        exac   = config['references']['EXAC'],
        chunk  = "{region}",
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_scatter_genotype_refinement", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_scatter_genotype_refinement", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_scatter_genotype_refinement", cluster).lower().rstrip('g')) - 2,
    threads: int(allocated("threads", "gatk_germline_scatter_genotype_refinement", cluster))
    container: config['images']['genome-seek']
    envmodules:
        config['tools']['gatk4'],
        config['tools']['bcftools']
    shell: """
    # Calculate genotype posterior probabilities
    # in relevant populations to further refine
    # the genotype assignments.
    gatk --java-options '-Xmx{params.memory}g' CalculateGenotypePosteriors \\
        --reference {params.genome} \\
        --use-jdk-inflater --use-jdk-deflater \\
        -V {input.vcf} \\
        -supporting {params.onekg} \\
        -supporting {params.exac} \\
        -O {output.tmp} \\
        -L {params.chunk}

    # Unphase all genotypes and sort by allele
    # frequency, example: (1|0 becomes 0/1)
    bcftools +setGT \\
        {output.tmp} \\
        -O z \\
        -o {output.vcf} \\
        -- -t a -n u

    # Create an tabix index for the VCF file
    tabix -p vcf {output.vcf}
    """


rule gatk_germline_gather_genotype_refinement:
    """
    Data-processing step to merge the chunked (chr:start-stop) genotype
    refinements into a single VCF file.
    @Input:
        Genotype refined VCF across all regions.
        (gather-per-cohort) ~ singleton.
    @Output:
        Genotype refined VCF file.
    """
    input:
        tmps = expand(
            join(workpath, "haplotypecaller", "VCFs", "gtype_temp_chunks", "snps_and_indels_recal_refinement_variants_{region}.vcf.gz"),
            region=regions
        ),
        vcfs = expand(
            join(workpath, "haplotypecaller", "VCFs", "gtype_fixed_chunks", "snps_and_indels_recal_refinement_variants_{region}.GTfix.vcf.gz"),
            region=regions
        ),
    output:
        t_lsl = join(workpath, "haplotypecaller", "VCFs", "gtype_temp_chunks", "snps_and_indels_recal_refinement.region.vcfs.list"),
        t_vcf = join(workpath, "haplotypecaller", "VCFs", "snps_and_indels_recal_refinement_variants.vcf.gz"),
        v_lsl = join(workpath, "haplotypecaller", "VCFs", "gtype_fixed_chunks", "snps_and_indels_recal_refinement.region.GTfix.vcf.list"),
        v_vcf = join(workpath, "haplotypecaller", "VCFs", "snps_and_indels_recal_refinement_variants.GTfix.vcf.gz"),
    params:
        rname  = "gatk_gl_gather_gtype_refine",
        t_dir = join(workpath, "haplotypecaller", "VCFs", "gtype_temp_chunks"),
        v_dir = join(workpath, "haplotypecaller", "VCFs", "gtype_fixed_chunks"),
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory = lambda _: int(
            int(allocated("mem", "gatk_germline_gather_genotype_refinement", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_germline_gather_genotype_refinement", cluster)) 
        )-2 if run_mode == "uge" \
        else int(allocated("mem", "gatk_germline_gather_genotype_refinement", cluster).lower().rstrip('g')) - 2,
    threads: int(allocated("threads", "gatk_germline_gather_genotype_refinement", cluster))
    container: config['images']['genome-seek']
    envmodules:
        config['tools']['gatk4']
    shell: """
    # Get a list of region VCF files to
    # merge into a single VCF. Avoids an
    # ARG_MAX issue which will limit max
    # length of a command.
    find {params.t_dir}/ \\
        -maxdepth 1 \\
        -type f \\
        -iname '*.vcf.gz' \\
    > {output.t_lsl}
    # Merge the chunked (chr:start-stop) genotype
    # refinements into a single VCF file.
    gatk --java-options '-Xmx{params.memory}g' MergeVcfs \\
        --OUTPUT {output.t_vcf} \\
        --INPUT  {output.t_lsl}

    # Get a list of region VCF files to
    # merge into a single VCF. Avoids an
    # ARG_MAX issue which will limit max
    # length of a command.
    find {params.v_dir}/ \\
        -maxdepth 1 \\
        -type f \\
        -iname '*.vcf.gz' \\
    > {output.v_lsl}
    # Merge the chunked (chr:start-stop) genotype
    # refinements into a single VCF file.
    gatk --java-options '-Xmx{params.memory}g' MergeVcfs \\
        --OUTPUT {output.v_vcf} \\
        --INPUT  {output.v_lsl}
    """
