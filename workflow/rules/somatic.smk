# Functions and rules for processing data
from scripts.common import (
    abstract_location, 
    allocated, 
    joint_option
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


def get_normal_sorted_bam(wildcards):
    """
    Returns a tumor samples paired normal
    See config['pairs'] for tumor, normal pairs.
    """
    normal = tumor2normal[wildcards.name]
    if normal:
        # Runs in a tumor, normal mode
        return join(workpath, "BAM", "{0}.sorted.bam".format(normal))
    else:
        # Runs in tumor-only mode
        return []


def get_normal_pileup_table(wildcards):
    """
    Returns a tumor samples paired normal pileup
    See config['pairs'] for tumor, normal pairs.
    """
    tumor = wildcards.name
    normal = tumor2normal[tumor]
    if normal:
        # Runs in a tumor, normal mode
        # Paired normal pileup contains
        # the tumor sample's name in file
        return join(workpath, "mutect2", "somatic", "{0}.normalPileup.table".format(tumor))
    else:
        # Runs in tumor-only mode
        return []


def get_somatic_tn_callers(wildcards):
    """Returns somatic variants found with tumor-normal variant
    callers. For tumor-normal samples, extra somatic callers 
    (i.e. MuSE and Strelka and DeepSomatic) are run. Tumor-only 
    samples return an empty list (rule already has reference 
    in input section). See config['pairs'] for tumor, normal pairs.
    """
    tumor = wildcards.name
    normal = tumor2normal[tumor]
    if normal:
        # Callers = MuSE, Strelka
        return [
            join(workpath, caller, "somatic", "{0}.{1}.filtered.norm.vcf".format(tumor, caller)) \
            for caller in tn_somatic_callers
        ]
    else:
        # No paired normal, return nothing
        return []


# Data processing rules for calling somatic variants
rule octopus_somatic:
    """
    Data-processing step to call somatic variants. Octopus is a bayesian 
    haplotype-based mutation caller. Octopus takes inspiration from particle 
    filtering by constructing a tree of haplotypes and dynamically pruning 
    and extending the tree based on haplotype posterior probabilities in a 
    sequential manner. This rule is scattered across genomic intervals or 
    chunks to reduce its overall runtime. More information about Octopus 
    can be found here: https://github.com/luntergroup/octopus
    @Input:
        Realigned, recalibrated BAM file (scatter-per-sample-per-chunk)
    @Output:
        Somatic variants in VCF format  
    """
    input:
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        vcf = join(workpath, "octopus", "somatic", "chunks", "{region}", "{name}.vcf.gz"),
    params: 
        genome  = config['references']['GENOME'],
        rname   = "octosomatic",
        chunk   = "{region}",
        tumor   = "{name}",
        wd      = workpath,
        tmpdir  = join("octopus", "somatic", "chunks", "{region}", "{name}_tmp"),
        tmppath = join(workpath, "octopus", "somatic", "chunks", "{region}", "{name}_tmp"),
        s_model = config['references']['OCTOPUS_SOMATIC_FOREST_MODEL'],
        g_model = config['references']['OCTOPUS_GERMLINE_FOREST_MODEL'],
        error   = config['references']['OCTOPUS_ERROR_MODEL'],
        # Building optional argument for paired normal 
        normal_option = lambda w: "--normal-sample {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "octopus_somatic", cluster))
    container: config['images']['octopus']
    shell: """
    mkdir -p '{params.tmppath}'
    octopus --threads {threads} \\
        -C cancer \\
        --working-directory {params.wd} \\
        --temp-directory-prefix {params.tmpdir} \\
        -R {params.genome} \\
        -I {input.normal} {input.tumor} {params.normal_option} \\
        -o {output.vcf} \\
        --forest-model {params.g_model} \\
        --somatic-forest-model {params.s_model} \\
        --annotations AC AD DP \\
        -T {params.chunk}
    """


rule octopus_merge:
    """
    Data-processing step to merge scattered variant calls from Octopus. Octopus 
    is scattered across genomic intervals or chunks to reduce its overall runtime. 
    @Input:
        Somatic variants in VCF format (gather-per-sample-per-chunks)
    @Output:
        Per sample somatic variants in VCF format  
    """
    input:
        vcfs = expand(join(workpath, "octopus", "somatic", "chunks", "{region}", "{{name}}.vcf.gz"), region=regions),
        bed = provided(join(workpath, "references", "wes_regions_50bp_padded.bed"), run_wes),
    output:
        lsl  = join(workpath, "octopus", "somatic", "{name}.list"),
        raw  = join(workpath, "octopus", "somatic", "{name}.octopus.raw.vcf"),
        vcf  = join(workpath, "octopus", "somatic", "{name}.octopus.unfiltered.vcf"),
        chroms  = temp(join(workpath, "octopus", "somatic", "{name}.chroms")),
        starts  = temp(join(workpath, "octopus", "somatic", "{name}.starts")),
        stops   = temp(join(workpath, "octopus", "somatic", "{name}.stops")),
        bed     = temp(join(workpath, "octopus", "somatic", "{name}.list.bed")),
        sortbed = temp(join(workpath, "octopus", "somatic", "{name}.sort.bed")),
        sortlsl = join(workpath, "octopus", "somatic", "{name}.sort.list"),
    params: 
        genome = config['references']['GENOME'],
        rname  = "octomerge",
        tumor  = "{name}",
        octopath = join(workpath, "octopus", "somatic", "chunks"),
        # Building option for WES, if WES use padded 
        # WES BED file as regions file
        wes_regions_option = lambda _: "--regions-file {0}".format(
            join(workpath, "references", "wes_regions_50bp_padded.bed"),
        ) if run_wes else '',
    threads: 
        int(allocated("threads", "octopus_merge", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['bcftools'],
    shell: """
    # Create list of chunks to merge
    find {params.octopath} -iname '{params.tumor}.vcf.gz' \\
        > {output.lsl}
    
    # Sort list of chunks to merge
    awk -F '/' '{{print $(NF-1)}}' {output.lsl} \\
        | awk -F ':' '{{print $(1)}}' > {output.chroms}
    awk -F '/' '{{print $(NF-1)}}' {output.lsl} \\
        | awk -F ':' '{{print $(2)}}' \\
        | awk -F '-' '{{print $(NF-1)}}' > {output.starts}
    awk -F '/' '{{print $(NF-1)}}' {output.lsl} \\
        | awk -F ':' '{{print $(2)}}' \\
        | awk -F '-' '{{print $(NF)}}' > {output.stops}
    paste {output.chroms} \\
        {output.starts} \\
        {output.stops} \\
        {output.lsl} \\
        > {output.bed}
    bedtools sort \\
        -i {output.bed} \\
        -faidx {params.genome}.fai \\
        > {output.sortbed}
    cut -f4 {output.sortbed} > {output.sortlsl}

    # Merge octopus chunk calls,
    # contains both germline and
    # somatic variants
    bcftools concat \\
        --threads {threads} \\
        -d exact \\
        -a \\
        -f {output.sortlsl} \\
        -o {output.raw} \\
        -O v {params.wes_regions_option}
    # Filter Octopus callset for 
    # variants with SOMATIC tag
    grep -E "#|CHROM|SOMATIC" {output.raw} \\
        > {output.vcf}
    """


rule octopus_filter:
    """
    Data-processing step to merge scattered variant calls from Octopus. Octopus 
    is scattered across genomic intervals or chunks to reduce its overall runtime. 
    @Input:
        Somatic variants in VCF format (gather-per-sample-per-chunks)
    @Output:
        Per sample somatic variants in VCF format  
    """
    input:
        vcf  = join(workpath, "octopus", "somatic", "{name}.octopus.unfiltered.vcf"),
    output:
        vcfa    = join(workpath, "octopus", "somatic", "{name}.octopus.PASS.vcf"),
        vcfb    = join(workpath, "octopus", "somatic", "{name}.octopus.SOMATIC.vcf"),
        vcfsort = join(workpath, "octopus", "somatic", "{name}.octopus.vcf"),
    params: 
        genome = config['references']['GENOME'],
        rname  = "octofilter",
        tumor  = "{name}",
        octopath = join(workpath, "octopus", "somatic", "chunks"),
        # Building optional argument for paired normal
        bcftools_filter_i_option = lambda w: 'FMT/FT[1:0]=="PASS"' \
            if tumor2normal[w.name] else 'FMT/FT[0:0]=="PASS"',
    threads: 
        int(allocated("threads", "octopus_filter", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['bcftools']
    shell: """
    bcftools filter \\
        -o {output.vcfa} \\
        -O v \\
        -i 'FILTER=="PASS"' \\
        {input.vcf}
    
    awk -F '\\t' -v OFS='\\t' \\
        '( $1 ~ /^#/ ) || ( $8 ~ /SOMATIC/ ) {{print}}' \\
        {output.vcfa} \\
        > {output.vcfb}

    bcftools sort \\
        -o {output.vcfsort} \\
        -O v \\
        {output.vcfb}
    """


rule octopus_germline:
    """
    @NOTE: This rule is no longer being run in the pipeline. It has been removed in 
    favor of calling germline variants with DeepVariant, see "germline.smk", AND calling 
    germline variants while calling somatic variants by providing the --forest-model 
    option at runtime. 

    Data-processing step to call germline variants in normals. Many somatic downstream 
    tools take germline calls as part of their input. This rule should only be run with 
    normal samples in tumor, normal pairs. This is why it is included in the somatic.smk
    file. Germline calling with Octopus is significantly faster than calling somatic 
    variants, and as so, there is no need to scatter the calls. More information about 
    Octopus canbe found here: https://github.com/luntergroup/octopus
    @Input:
        Realigned, recalibrated BAM file of a normal
    @Output:
        Germline variants in VCF format  
    """
    input:
        normal  = join(workpath, "BAM", "{name}.recal.bam"),
    output:
        vcf  = join(workpath, "octopus", "germline", "{name}.octopus.vcf"),
    params: 
        genome = config['references']['GENOME'],
        rname  = "octogermline",
        normal = "{name}",
        wd = workpath,
        tmpdir = join("octopus", "germline", "{name}_tmp"),
        model = config['references']['OCTOPUS_GERMLINE_FOREST_MODEL'],
        error = config['references']['OCTOPUS_ERROR_MODEL'],
        # Regions to evaluate: chr1 to chrM
        chroms = "{0}".format(" ".join(config['references']['CHR_CHUNKS']))
    threads: 
        int(allocated("threads", "octopus_germline", cluster))
    container: config['images']['octopus']
    shell: """
    octopus --threads {threads} \\
        --working-directory {params.wd} \\
        --temp-directory-prefix {params.tmpdir} \\
        -R {params.genome} \\
        -I {input.normal} \\
        -o {output.vcf} \\
        --forest-model {params.model} \\
        --annotations AC AD AF DP \\
        -T {params.chroms}
    """


rule gatk_scatter_mutect2:
    """Data-processing step to call somatic short mutations via local assembly 
    of haplotypes. Short mutations include single nucleotide (SNA) and insertion 
    and deletion (indel) alterations. The caller uses a Bayesian somatic genotyping 
    model that differs from the original MuTect (which is depreicated). Mutect2 is 
    scatter per chromosome to decrease its overall runtime. More information about 
    Mutect2 can be found on the Broad's website: https://gatk.broadinstitute.org/
    @Input:
        Realigned, recalibrated BAM file (scatter-per-sample-per-chrom)
    @Output:
        Per sample somatic variants in VCF format  
    """
    input: 
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        vcf   = join(workpath, "mutect2", "chrom_split", "{name}.{chrom}.vcf"),
        orien = join(workpath, "mutect2", "chrom_split", "{name}.{chrom}.f1r2.tar.gz"),
        stats = join(workpath, "mutect2", "chrom_split", "{name}.{chrom}.vcf.stats")
    params:
        tumor  = '{name}',
        chrom  = '{chrom}',
        rname  = 'mutect2',
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        germsource = config['references']['GNOMAD'],
        # Building optional argument for paired normal
        normal_option = lambda w: "-normal {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
        i_option = lambda w: "-I {0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "gatk_scatter_mutect2", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    gatk Mutect2 \\
        -R {params.genome} \\
        -I {input.tumor} {params.i_option} {params.normal_option} \\
        --panel-of-normals {params.pon} \\
        --germline-resource {params.germsource} \\
        -L {params.chrom} \\
        -O {output.vcf} \\
        --f1r2-tar-gz {output.orien} \\
        --independent-mates
    """


rule gatk_gather_mutect2:
    """Data-processing step to gather chromosome chunked somatic mutation calls 
    from Mutect2. This step uses CombineVariants from GATK3 since it has not been
    ported to GATK4 and picard MergeVcfs does not merge genotypes (not functionally
    eqiuvalent). Mutect2 is scatter per chromosome to decrease its overall runtime.
    @Input:
        Per sample somatic variants in VCF format (gather-per-sample-per-chrom)
    @Output:
        Gathered sample somatic variants in VCF format  
    """
    input: 
        vcf = expand(join(workpath, "mutect2", "chrom_split", "{{name}}.{chrom}.vcf"), chrom=chunks),
        bed = provided(wes_bed_file, run_wes),
    output:
        vcf = temp(join(workpath, "mutect2", "somatic", "{name}.mutect2.tmp.vcf")),
    params:
        tmpdir = tmpdir,
        rname  = 'merge_mutect2',
        genome = config['references']['GENOME'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "gatk_gather_mutect2", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_gather_mutect2", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "gatk_gather_mutect2", cluster).lower().rstrip('g'),
        # Building optional argument for paired normal
        multi_variant_option = joint_option(
            '--variant',
            expand(join(workpath, "mutect2", "chrom_split", "{{name}}.{chrom}.vcf"), chrom=chunks),
        ),
        # Building option for WES, if WES run
        # build option for intervals file
        wes_intervals_option = lambda _: "--intervals {0}".format(
            wes_bed_file,
        ) if run_wes else '',
    threads: 
        max(int(allocated("threads", "gatk_gather_mutect2", cluster))-1, 1)
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk3']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    java -Xmx{params.memory}g -Djava.io.tmpdir=${{tmp}} \
        -XX:ParallelGCThreads={threads} -jar $GATK_JAR -T CombineVariants \\
        --use_jdk_inflater --use_jdk_deflater \\
        -R {params.genome} \\
        --filteredrecordsmergetype KEEP_UNCONDITIONAL \\
        --assumeIdenticalSamples \\
        -o {output.vcf} {params.wes_intervals_option} \\
        {params.multi_variant_option}
    """


rule gatk_learnReadOrientationModel:
    """Data-processing step to learn read orientation artifacts for Mutect2's bias
    filter. If you suspect any of your samples of substitution errors that occur on
    a single strand before sequencing you should definitely use Mutect2's orientation 
    bias filter. This applies to all FFPE tumor samples and samples sequenced with 
    Illumina Novaseq machines, among others. You can run the filter even when you're 
    not suspicious. It won't hurt accuracy and the CPU cost is now quite small. 
    In this step, we are also gathering the chromosome chunked stats output with 
    MergeMutectStats to save overhead for job submission. FilterMutectCalls will 
    need an aggregated stats file at run time. More information about Mutect2 can
    be found on the Broad's website: https://gatk.broadinstitute.org/
    @Input:
        Mutect2 --f1r2-tar-gz output (gather-per-sample-per-chrom)
        Mutect2 --stats output       (gather-per-sample-per-chrom)
    @Output:
        Per sample stats file for FilterMutectCalls
        Per sample read orientation model for GetPileupSummaries -> FilterMutectCalls
    """
    input: 
        orien = expand(join(workpath, "mutect2", "chrom_split", "{{name}}.{chrom}.f1r2.tar.gz"), chrom=chunks),
        stats = expand(join(workpath, "mutect2", "chrom_split", "{{name}}.{chrom}.vcf.stats"),  chrom=chunks),
    output:
        orien = join(workpath, "mutect2", "somatic", "{name}.read-orientation-model.tar.gz"),
        stats = join(workpath, "mutect2", "somatic", "{name}.vcf.stats"),
    params:
        tumor  = '{name}',
        rname  = 'learnReadOrien',
        genome = config['references']['GENOME'],
        # Building multi arguments for chrom chunks
        multi_stats_option = joint_option(
            '-stats',
            expand(join(workpath, "mutect2", "chrom_split", "{{name}}.{chrom}.vcf.stats"), chrom=chunks)
        ),
        multi_orien_option = joint_option(
            '--input',
            expand(join(workpath, "mutect2", "chrom_split", "{{name}}.{chrom}.f1r2.tar.gz"), chrom=chunks)
        ),
    threads: 
        int(allocated("threads", "gatk_learnReadOrientationModel", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    # Gather Mutect2 stats
    gatk MergeMutectStats \\
        {params.multi_stats_option} \\
        -O {output.stats} &
    # Learn read orientaion model
    # for artifact filtering 
    gatk LearnReadOrientationModel \\
        --output {output.orien} \\
        {params.multi_orien_option} &
    wait
    """


rule gatk_tumorPileup:
    """Data-processing step to summarize read support for a set number of known 
    tumor variant sites in a tumor, normal pair to estimate contamination with 
    GATK4 CalculateContamination. This is a part of a series of steps to call 
    somatic mutations using GATK4 Mutect2. More information about Mutect2 can
    be found on the Broad's website: https://gatk.broadinstitute.org/
    @Input:
        Realigned, recalibrated BAM file for a tumor in a TN pair
    @Output:
        PileupSummaries of tumor for CalculateContamination -> FilterMutectCalls
    """
    input:
        tumor = join(workpath, "BAM", "{name}.recal.bam"),
    output:
        summary = join(workpath, "mutect2", "somatic", "{name}.tumorPileup.table"),
    params:
        tumor  = '{name}',
        rname  = 'tumorPileup',
        genome = config['references']['GENOME'],
        gsnp   = config['references']['1000GSNP'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "gatk_tumorPileup", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_tumorPileup", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "gatk_tumorPileup", cluster).lower().rstrip('g'),
    threads: 
        int(allocated("threads", "gatk_tumorPileup", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    gatk --java-options '-Xmx{params.memory}g' GetPileupSummaries \\
        -I {input.tumor} \\
        -V {params.gsnp} \\
        -L {params.gsnp} \\
        -O {output.summary}
    """


rule gatk_normalPileup:
    """Data-processing step to summarize read support for a set number of known 
    normal variant sites in a tumor, normal pair to estimate contamination with 
    GATK4 CalculateContamination. This is a part of a series of steps to call 
    somatic mutations using GATK4 Mutect2. This rule is only called if a given 
    sample has a paired normal. More information about Mutect2 can be found on 
    the Broad's website: https://gatk.broadinstitute.org/
    @Input:
        Realigned, recalibrated BAM file for a normal in a TN pair
    @Output:
        PileupSummaries of normal for CalculateContamination -> FilterMutectCalls
    """
    input: 
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        summary = join(workpath, "mutect2", "somatic", "{name}.normalPileup.table"),
    params:
        tumor  = '{name}',
        rname  = 'normalPileup',
        genome = config['references']['GENOME'],
        gsnp   = config['references']['1000GSNP'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "gatk_normalPileup", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_normalPileup", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "gatk_normalPileup", cluster).lower().rstrip('g'),
    threads: 
        int(allocated("threads", "gatk_normalPileup", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    gatk --java-options '-Xmx{params.memory}g' GetPileupSummaries \\
        -I {input.normal} \\
        -V {params.gsnp} \\
        -L {params.gsnp} \\
        -O {output.summary}
    """


rule gatk_contamination:
    """Data-processing step to calculate the fraction of reads coming from 
    cross-sample contamination. More information about Mutect2 can be found 
    on the Broad's website: https://gatk.broadinstitute.org/
    @Input:
        PileupSummaries of tumor for CalculateContamination 
        PileupSummaries of normal for CalculateContamination (optional)
    @Output:
        Estimated contamination table for FilterMutectCalls
    """
    input: 
        tumor  = join(workpath, "mutect2", "somatic", "{name}.tumorPileup.table"),
        normal = get_normal_pileup_table,
    output:
        summary = join(workpath, "mutect2", "somatic", "{name}.contamination.table"),
    params:
        tumor  = '{name}',
        rname  = 'calcontaim',
        genome = config['references']['GENOME'],
        # Building optional argument for paired normal
        normal_option = lambda w: "--matched-normal {0}".format(
            join(workpath, "mutect2", "somatic", "{0}.normalPileup.table".format(w.name))
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "gatk_contamination", cluster))
    container: config['images']['genome-seek']
    envmodules: config['tools']['gatk4']
    shell: """
    gatk CalculateContamination \\
        -I {input.tumor} {params.normal_option} \\
        -O {output.summary}
    """


rule gatk_filter_mutect2:
    """Data-processing step to filter somatic SNVs and indels called by Mutect2 
    using the learned read orientation model. More information about Mutect2 can 
    be found on the Broad's website: https://gatk.broadinstitute.org/
    @Input:
        Per sample somatic variants from CombineVariants
        Per sample stats file from MergeMutectStats
        Per sample read orientation from LearnReadOrientationModel
        Estimated contamination table from CalculateContamination 
    @Output:
        Per sample VCF of filtered somatic SNVs and indels
    """
    input:
        vcf     = join(workpath, "mutect2", "somatic", "{name}.mutect2.tmp.vcf"),
        orien   = join(workpath, "mutect2", "somatic", "{name}.read-orientation-model.tar.gz"),
        stats   = join(workpath, "mutect2", "somatic", "{name}.vcf.stats"),
        summary = join(workpath, "mutect2", "somatic", "{name}.contamination.table"),
    output:
        vcf   = join(workpath, "mutect2", "somatic", "{name}.mutect2.vcf"),
    params:
        rname  = 'filtmutect2',
        genome = config['references']['GENOME'],
    threads: 
        int(allocated("threads", "gatk_filter_mutect2", cluster))
    container: config['images']['genome-seek']
    envmodules:
        config['tools']['gatk4'],
        config['tools']['vcftools'],
    shell: """
    # Mutect2 orien bias filter,
    # removing the contamination filter 
    # option due to causing low recall
    # with SEQC2 truth set
    gatk FilterMutectCalls \\
        -R {params.genome} \\
        -V {input.vcf} \\
        --ob-priors {input.orien} \\
        -O {output.vcf} \\
        --stats {input.stats} 
    """


rule hmftools_sage:
    """Data-processing step to call somatic variants in TO and TN samples 
    using hmftools sage. HMF Tools is a suite of tools the Hartwig Medical 
    Foundation developed to analyze genomic data. Sage can be run with WES
    data using the same set of options for WGS. At the current moment, sage
    does not have an option to restrict variant calling to specific regions.
    It does have an -high_depth_mode option; however, the authors state it 
    should only be used for small targeted panels. In the 'somatic_selectvar'
    rule, any variants outside the padded regions/capture-kit BED file are
    removed in WES data. For more information about hmftools visit github:
    https://github.com/hartwigmedical/hmftools
    @Input:
        Sorted BAM file (scatter-per-tumor-sample)
    @Output:
        Per sample somatic variants in VCF format
    """
    input:
        tumor  = join(workpath, "BAM", "{name}.sorted.bam"),
        normal = get_normal_sorted_bam
    output:
        vcf  = join(workpath, "sage", "somatic", "{name}.sage.vcf"),
    params:
        rname     = 'hmfsage',
        tumor     = '{name}',
        genome       = config['references']['GENOME'],
        amber_jar    = config['references']['HMFTOOLS_SAGE_JAR'],
        ref_version  = config['references']['HMFTOOLS_SAGE_REF_VERSION'],
        hotspots     = config['references']['HMFTOOLS_SAGE_HOTSPOTS'],
        panel        = config['references']['HMFTOOLS_SAGE_PANEL'],
        high_conf    = config['references']['HMFTOOLS_SAGE_HIGH_CONF'],
        ensembl_data = config['references']['HMFTOOLS_SAGE_ENSEMBL_DATA'],
        # For UGE/SGE clusters memory is allocated 
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "hmftools_sage", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "hmftools_sage", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "hmftools_sage", cluster).lower().rstrip('g'),
        # Building optional argument for paired normal
        normal_name = lambda w: "-reference {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
        normal_bam = lambda w: "-reference_bam {0}.sorted.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "hmftools_sage", cluster)),
    container: config['images']['genome-seek_cnv']
    envmodules: config['tools']['rlang']
    shell: """
    # Call somatic variants with hmftools
    # Somatic Alterations in Genome (SAGE)
    java -Xmx{params.memory}g -cp {params.amber_jar} \\
        com.hartwig.hmftools.sage.SageApplication \\
            -threads {threads} \\
            -tumor {params.tumor} {params.normal_name} \\
            -tumor_bam {input.tumor}  {params.normal_bam} \\
            -ref_genome_version {params.ref_version} \\
            -ref_genome {params.genome} \\
            -hotspots {params.hotspots} \\
            -panel_bed {params.panel} \\
            -high_confidence_bed {params.high_conf} \\
            -ensembl_data_dir {params.ensembl_data} \\
            -output_vcf {output.vcf}
    """


rule clairs_tumor_only:
    """Data-processing step to call somatic variants in tumor-only samples using 
    ClairS. ClairS is a deep-learning based variant caller that uses an ensembl
    of two neural networks to call somatic variants. ClairS-TO is unique in that
    it can call somatic variants without a matched normal. More information about 
    ClairS-TO can be found here: https://github.com/HKU-BAL/ClairS-TO
    @Input:
        Realigned, recalibrated BAM file (scatter-per-tumor-sample)
    @Output:
        Per sample somatic variants in VCF format
    """
    input:
        tumor = join(workpath, "BAM", "{name}.recal.bam"),
        bed   = provided(join(workpath, "references", "wes_regions_50bp_padded.bed"), run_wes),
    output:
        snps   = join(workpath, "clairs", "somatic", "{name}", "snv.vcf.gz"),
        indels = join(workpath, "clairs", "somatic", "{name}", "indel.vcf.gz"),
        tmp    = join(workpath, "clairs", "somatic", "{name}", "clairs_snps_indels.vcf"),
        vcf    = join(workpath, "clairs", "somatic", "{name}.clairs.vcf"),
    params:
        rname   = 'clairs_to',
        tumor   = '{name}',
        genome  = config['references']['GENOME'],
        outdir = join(workpath, "clairs", "somatic", "{name}"),
        wes_region_option = lambda _: "--bed_fn {0}".format(
            join(workpath, "references", "wes_regions_50bp_padded.bed"),
        ) if run_wes else '',
    threads: 
        # ClairS-TO over utilizes threads,
        # testing has shown it over utilizes
        # around 50% of the threads allocated
        max(int(int(allocated("threads", "clairs_tumor_only", cluster))/2.0), 2),
    container: config['images']['clairs-to']
    envmodules: config['tools']['rlang']
    shell: """
    # Call somatic variants with ClairS-TO,
    # run in isolated sample directory to
    # collisions in file names
    /opt/bin/run_clairs_to \\
        --tumor_bam_fn {input.tumor}  \\
        --ref_fn {params.genome} \\
        --threads {threads} \\
        --platform  ilmn {params.wes_region_option} \\
        --output_dir {params.outdir} \\
        --conda_prefix /opt/micromamba/envs/clairs-to
    
    # Concatenate SNPs and Indels
    bcftools concat \\
        -a \\
        -O v \\
        -o {output.tmp} \\
        {output.snps} \\
        {output.indels} 
    
    # Filter for PASS variants
    bcftools view \\
        -f 'PASS' \\
        -O v \\
        -o {output.vcf} \\
        {output.tmp}
    """


rule deepsomatic_make_examples:
    """
    Data processing step to call somatic variants using deep neural 
    network. The make_examples step prepares the input data for the
    deepsomatic's CNN. DeepSomatic is an extension of deep learning-
    based variant caller Deepvariant. It is composed of multiple steps
    that takes aligned reads (in BAM or CRAM format), produces pileup
    image tensors from them, classifies each tensor using a convolutional 
    neural network, and finally reports the results in a standard VCF or
    gVCF file. This rule is the first step in the deepsomatic pipeline:
     * 1. make_examples        (CPU, parallelizable with gnu-parallel)
       2. call_variants        (GPU, use a GPU node)
       3. postprocess_variants (CPU, single-threaded)
    Running deepsomatic in a single step using run_deepsomatic is not 
    optimal for large-scale projects as it will consume resources very
    inefficently. As so, it is better to run the 1st/3rd step on a 
    compute node and run the 2nd step on a GPU node.
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Flag file to indicate success of make_examples 
    """
    input: 
        tumor  = join(workpath, "BAM", "{name}.sorted.bam"),
        normal = get_normal_sorted_bam
    output:
        success = join(workpath, "deepsomatic", "mk_examples", "{name}.make_examples.success"),
    params: 
        rname  = "ds_mkexamples",
        genome = config['references']['GENOME'],
        tmpdir = tmpdir,
        nshards = int(allocated("threads", "deepsomatic_make_examples", cluster))-1,
        example = lambda w: join(workpath, "deepsomatic", "mk_examples", "{0}.make_examples.tfrecord@{1}.gz".format(
            w.name,
            int(allocated("threads", "deepsomatic_make_examples", cluster))
        )),
        # TODO: add option --ffpe option to pipeline, that
        # selects either the ffpe_wgs or ffpe_wes checkpoints.
        # Building option for checkpoint file (assumes TN-pairs and
        # non-FFPE samples), where:
        #  @WES = "/opt/models/deepsomatic/wes"
        #  @WGS = "/opt/models/deepsomatic/wgs"
        ckpt = lambda _: "/opt/models/deepsomatic/wes" if run_wes else "/opt/models/deepsomatic/wgs",
        # Call variants within regions BED 
        # file created from WES capture kit
        wes_region_option = lambda _: "--regions {0}".format(
            join(workpath, "references", "wes_regions_50bp_padded.bed"),
        ) if run_wes else '',
        # Get tumor and normal sample names 
        tumor  = '{name}',
        # Building option for the paired normal sorted bam
        normal_bam_option = lambda w: "--reads_normal {0}.sorted.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
        # Building option for the normal sample name
        normal_name_option = lambda w: "--sample_name_normal {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: int(allocated("threads", "deepsomatic_make_examples", cluster))
    container: config['images']['deepsomatic']
    envmodules: config['tools']['deepsomatic']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    echo "Using tmpdir: ${{tmp}}"
    export TMPDIR="${{tmp}}"

    # Run DeepSomatic make_examples and
    # parallelize it using gnu-parallel.
    # Exporting OpenBLAS variable to avoid
    # any issues related to RLIMIT_NPROC,
    # that issue will occur when the number
    # of shards is >= 32. Exporting this
    # OpenBLAS variable to also controls
    # the number of threads in a thread
    # pool. By setting this variable to 1,
    # work is done in the thread that ran
    # the operation, rather than disbatching
    # the work to a thread pool. If this 
    # variable is not provided, it can lead
    # to issues related to nested parallelism.
    # See this Github issue for more info:
    # https://github.com/google/deepsomatic/issues/28 
    export OPENBLAS_NUM_THREADS=1
    time seq 0 {params.nshards} \\
        | parallel \\
            --eta \\
            -q \\
            --halt 2 \\
            --line-buffer \\
            make_examples_somatic \\
                --mode calling \\
                --ref {params.genome} \\
                --reads_tumor {input.tumor}  {params.normal_bam_option} \\
                --sample_name_tumor {params.tumor} {params.normal_name_option} \\
                --examples {params.example} \\
                --checkpoint "{params.ckpt}" {params.wes_region_option} \\
                --vsc_max_fraction_indels_for_non_target_sample "0.5" \\
                --vsc_max_fraction_snps_for_non_target_sample "0.5" \\
                --vsc_min_fraction_indels "0.05" \\
                --vsc_min_fraction_snps "0.029" \\
                --task {{}} \\
    && touch {output.success}
    """

if use_gpus:
    # Use GPU-acceleration to speed up 
    # the second step in deepsomatic
    rule deepsomatic_call_variants_gpu:
        """
        Data processing step to call somatic variants using deep neural 
        network. The make_examples step prepares the input data for the
        deepsomatic's CNN. DeepSomatic is an extension of deep learning-
        based variant caller Deepvariant. It is composed of multiple steps
        that takes aligned reads (in BAM or CRAM format), produces pileup
        image tensors from them, classifies each tensor using a convolutional 
        neural network, and finally reports the results in a standard VCF or
        gVCF file. This rule is the first step in the deepsomatic pipeline:
           1. make_examples        (CPU, parallelizable with gnu-parallel)
         * 2. call_variants        (GPU, use a GPU node)
           3. postprocess_variants (CPU, single-threaded)
        Running deepsomatic in a single step using run_deepsomatic is not 
        optimal for large-scale projects as it will consume resources very
        inefficently. As so, it is better to run the 1st/3rd step on a 
        compute node and run the 2nd step on a GPU node.
        NOTE: When deepsomatic is run on a GPU/TPU, it will scatter the
        writing of the output *.call_variants.tfrecord.gz across a pool
        of processes (by default, --writer_threads 16). This causes causes
        the final output file to be different if you are running DeepSomatic
        on a CPU versus GPU.
        @Input:
            Flag file to indicate success of make_examples (scatter)
        @Output:
            Flag file to indicate success of call_variants, 
            actually produces (given 16 writer threads):
              {name}.call_variants-00000-of-00016.tfrecord.gz, ...
              {name}.call_variants-00015-of-00016.tfrecord.gz
        """
        input: 
            success = join(workpath, "deepsomatic", "mk_examples", "{name}.make_examples.success"),
        output:
            success = join(workpath, "deepsomatic", "call_variants", "{name}.cv.success"),
        params: 
            rname  = "ds_callvars_gpu",
            genome = config['references']['GENOME'],
            # Singularity options
            sif = config['images']['deepsomatic_gpu'],
            bindpaths = ','.join(bindpath),
            tmpdir = tmpdir,
            # NOTE: There BE dragons here!
            # We need allocation info from make_examples rule
            # to determine the number of shards that were
            # used in the make_examples step, this is used
            # to resolve a dependency file of this rule,
            # which is the examples tf record file produced by 
            # make_examples. This file gets passed to the
            # --examples option of call_variants. 
            example = lambda w: join(workpath, "deepsomatic", "mk_examples", "{0}.make_examples.tfrecord@{1}.gz".format(
                w.name,
                int(allocated("threads", "deepsomatic_make_examples", cluster))
            )),
            callvar = join(workpath, "deepsomatic", "call_variants", "{name}.call_variants.tfrecord.gz"),
            # TODO: add option --ffpe option to pipeline, that
            # selects either the ffpe_wgs or ffpe_wes checkpoints.
            # Building option for checkpoint file (assumes TN-pairs and
            # non-FFPE samples), where:
            #  @WES = "/opt/models/deepsomatic/wes"
            #  @WGS = "/opt/models/deepsomatic/wgs"
            ckpt = lambda _: "/opt/models/deepsomatic/wes" if run_wes else "/opt/models/deepsomatic/wgs",
        threads: max(int(allocated("threads", "deepsomatic_call_variants_gpu", cluster)) - 2, 4),
        envmodules: config['tools']['deepsomatic']
        shell: """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        echo "Using tmpdir: ${{tmp}}"
        export TMPDIR="${{tmp}}"

        # Run DeepSomatic call_variants
        # using a GPU acceleration
        singularity exec \\
            -c \\
            --nv  \\
            -B {params.bindpaths},${{tmp}}:/tmp \\
            {params.sif} /bin/bash -c \\
        'time call_variants \\
            --outfile {params.callvar} \\
            --examples {params.example} \\
            --checkpoint {params.ckpt} \\
            --writer_threads {threads}'
        touch "{output.success}"
        """
else:
    # Use CPU accelerated version for the
    # second step in deepsomatic
    rule deepsomatic_call_variants_cpu:
        """
        Data processing step to call somatic variants using deep neural 
        network. The make_examples step prepares the input data for the
        deepsomatic's CNN. DeepSomatic is an extension of deep learning-
        based variant caller Deepvariant. It is composed of multiple steps
        that takes aligned reads (in BAM or CRAM format), produces pileup
        image tensors from them, classifies each tensor using a convolutional 
        neural network, and finally reports the results in a standard VCF or
        gVCF file. This rule is the first step in the deepsomatic pipeline:
           1. make_examples        (CPU, parallelizable with gnu-parallel)
         * 2. call_variants        (CPU, multi-threaded)
           3. postprocess_variants (CPU, single-threaded)
        Running deepsomatic in a single step using run_deepsomatic is not 
        optimal for large-scale projects as it will consume resources very
        inefficently. As so, it is better to run the 1st/3rd step on a 
        compute node and run the 2nd step on a GPU node. 
        NOTE: There be dragens here! When deepsomatic is run on a GPU/TPU,
        it will scatter the writing of the output *.call_variants.tfrecord.gz
        across a pool of processes (by default, --writer_threads 16). This 
        causes causes the final output file to be different if you are 
        running DeepSomatic on a CPU versus GPU.
        @Input:
            Flag file to indicate success of make_examples (scatter)
        @Output:
            Flag file to indicate success of call_variants,
            actually produces:
              {name}.call_variants.tfrecord.gz
        """
        input: 
            success = join(workpath, "deepsomatic", "mk_examples", "{name}.make_examples.success"),
        output:
            success = join(workpath, "deepsomatic", "call_variants", "{name}.cv.success"),
        params: 
            rname  = "ds_callvars_cpu",
            genome = config['references']['GENOME'],
            tmpdir = tmpdir,
            # NOTE: There BE dragons here!
            # We need allocation info from make_examples rule
            # to determine the number of shards that were
            # used in the make_examples step, this is used
            # to resolve a dependency file of this rule,
            # which is the examples tf record file produced by 
            # make_examples. This file gets passed to the
            # --examples option of call_variants. 
            example = lambda w: join(workpath, "deepsomatic", "mk_examples", "{0}.make_examples.tfrecord@{1}.gz".format(
                w.name,
                int(allocated("threads", "deepsomatic_make_examples", cluster))
            )),
            callvar = join(workpath, "deepsomatic", "call_variants", "{name}.call_variants.tfrecord.gz"),
            # TODO: add option --ffpe option to pipeline, that
            # selects either the ffpe_wgs or ffpe_wes checkpoints.
            # Building option for checkpoint file (assumes TN-pairs and
            # non-FFPE samples), where:
            #  @WES = "/opt/models/deepsomatic/wes"
            #  @WGS = "/opt/models/deepsomatic/wgs"
            ckpt = lambda _: "/opt/models/deepsomatic/wes" if run_wes else "/opt/models/deepsomatic/wgs",
        threads: int(allocated("threads", "deepsomatic_call_variants_cpu", cluster))
        container: config['images']['deepsomatic']
        envmodules: config['tools']['deepsomatic']
        shell: """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT
        echo "Using tmpdir: ${{tmp}}"
        export TMPDIR="${{tmp}}"

        # Run CPU DeepSomatic call_variants
        time call_variants \\
            --outfile {params.callvar} \\
            --examples {params.example} \\
            --checkpoint {params.ckpt}
        touch "{output.success}"
        """


rule deepsomatic_postprocess_variants:
    """
    Data processing step to call somatic variants using deep neural 
    network. The make_examples step prepares the input data for the
    deepsomatic's CNN. DeepSomatic is an extension of deep learning-
    based variant caller Deepvariant. It is composed of multiple steps
    that takes aligned reads (in BAM or CRAM format), produces pileup
    image tensors from them, classifies each tensor using a convolutional 
    neural network, and finally reports the results in a standard VCF or
    gVCF file. This rule is the first step in the deepsomatic pipeline:
       1. make_examples        (CPU, parallelizable with gnu-parallel)
       2. call_variants        (GPU, use a GPU node)
     * 3. postprocess_variants (CPU, single-threaded)
    Running deepsomatic in a single step using run_deepsomatic is not 
    optimal for large-scale projects as it will consume resources very
    inefficently. As so, it is better to run the 1st/3rd step on a 
    compute node and run the 2nd step on a GPU node.
    NOTE: There be dragens here! Deepsomatic will generate a different 
    set of output files at the call_variants steps if it is run on a 
    CPU versus a GPU. A flag file is used to indicate this step was
    successful and the actual sharded/non-shared file (which is input
    to this step) is resolved in the params section. Please note this
    file will not actually exist if call_variants was run with a GPU.
    Looking at their source code, it appears deepsomatic has some logic
    to detect if a sharded writer was used in the previous step, and it
    will read in the set of sharded call_variants files without issues. 
    @Input:
        Per-sample call_variants tensorflow records file (scatter)
    @Output:
        Single-sample VCF file with called variants
    """
    input: 
        success = join(workpath, "deepsomatic", "call_variants", "{name}.cv.success"),
    output:
        vcf  = join(workpath, "deepsomatic", "somatic", "{name}.deepsomatic.vcf"),
    params: 
        rname   = "ds_postprovars",
        genome  = config['references']['GENOME'],
        tmpdir  = tmpdir,
        callvar = join(workpath, "deepsomatic", "call_variants", "{name}.call_variants.tfrecord.gz"),
    threads: int(allocated("threads", "deepsomatic_postprocess_variants", cluster))
    container: config['images']['deepsomatic']
    envmodules: config['tools']['deepsomatic']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    echo "Using tmpdir: ${{tmp}}"
    export TMPDIR="${{tmp}}"

    # Run DeepSomatic postprocess_variants
    time postprocess_variants \\
        --ref {params.genome} \\
        --infile {params.callvar} \\
        --outfile {output.vcf} \\
        --process_somatic=true \\
        --cpus={threads}
    """


rule muse:
    """Data-processing step to call somatic mutations with MuSE. This tool is 
    unique in accounting for tumor heterogeneity using a sample-specific error 
    model that improves sensitivity and specificity in mutation calling from 
    sequencing data.
    More information about MuSE can be found here: 
    https://github.com/wwylab/MuSE
    @Input:
        Realigned, recalibrated BAM file for a normal in a TN pair
    @Output:
        Per sample VCF of somatic variants
    """
    input: 
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        txt    = join(workpath, "muse", "somatic", "{name}.MuSE.txt"),
        vcf    = temp(join(workpath, "muse", "somatic", "{name}.muse.tmp.vcf")),
        header = temp(join(workpath, "muse", "somatic", "{name}.samples")),
        final  = join(workpath, "muse", "somatic", "{name}.muse.vcf"),
    params:
        tumor  = join(workpath, "muse", "somatic", "{name}"),
        rename = "{name}",
        rname  = 'muse',
        genome = config['references']['GENOME'],
        dbsnp  = config['references']['DBSNP'],
        # Building optional argument for paired normal
        normal_option = lambda w: "{0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
        # Creating optional reheader for paired normal,
        # resolves to "\nNORMAL\t${normalName}"
        normal_header = lambda w: "\\nNORMAL\\t{0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
        # Building option for WGS/WES,
        # if WES then "muse sump -E"
        # else (WGS) set "muse sump -G"
        muse_sump_option = lambda _: "-E" if run_wes else '-G',
    threads: 
        # MuSE over-alllocates threads,
        # see this issue for more info:
        # https://github.com/wwylab/MuSE/issues/8
        max(int(allocated("threads", "muse", cluster)) - 9, 4),
    container: config['images']['genome-seek_somatic']
    envmodules:
        config['tools']['muse'],
        config['tools']['bcftools']
    shell: """
    # Prefilter and calculate position
    # specific summary statistics 
    MuSE call \\
        -n {threads} \\
        -f {params.genome} \\
        -O {params.tumor} \\
        {input.tumor} {params.normal_option} 
    # Calculate cutoffs from a 
    # sample specific error model
    MuSE sump \\
        -n {threads} \\
        {params.muse_sump_option} \\
        -I {output.txt} \\
        -O {output.vcf} \\
        -D {params.dbsnp}
    # Renaming TUMOR/NORMAL in VCF 
    # with real sample names
    echo -e "TUMOR\\t{params.rename}{params.normal_header}" \\
    > {output.header} 
    bcftools reheader \\
        -o {output.final} \\
        -s {output.header} \\
        {output.vcf}
    """


rule strelka:
    """Data-processing step to call somatic mutations with Strelka. This tool is 
    optimized for rapid clinical analysis of germline variation in small cohorts 
    and somatic variation in tumor/normal sample pairs. More information about 
    strelka can be found here: https://github.com/Illumina/strelka
    @Input:
        Realigned, recalibrated BAM file for a normal in a TN pair
    @Output:
        Per sample VCF of somatic variants
    """
    input: 
        tumor  = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        snps   = join(workpath, "strelka", "{name}", "results", "variants", "somatic.snvs.vcf.gz"),
        indels = join(workpath, "strelka", "{name}", "results", "variants", "somatic.indels.vcf.gz"),
        vcf    = temp(join(workpath, "strelka", "{name}", "{name}.tmp.vcf")),
        header = temp(join(workpath, "strelka", "{name}", "{name}.samples")),
        rehead  = temp(join(workpath, "strelka", "somatic", "{name}.rehead.vcf")),
    params:
        tmpdir = tmpdir,
        tumor  = '{name}',
        rname  = 'strelka',
        purple_jar = config['references']['HMFTOOLS_PURPLE_JAR'],
        outdir = join(workpath, "strelka", "{name}"),
        workflow = join(workpath, "strelka", "{name}", "runWorkflow.py"),
        genome   = config['references']['GENOME'],
        pon      = config['references']['PON'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "strelka", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "strelka", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "strelka", cluster).lower().rstrip('g'),
        # Building optional argument for paired normal
        normal_option = lambda w: "--normalBam {0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
        # Creating optional reheader for paired normal,
        # resolves to "\nNORMAL\t${normalName}"
        normal_header = lambda w: "\\nNORMAL\\t{0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
        # Building option for WGS/WES, if WES use padded 
        # WES BED file, else use manta default
        regions = lambda _: "{0}".format(
            join(workpath, "references", "wes_regions_50bp_padded.bed.gz"),
        ) if run_wes else config['references']['MANTA_CALLREGIONS'],
        # Building option for WES flag
        wes = lambda _: "--exome" if run_wes else "",
    threads: 
        max(int(allocated("threads", "strelka", cluster))-1, 1)
    container: config['images']['genome-seek_somatic']
    envmodules:
        config['tools']['strelka'],
        config['tools']['gatk3'],
        config['tools']['bcftools']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    # Delete previous attempts output
    # directory to ensure hard restart
    if [ -d "{params.outdir}" ]; then
        rm -rf "{params.outdir}"
    fi

    # Configure Strelka somatic workflow
    configureStrelkaSomaticWorkflow.py \\
        --referenceFasta {params.genome} \\
        --tumorBam {input.tumor} {params.normal_option} \\
        --runDir {params.outdir} \\
        --callRegions {params.regions} {params.wes}
    
    # Call somatic variants with Strelka
    echo "Starting Strelka workflow..."
    {params.workflow} \\
        -m local \\
        -j {threads} \\
        -g {params.memory} 
    
    # Combine and filter results
    echo "Running CombineVariants..."
    java -Xmx{params.memory}g -Djava.io.tmpdir=${{tmp}} \
        -XX:ParallelGCThreads={threads} -jar $GATK_JAR -T CombineVariants \\
        --use_jdk_inflater --use_jdk_deflater \\
        -R {params.genome} \\
        --variant {output.snps} \\
        --variant {output.indels} \\
        --assumeIdenticalSamples \\
        --filteredrecordsmergetype KEEP_UNCONDITIONAL \\
        -o {output.vcf}
    
    # Renaming TUMOR/NORMAL in 
    # VCF with real sample names
    echo -e "TUMOR\\t{params.tumor}{params.normal_header}" \\
    > {output.header} 
    bcftools reheader \\
        -o {output.rehead} \\
        -s {output.header} \\
        {output.vcf}
    """  


rule strelka_format:
    """Data-processing step to call somatic mutations with Strelka. This tool is 
    optimized for rapid clinical analysis of germline variation in small cohorts 
    and somatic variation in tumor/normal sample pairs. More information about 
    strelka can be found here: https://github.com/Illumina/strelka
    @Input:
        Realigned, recalibrated BAM file for a normal in a TN pair
    @Output:
        Per sample VCF of somatic variants
    """
    input: 
        rehead  = join(workpath, "strelka", "somatic", "{name}.rehead.vcf"),
    output:
        final  = join(workpath, "strelka", "somatic", "{name}.strelka.vcf"),
    params:
        tmpdir = tmpdir,
        tumor  = '{name}',
        rname  = 'strelka_format',
        purple_jar = config['references']['HMFTOOLS_PURPLE_JAR'],
        outdir = join(workpath, "strelka", "{name}"),
        workflow = join(workpath, "strelka", "{name}", "runWorkflow.py"),
        regions  = config['references']['MANTA_CALLREGIONS'],
        genome   = config['references']['GENOME'],
        pon      = config['references']['PON'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "strelka_format", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "strelka_format", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "strelka_format", cluster).lower().rstrip('g'),
        # Building optional argument for paired normal
        normal_option = lambda w: "--normalBam {0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
        # Creating optional reheader for paired normal,
        # resolves to "\nNORMAL\t${normalName}"
        normal_header = lambda w: "\\nNORMAL\\t{0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: 
        max(int(allocated("threads", "strelka_format", cluster))-1, 1)
    container: config['images']['genome-seek_cnv']
    envmodules:
        config['tools']['rlang'],
    shell: """
    # Adding AD annotation to VCF
    java -Xmx{params.memory}g -cp {params.purple_jar} \\
        com.hartwig.hmftools.purple.tools.AnnotateStrelkaWithAllelicDepth \\
        -in {input.rehead} \\
        -out {output.final}
    """


rule somatic_selectvar:
    """Data-processing step to post-process vcf file generated by all the
    somatic callers. This step takes the somatic calls from all the callers
    (assumes already re-headered if needed, i.e. strelka and muse), and then 
    runs bcftools norm to split multi-allelic sites AND gatk SelectVariants 
    to filter sites. For WES data, this step will also remove any variants
    that are outside the padded regions/capture-kit BED file.
    @Input:
        Per sample, per caller, VCF somatic variants
    @Output:
        Per sample, per caller, VCF filtered and normalized somatic SNVs and indels
    """
    input:
        vcf  = join(workpath, "{caller}", "somatic", "{name}.{caller}.vcf"),
        bed = provided(join(workpath, "references", "wes_regions_50bp_padded.bed.gz"), run_wes),
    output:
        norm = join(workpath, "{caller}", "somatic", "{name}.{caller}.norm.vcf"),
        ngz  = join(workpath, "{caller}", "somatic", "{name}.{caller}.norm.vcf.gz"),
        tbi  = join(workpath, "{caller}", "somatic", "{name}.{caller}.norm.vcf.gz.tbi"),
        filt = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.vcf"),
    params:
        rname  = 'somselect',
        genome = config['references']['GENOME'],
        pon    = config['references']['PON'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "somatic_selectvar", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "somatic_selectvar", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "somatic_selectvar", cluster).lower().rstrip('g'),
        # Building intervals option for WES
        wes_intervals_option = lambda _: "--intervals {0}".format(
            join(workpath, "references", "wes_regions_50bp_padded.bed.gz"),
        ) if run_wes else '',
    threads: 
        int(allocated("threads", "somatic_selectvar", cluster))
    container: config['images']['genome-seek_somatic']
    envmodules:
        config['tools']['gatk4'],
        config['tools']['bcftools'],
    shell: """
    # Normalize VCF prior to SelectVar
    # which needs multi-allelic sites
    # to be split prior to running 
    echo "Running bcftools norm..."
    bcftools norm \\
        -c w \\
        -m - \\
        -Ov \\
        --threads {threads} \\
        -f {params.genome} \\
        -o {output.norm} \\
        {input.vcf}
    # Index normalized VCF for GATK
    # SelectVariants, which requires
    # a tabix index-ed bgzip-ed VCF
    # if the --interval option is used
    bgzip -c \\
        {output.norm} \\
    > {output.ngz}
    tabix -p vcf \\
        {output.ngz}
    # Remove filtered sites and output
    # variants not called in the PON
    echo "Running SelectVariants..."
    gatk --java-options "-Xmx{params.memory}g" SelectVariants \\
        -R {params.genome} \\
        --variant {output.ngz} \\
        --discordance {params.pon} \\
        --exclude-filtered {params.wes_intervals_option} \\
        --output {output.filt}
    # Fix format number metadata, gatk 
    # SelectVariants converts Number
    # metadata incorrectly when it
    # it is set to Number=.
    sed -i 's/Number=R/Number=./g' \\
        {output.filt}
    """


rule somatic_merge_tumor:
    """Data-processing step to post-process vcf file generated by all the
    somatic callers. This step takes filtered tumor sample callsets from 
    each caller and intersects them so only variants found in at least 
    2 callers are retained. bcftools isec needs at least 2 files as input.
    @Input:
        Somatic variants found in the tumor sample (gather-across-somatic-callers)
    @Output:
        Variants found in at least 2 callers
    """
    input:
        tn_callers = get_somatic_tn_callers,   # i.e muse, strelka, deepsomatic
        octopus = join(workpath, "octopus", "somatic", "{name}.octopus.filtered.norm.vcf"),
        mutect2 = join(workpath, "mutect2", "somatic", "{name}.mutect2.filtered.norm.vcf"),
    output:
        merged = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.vcf.gz"),
    params:
        rname  = 'tumormerge',
        genome = config['references']['GENOME'],
        # For UGE/SGE clusters memory is allocated
        # per cpu, so we must calculate total mem
        # as the product of threads and memory
        memory    = lambda _: int(
            int(allocated("mem", "somatic_merge_tumor", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "somatic_merge_tumor", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "somatic_merge_tumor", cluster).lower().rstrip('g'),
        tmpdir = tmpdir,
        # Dynamically update the priority list
        # based on wether a sample is a tumor-only
        # or a tumor-normal
        priority_list = lambda w: "mutect2,octopus,muse,strelka" \
          if tumor2normal[w.name] else "mutect2,octopus",
        strelka_option = lambda w: "--variant:strelka {0}.strelka.filtered.norm.vcf".format(
            join(workpath, "strelka", "somatic", w.name)
        ) if tumor2normal[w.name] else "",
        muse_option = lambda w: "--variant:muse {0}.muse.filtered.norm.vcf".format(
            join(workpath, "muse", "somatic", w.name)
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "somatic_merge_tumor", cluster))
    container: config['images']['genome-seek_somatic']
    envmodules: config['tools']['gatk3']
    shell: """
    # Intersect somatic callset to find
    # variants in at least two callers
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    java -Xmx{params.memory}g -Djava.io.tmpdir=${{tmp}} \\
        -XX:ParallelGCThreads={threads} -jar $GATK_JAR -T CombineVariants \\
            -R {params.genome} \\
            --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \\
            --genotypemergeoption PRIORITIZE \\
            --rod_priority_list {params.priority_list} \\
            -o {output.merged} \\
            --variant:octopus {input.octopus} --variant:mutect2 {input.mutect2} {params.strelka_option} {params.muse_option} 
    """

rule somatic_sample_maf:
    """Data-processing step to convert the merged somatic calls from 
    each caller into a MAF file. This step takes filtered, norm, tumor 
    callset from all callers and annotates the variants with VEP/106
    and converts the resulting VCF file into MAF file format. vcf2maf 
    requires the input vcf file is NOT compressed.
    @Input:
        Somatic variants found in the tumor sample (scatter-per-sample)
    @Output:
        Annotated, merged MAF file contaning somatic callsets
    """
    input:
        vcf = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.vcf.gz"),
    output:
        vcf = temp(join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.temp.vcf")),
        vep = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.temp.vep.vcf"),
        maf = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.maf"),
    params:
        rname  = 'samplemaf',
        tumor  = '{name}',
        memory = allocated("mem", "somatic_sample_maf", cluster).rstrip('G'),
        vep_data    = config['references']['VEP_DATA'],
        vep_build   = config['references']['VEP_BUILD'],
        vep_species = config['references']['VEP_SPECIES'],
        ref_version = config['references']['VEP_REF_VERSION'],
        genome      = config['references']['GENOME'],
        # Building optional argument for paired normal
        normal_option = lambda w: "--normal-id {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "somatic_sample_maf", cluster))
    container: config['images']['vcf2maf']
    shell: """
    # vcf2maf needs an uncompressed VCF file
    zcat {input.vcf} \\
    > {output.vcf}
    # Run VEP and convert VCF into MAF file
    vcf2maf.pl \\
        --input-vcf {output.vcf} \\
        --output-maf {output.maf} \\
        --vep-path ${{VEP_HOME}} \\
        --vep-data {params.vep_data} \\
        --cache-version {params.ref_version} \\
        --ref-fasta {params.genome} \\
        --vep-forks {threads} \\
        --tumor-id {params.tumor} {params.normal_option} \\
        --ncbi-build {params.vep_build} \\
        --species {params.vep_species} \\
        --retain-info set
    """


rule somatic_cohort_maf:
    """Data-processing step to merge the per-sample MAF files into a 
    single MAF file for all samples, a cohort MAF file. 
    @Input:
        Somatic tumor MAF files (gather-per-sample)
    @Output:
        Merged somatic tumor MAF, cohort-level, with all call sets
    """
    input: 
        mafs = expand(join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.maf"), name=tumors),
    output: 
        maf  = join(workpath, "merged", "somatic", "cohort_somatic_variants.maf"),
    params:
        rname = 'cohortmaf',
        memory = allocated("mem", "somatic_cohort_maf", cluster).rstrip('G'),
    threads: 
        int(allocated("threads", "somatic_cohort_maf", cluster))
    container: config['images']['vcf2maf']
    shell: """
    echo "Combining MAFs..."
    head -2 {input.mafs[0]} > {output.maf}
    awk 'FNR>2 {{print}}' {input.mafs} >> {output.maf}
    """


rule somatic_cohort_maftools:
    """Data-processing step to run maftools on merged cohort-level 
    MAF file, produces summarized plots like an oncoplot. 
    @Input:
        Cohort-level somatic MAF file (indirect-gather-due-to-aggregation) 
    @Output:
        TCGA comparsion plot
        Top 20 genes by Vaf plot
        MAF Summary plot
        Oncoplot
    """
    input: 
        maf  = join(workpath, "merged", "somatic", "cohort_somatic_variants.maf"),
    output:
        tcga  = join(workpath, "merged", "somatic", "cohort_tcga_comparison.pdf"),
        gvaf  = join(workpath, "merged", "somatic", "cohort_genes_by_VAF.pdf"),
        summary  = join(workpath, "merged", "somatic", "cohort_maf_summary.pdf"),
        oncoplot = join(workpath, "merged", "somatic", "cohort_oncoplot.pdf"),
    params:
        rname = 'maftools',
        wdir  = join(workpath, "merged", "somatic"),
        memory = allocated("mem", "somatic_cohort_maftools", cluster).rstrip('G'),
        script = join("workflow", "scripts", "maftools.R"),
    threads: 
        int(allocated("threads", "somatic_cohort_maftools", cluster)),
    container: config['images']['genome-seek_somatic']
    envmodules: config['tools']['rlang']
    shell: """
    Rscript {params.script} \\
        {params.wdir} \\
        {input.maf} \\
        {output.summary} \\
        {output.oncoplot}
    """


rule somatic_sample_sigprofiler:
    """Data-processing step to run sigprofiler on each sample's MAF file, 
    produces summarized sample protrait plot. SigProfiler is scattered per
    sample to avoid issues where a failure with one specific sample in a
    multi-sample MAF file prevents the final plot from being generated.
    @IMPORTANT NOTICE
        At the current moment, the docker image for this tool only contains 
    the reference files for GRCh38/hg38. The SigProfilerMatrixGenerator step 
    of SigProfiler requires reference files; however, the tool installs these
    reference filed in the python package site installation location. This is 
    not ideal for several reasons. I have reached out to the author of the tool 
    about decoupling the reference file installation from the python-site 
    installation location. This would allow also to include the reference 
    files in the resource bundle, AND it would avoid the problem of rebuilding 
    the docker image to support new reference genomes, AND issue related to the 
    docker image's size. Here is a link to the Github issue reference above:
    https://github.com/AlexandrovLab/SigProfilerMatrixGenerator/issues/101
    For more information about SigProfiler, please visit Github: 
    https://github.com/AlexandrovLab/SigProfilerPlotting
    @Input:
        Single sample somatic MAF file (scatter-per-sample) 
    @Output:
        SigProfiler Sample Portrait Plot (pdf)
    """
    input: 
        maf = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.maf"),
    output:
        pdf = join(workpath, "sigprofiler", "sample_portrait_{name}.pdf"),
    params:
        rname = 'samplesigpro',
        sample = '{name}',
        wdir   = join(workpath, "sigprofiler", "{name}"),
        odir   = join(workpath, "sigprofiler"),
        memory = allocated("mem", "somatic_sample_sigprofiler", cluster).rstrip('G'),
        script = join("workflow", "scripts", "sigprofiler.py"),
        genome = config['references']['SIGPROFILER_GENOME']
    threads: 
        int(allocated("threads", "somatic_sample_sigprofiler", cluster)),
    container: config['images']['sigprofiler']
    shell: """
    # SigProfiler input directory must
    # only contain input MAF
    mkdir -p "{params.wdir}"
    ln -sf {input.maf} {params.wdir}
    python3 {params.script} \\
        -i {params.wdir}/ \\
        -o {params.odir}/ \\
        -p {params.sample} \\
        -r {params.genome}
    """


rule somatic_cohort_sigprofiler:
    """Data-processing step to merge each samples SigProfiler Portrait
    plots into one PDF file. SigProfiler is scattered per sample to avoid 
    issues where a failure with one specific sample in a multi-sample MAF 
    file prevents the final plot from being generated. This way of running
    the tool is more fault-tolerant. For more information about SigProfiler, 
    please visit their repo: https://github.com/AlexandrovLab/SigProfilerPlotting
    @Input:
        SigProfiler Sample Portrait Plots (gather-per-sample) 
    @Output:
        Merged SigProfiler Portrait Plot
    """
    input: 
        pdfs = expand(join(workpath, "sigprofiler", "sample_portrait_{name}.pdf"), name=tumors),
    output:
        pdf = join(workpath, "sigprofiler", "merged_sigprofiler.pdf"),
    params:
        rname = 'mergesigpro',
        memory = allocated("mem", "somatic_cohort_sigprofiler", cluster).rstrip('G'),
    threads: 
        int(allocated("threads", "somatic_cohort_sigprofiler", cluster)),
    container: config['images']['sigprofiler']
    shell: """
    # Merge SigProfiler PDFs
    pdfunite {input.pdfs} \\
        {output.pdf}
    """
