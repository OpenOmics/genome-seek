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
        genome = config['references']['GENOME'],
        rname  = "octosomatic",
        chunk = "{region}",
        tumor = "{name}",
        wd = workpath,
        tmpdir = join("octopus", "somatic", "chunks", "{region}", "{name}_tmp"),
        tmppath = join(workpath, "octopus", "somatic", "chunks", "{region}", "{name}_tmp"),
        model = config['references']['OCTOPUS_SOMATIC_FOREST_MODEL'],
        error = config['references']['OCTOPUS_ERROR_MODEL'],
        # Building optional argument for paired normal 
        normal_option = lambda w: "--normal-sample {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "octopus_somatic", cluster))
    container: 
        config['images']['octopus']
    shell: """
    mkdir -p '{params.tmppath}'
    octopus --threads {threads} \\
        -C cancer \\
        --working-directory {params.wd} \\
        --temp-directory-prefix {params.tmpdir} \\
        -R {params.genome} \\
        -I {input.normal} {input.tumor} {params.normal_option} \\
        -o {output.vcf} \\
        --somatic-forest-model {params.model} \\
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
    output:
        lsl  = join(workpath, "octopus", "somatic", "{name}.list"),
        vcf  = join(workpath, "octopus", "somatic", "{name}.octopus.vcf"),
    params: 
        genome = config['references']['GENOME'],
        rname  = "octomerge",
        tumor  = "{name}",
        octopath = join(workpath, "octopus", "somatic", "chunks")
    threads: 
        int(allocated("threads", "octopus_merge", cluster))
    envmodules:
        config['tools']['bcftools']
    shell: """
    # Create list of chunks to merge
    find {params.octopath} -iname '{params.tumor}.vcf.gz' \\
        > {output.lsl}
    # Merge octopus chunk calls 
    bcftools concat \\
        --threads {threads} \\
        -d exact \\
        -a \\
        -f {output.lsl} \\
        -o {output.vcf} \\
        -O v
    """


rule octopus_germline:
    """
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
        int(allocated("threads", "octopus_normal", cluster))
    container: 
        config['images']['octopus']
    shell: """
    octopus --threads {threads} \\
        --working-directory {params.wd} \\
        --temp-directory-prefix {params.tmpdir} \\
        -R {params.genome} \\
        -I {input.normal} \\
        -o {output.vcf} \\
        --forest-model {params.model} \\
        --annotations AC AD AF \\
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
    envmodules:
        config['tools']['gatk4']
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
    output:
        vcf = temp(join(workpath, "mutect2", "somatic", "{name}.mutect2.tmp.vcf")),
    params:
        rname  = 'merge_mutect2',
        genome = config['references']['GENOME'],
        memory = allocated("mem", "gatk_gather_mutect2", cluster).lower().rstrip('g'),
        # Building optional argument for paired normal
        multi_variant_option = joint_option(
            '--variant',
            expand(join(workpath, "mutect2", "chrom_split", "{{name}}.{chrom}.vcf"), chrom=chunks),
        ),
    threads: 
        int(allocated("threads", "gatk_gather_mutect2", cluster))
    envmodules:
        config['tools']['gatk3']
    shell: """
    GATK -m {params.memory}G CombineVariants \\
        -R {params.genome} \\
        --filteredrecordsmergetype KEEP_UNCONDITIONAL \\
        --assumeIdenticalSamples \\
        -o {output.vcf} \\
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
    envmodules:
        config['tools']['gatk4']
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
        memory = allocated("mem", "gatk_tumorPileup", cluster).lower().rstrip('g'),
    threads: 
        int(allocated("threads", "gatk_tumorPileup", cluster))
    envmodules:
        config['tools']['gatk4']
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
        memory = allocated("mem", "gatk_tumorPileup", cluster).lower().rstrip('g'),
    threads: 
        int(allocated("threads", "gatk_normalPileup", cluster))
    envmodules:
        config['tools']['gatk4']
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
    envmodules:
        config['tools']['gatk4']
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
    envmodules:
        config['tools']['gatk4'],
        config['tools']['vcftools'],
    shell: """
    # Mutect2 orien bias filter
    gatk FilterMutectCalls \\
        -R {params.genome} \\
        -V {input.vcf} \\
        --ob-priors {input.orien} \\
        --contamination-table {input.summary} \\
        -O {output.vcf} \\
        --stats {input.stats} 
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
    threads: 
        int(allocated("threads", "muse", cluster))
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
        -G \\
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
        final  = join(workpath, "strelka", "somatic", "{name}.strelka.vcf"),
    params:
        tumor  = '{name}',
        rname  = 'strelka',
        outdir = join(workpath, "strelka", "{name}"),
        workflow = join(workpath, "strelka", "{name}", "runWorkflow.py"),
        regions  = config['references']['MANTA_CALLREGIONS'],
        genome   = config['references']['GENOME'],
        pon      = config['references']['PON'],
        memory   = allocated("mem", "strelka", cluster).rstrip('G'),
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
        int(allocated("threads", "strelka", cluster))
    envmodules:
        config['tools']['strelka'],
        config['tools']['gatk3'],
        config['tools']['bcftools']
    shell: """
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
        --callRegions {params.regions}
    # Call somatic variants with Strelka
    echo "Starting Strelka workflow..."
    {params.workflow} \\
        -m local \\
        -j {threads} \\
        -g {params.memory} 
    # Combine and filter results
    echo "Running CombineVariants..."
    GATK -m {params.memory}G CombineVariants \\
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
        -o {output.final} \\
        -s {output.header} \\
        {output.vcf}
    """


rule somatic_selectvar:
    """Data-processing step to post-process vcf file generated by all the
    somatic callers. This step takes the somatic calls from all the callers
    (assumes already re-headered if needed, i.e. strelka and muse), and then 
    runs gatk SelectVariants to filter sites AND bcftools norm to split multi-
    allelic sites. 
    @Input:
        Per sample, per caller, VCF somatic variants
    @Output:
        Per sample, per caller, VCF filtered and normalized somatic SNVs and indels
    """
    input:
        vcf  = join(workpath, "{caller}", "somatic", "{name}.{caller}.vcf"),
    output:
        filt = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.vcf"),
        norm = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.vcf"),
    params:
        rname  = 'somselect',
        genome = config['references']['GENOME'],
        pon    = config['references']['PON'],
        memory = allocated("mem", "somatic_selectvar", cluster).rstrip('G'),
    threads: 
        int(allocated("threads", "somatic_selectvar", cluster))
    envmodules:
        config['tools']['gatk4'],
        config['tools']['bcftools'],
    shell: """
    # Remove filtered sites and output
    # variants not called in the PON
    echo "Running SelectVariants..."
    gatk --java-options "-Xmx{params.memory}g" SelectVariants \\
        -R {params.genome} \\
        --variant {input.vcf} \\
        --discordance {params.pon} \\
        --exclude-filtered \\
        --output {output.filt}
    # Normalize VCF 
    echo "Running bcftools norm..."
    bcftools norm \\
        -c w \\
        -m - \\
        -Ov \\
        --threads {threads} \\
        -f {params.genome} \\
        -o {output.norm} \\
        {output.filt}
    """


rule somatic_split_tumor:
    """Data-processing step to post-process vcf file generated by all the
    somatic callers. This step takes a re-header, filtered, and normalized
    somatic vcf file and filters the callset to only contain sites in the 
    tumor sample (needed due to mixed sites in the tumor-normal callset).
    @Input:
        Per sample, per caller somatic variants
    @Output:
        Per sample, per caller somatic variants only found in the tumor sample
    """
    input:
        vcf = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.vcf"),
    output:
        tmp    = temp(join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.tumor.temp.vcf.gz")),
        header = temp(join(workpath, "{caller}", "somatic", "{name}.{caller}.tumor.samples")),
        tumor  = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.tumor.vcf.gz"),
        tbi    = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.tumor.vcf.gz.tbi"),
    params:
        rname  = 'tumorsplit',
        sample = '{name}',
        rename = '{caller}:{name}',
        memory = allocated("mem", "somatic_split_tumor", cluster).rstrip('G'),
    threads: 
        int(allocated("threads", "somatic_split_tumor", cluster))
    envmodules:
        config['tools']['bcftools'],
    shell: """
    # Filter call set to tumor sites
    bcftools view \\
        -c1 \\
        -Oz \\
        -s '{params.sample}' \\
        -o {output.tmp} \\
        {input.vcf}
    # Renaming sample name in VCF 
    # to contain caller name
    echo -e "{params.sample}\\t{params.rename}" \\
    > {output.header} 
    bcftools reheader \\
        -o {output.tumor} \\
        -s {output.header} \\
        {output.tmp}
    # Create an VCF index for intersect
    bcftools index \\
        -f \\
        --tbi \\
        {output.tumor} 
    """


rule somatic_split_normal:
    """Data-processing step to post-process vcf file generated by all the
    somatic callers. This step takes a re-header, filtered, and normalized
    somatic vcf file and filters the callset to only contain sites in the 
    normal sample (needed due to mixed sites in the tumor-normal callset).
    @Input:
        Per sample, per caller somatic variants
    @Output:
        Per sample, per caller somatic variants only found in the normal sample
    """
    input:
        vcf = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.vcf"),
    output:
        tmp    = temp(join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.normal.temp.vcf.gz")),
        header = temp(join(workpath, "{caller}", "somatic", "{name}.{caller}.normal.samples")),
        normal = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.normal.vcf.gz"),
        tbi    = join(workpath, "{caller}", "somatic", "{name}.{caller}.filtered.norm.normal.vcf.gz.tbi"),
    params:
        rname  = 'normalsplit',
        sample = lambda w: "{0}".format(tumor2normal[w.name]),
        rename = lambda w: "{0}:{1}".format(w.caller, tumor2normal[w.name]),
        memory = allocated("mem", "somatic_split_normal", cluster).rstrip('G'),
    threads: 
        int(allocated("threads", "somatic_split_normal", cluster))
    envmodules:
        config['tools']['bcftools'],
    shell: """
    # Filter call set to normal sites
    bcftools view \\
        --force-samples \\
        -c1 \\
        -Oz \\
        -s '{params.sample}' \\
        -o {output.tmp} \\
        {input.vcf}
    # Renaming sample name in VCF 
    # to contain caller name
    echo -e "{params.sample}\\t{params.rename}" \\
    > {output.header} 
    bcftools reheader \\
        -o {output.normal} \\
        -s {output.header} \\
        {output.tmp}
    # Create an VCF index for intersect
    bcftools index \\
        -f \\
        --tbi \\
        {output.normal}
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
        tumors = expand(join(workpath, "{caller}", "somatic", "{{name}}.{caller}.filtered.norm.tumor.vcf.gz"), caller=somatic_callers),
    output:
        lsl    = join(workpath, "merged", "somatic", "{name}.intersect.lsl"),
        merged = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.tumor.vcf.gz"),
        tbi = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.tumor.vcf.gz.tbi"),
    params:
        rname  = 'tumormerge',
        memory = allocated("mem", "somatic_merge_tumor", cluster).rstrip('G'),
        isec_dir = join(workpath, "merged", "somatic", "intersect"),
    threads: 
        int(allocated("threads", "somatic_merge_tumor", cluster))
    envmodules:
        config['tools']['bcftools'],
    shell: """
    # Delete previous attempts output
    # directory to ensure hard restart
    if [ -d "{params.isec_dir}" ]; then
        rm -rf "{params.isec_dir}"
    fi
    # Intersect somatic callset to find
    # variants in at least two callers
    bcftools isec \\
        -Oz \\
        -n=2 \\
        -c none \\
        -p {params.isec_dir} \\
        {input.tumors} 
    # Create list of files to merge 
    find {params.isec_dir} \\
        -name '*.vcf.gz' \\
        | sort \\
    > {output.lsl} 
    # Merge variants found in at 
    # least two callers 
    bcftools merge \\
        -Oz \\
        -o {output.merged} \\
        -l {output.lsl}
    # Create an VCF index for merge
    bcftools index \\
        -f \\
        --tbi \\
        {output.merged}
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
        vcf = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.tumor.vcf.gz"),
    output:
        vcf = temp(join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.tumor.temp.vcf")),
        vep = temp(join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.tumor.temp.vep.vcf")),
        maf = join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.tumor.maf"),
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
    container: 
        config['images']['vcf2maf']
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
        --species {params.vep_species}
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
        mafs = expand(join(workpath, "merged", "somatic", "{name}.merged.filtered.norm.tumor.maf"), name=tumors),
    output: 
        maf  = join(workpath, "merged", "somatic", "cohort_somatic_variants.maf"),
    params:
        rname = 'cohortmaf',
        memory = allocated("mem", "somatic_cohort_maf", cluster).rstrip('G'),
    threads: 
        int(allocated("threads", "somatic_cohort_maf", cluster))
    container:
        config['images']['vcf2maf']
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
    envmodules:
        config['tools']['rlang']
    shell: """
    Rscript {params.script} \\
        {params.wdir} \\
        {input.maf} \\
        {output.summary} \\
        {output.oncoplot}
    """
