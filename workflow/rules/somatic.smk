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
        return join(workpath, "mutect2", "{0}.normalPileup.table".format(tumor))
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
        model = config['references']['OCTOPUS_FOREST_MODEL'],
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
    octopus --threads {threads} \\
        -C cancer \\
        --working-directory {params.wd} \\
        --temp-directory-prefix {params.tmpdir} \\
        -R {params.genome} \\
        -I {input.normal} {input.tumor} {params.normal_option} \\
        -o {output.vcf} \\
        --somatic-forest-model {params.model} \\
        --sequence-error-model {params.error} \\
        --annotations AC AD AF \\
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
        vcf  = join(workpath, "octopus", "somatic", "{name}.vcf"),
        lsl  = join(workpath, "octopus", "somatic", "{name}.list"),
        norm = join(workpath, "octopus", "somatic", "{name}.norm.vcf.gz"),
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
    # Normalize Octopus VCF
    bcftools norm \\
        -c w \\
        -m - \\
        -Oz \\
        --threads {threads} \\
        -f {params.genome} \\
        -o {output.norm} \\
        {output.vcf}
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
        vcf = join(workpath, "mutect2", "{name}_mutect2.vcf"),
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
        orien = join(workpath, "mutect2", "{name}.read-orientation-model.tar.gz"),
        stats = join(workpath, "mutect2", "{name}.vcf.stats"),
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
        summary = join(workpath, "mutect2", "{name}.tumorPileup.table"),
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
        summary = join(workpath, "mutect2", "{name}.normalPileup.table"),
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
        tumor  = join(workpath, "mutect2", "{name}.tumorPileup.table"),
        normal = get_normal_pileup_table,
    output:
        summary = join(workpath, "mutect2", "{name}.contamination.table"),
    params:
        tumor  = '{name}',
        rname  = 'calcontaim',
        genome = config['references']['GENOME'],
        # Building optional argument for paired normal
        normal_option = lambda w: "--matched-normal {0}".format(
            join(workpath, "mutect2", "{0}.normalPileup.table".format(w.name))
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
    using the learned read orientation model. VCFtools is used to further filter 
    the VCF file remove any sites that have a filter tag. More information about 
    Mutect2 can be found on the Broad's website: https://gatk.broadinstitute.org/
    @Input:
        Per sample somatic variants from CombineVariants
        Per sample stats file from MergeMutectStats
        Per sample read orientation from LearnReadOrientationModel
        Estimated contamination table from CalculateContamination 
    @Output:
        Per sample VCF of filtered somatic SNVs and indels
        Per sample VCF of final, filtered somatic SNVs and indels
    """
    input:
        vcf   = join(workpath, "mutect2", "{name}_mutect2.vcf"),
        orien = join(workpath, "mutect2", "{name}.read-orientation-model.tar.gz"),
        stats = join(workpath, "mutect2", "{name}.vcf.stats"),
        summary = join(workpath, "mutect2", "{name}.contamination.table"),
    output:
        vcf   = join(workpath, "mutect2", "{name}.filtered.vcf"),
        final = join(workpath, "mutect2", "{name}.filtered.final.vcf"),
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
    # Remove remaining sites 
    # with a filter tag 
    vcftools \\
        --vcf {output.vcf} \\
        --recode \\
        --remove-filtered-all \\
        --stdout \\
    > {output.final}
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
        txt = join(workpath, "muse", "{name}.MuSE.txt"),
        vcf = join(workpath, "muse", "{name}.muse.vcf"),
    params:
        tumor  = join(workpath, "muse", "{name}"),
        rname  = 'muse',
        genome = config['references']['GENOME'],
        dbsnp  = config['references']['DBSNP'],
        # Building optional argument for paired normal
        normal_option = lambda w: "{0}.recal.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "muse", cluster))
    envmodules:
        config['tools']['muse']
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
        vcf    = join(workpath, "strelka", "{name}", "{name}.vcf"),
        filt   = join(workpath, "strelka", "{name}", "{name}.filtered.vcf"),
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
    threads: 
        int(allocated("threads", "strelka", cluster))
    envmodules:
        config['tools']['strelka'],
        config['tools']['gatk3']
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
    
    echo "Running SelectVariants..."
    GATK -m {params.memory}G SelectVariants \\
        -R {params.genome} \\
        --variant {output.vcf} \\
        --discordance {params.pon} \\
        --excludeFiltered \\
        -o {output.filt}
    """
