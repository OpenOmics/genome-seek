# Functions and rules for processing data
from scripts.common import abstract_location, allocated, joint_option


# Data processing rules following GATK best practices for
# base quality recalibration (BQSR) and indel realignment
rule gatk_realign:
    """
    Data-processing step to realign indels. This step has been removed 
    in GATK/4.X as it has been integrated into the GATK variant callers; 
    however, regardless of how you perform realignment, it is still 
    preferred. More discussion on this topic can be found here:
    https://github.com/broadinstitute/gatk/issues/3104
    """
    input:
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai")
    output:
        bam = temp(join(workpath, "BAM", "{name}.realign.bam")),
        intervals = join(workpath, "BAM", "{name}.intervals"),
    params: 
        genome = config['references']['GENOME'],
        knowns = joint_option('-known', config['references']['KNOWNINDELS']),
        memory = allocated("mem", "gatk_realign", cluster).lower().rstrip('g'),
        rname = "realign"
    threads: 
        int(allocated("threads", "gatk_realign", cluster))
    envmodules: 
        config['tools']['gatk3'],
    shell: """
    java -Xmx{params.memory}g -jar ${{GATK_JAR}} -T RealignerTargetCreator \\
        -I {input.bam} \\
        -R {params.genome} \\
        -o {output.intervals} \\
        {params.knowns}
    
    java -Xmx{params.memory}g -jar ${{GATK_JAR}} -T IndelRealigner \\
        -R {params.genome} \\
        -I {input.bam} \\
        {params.knowns} \\
        --use_jdk_inflater \\
        --use_jdk_deflater \\
        -targetIntervals {output.intervals} \\
        -o {output.bam}
    """


rule gatk_scatter_recal:
    """
    Scatters base quality recalibration (BQSR), part of the GATK Best Practices.
    The idea is that each sequencer/run will have systematic biases.
    BQSR learns about these biases using known sites of variation (common SNPs),
    and uses it to adjust base quality on all sites, including novel sites of
    variation.  Since base quality is taken into account during variant calling,
    this will help pick up real variants in low depth or otherwise noisy loci.
    @Input:
        Aligned reads in BAM format (scatter-per-sample-per-chrset)
    @Output:
        Recalibration table for BQSR
    """
    input:
        bam = join(workpath, "BAM", "{name}.realign.bam"),
    output:
        recal =  join(workpath, "BAM", "{name}_{recal}_data.grp")
    params: 
        genome = config['references']['GENOME'],
        knowns = joint_option('--known-sites', config['references']['KNOWNRECAL']),
        intervals = lambda w: joint_option('-L', config['references']['NRECALS'][w.recal]),
        memory = allocated("mem", "gatk_scatter_recal", cluster).lower().rstrip('g'),
        rname = '{recal}'
    envmodules:
        config['tools']['gatk4'],
    threads: int(allocated("threads", "gatk_scatter_recal", cluster))
    shell: """
    gatk --java-options '-Xmx{params.memory}g' BaseRecalibrator \\
        --input {input.bam} \\
        --reference {params.genome} \\
        {params.knowns} \\
        --output {output.recal} \\
        {params.intervals}
    """


rule gatk_gather_recal:
    """
    Gathers base quality recalibration (BQSR), part of the GATK Best Practices.
    The BaseRecalibrator is scatter on sets of chromosomes for a given sample 
    to decrease the overall runtime of this step, since it is slow. This rule 
    gathers the per-sample-per-chrset BaseRecalibrator results prior to applying 
    the base recalibration with ApplyBQSR.  
    @Input:
        Recalibration tables for BQSR (gather-per-sample-per-chrset)
    @Output:
        Merged recalibration table for ApplyBQSR 
    """
    input:
        grps  =  expand(join(workpath, "BAM", "{{name}}_{recal}_data.grp"), recal=recals),
    output:
        lsl   = join(workpath, "BAM", "{name}.recals.list"),
        recal = join(workpath, "BAM", "{name}_gathered_recal_data.grp")
    params: 
        sample = "{name}",
        bams   = join(workpath, "BAM"),
        memory = allocated("mem", "gatk_gather_recal", cluster).lower().rstrip('g'),
        rname = 'gather_recal'
    envmodules:
        config['tools']['gatk4'],
    threads: int(allocated("threads", "gatk_gather_recal", cluster))
    shell: """
    # Create GatherBQSR list
    find {params.bams} -iname '{params.sample}_recal*_data.grp' \\
        > {output.lsl}
    # Gather per sample BQSR results
    gatk --java-options '-Xmx{params.memory}g' GatherBQSRReports \\
        --use-jdk-inflater --use-jdk-deflater \\
        -I {output.lsl} \\
        -O {output.recal}
    """
