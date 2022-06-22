# Functions and rules for processing data
from scripts.common import abstract_location, allocated, joint_option


# Data processing rules following GATK best practices for
# base quality recalibration (BQSR) and indel realignment
rule gatk_realign:
    """
    Data-processing step to realign indels. This step has been removed 
    in GATK/4.X as it has been integrated into the GATK variant callers; 
    however, However, regardless of how you perform realignment, it is 
    still preferred. More discussion on this topic can be found here:
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


# rule gatk_recal:
#     """
#     Base quality recalibration (BQSR), part of the GATK Best Practices.
#     The idea is that each sequencer/run will have systematic biases.
#     BQSR learns about these biases using known sites of variation (common SNPs),
#     and uses it to adjust base quality on all sites, including novel sites of
#     variation.  Since base quality is taken into account during variant calling,
#     this will help pick up real variants in low depth or otherwise noisy loci.
#     @Input:
#         Aligned reads in BAM format (scatter)
#     @Output:
#         Aligned reads in BAM format, with altered quality scores
#     """
#     input:
#         bam = join(output_bamdir, "preprocessing", "{samples}.raw_map.bam"),
#         bai = join(output_bamdir, "preprocessing", "{samples}.raw_map.bai"),
#     output:
#         bam = join(input_bamdir, "{samples}.input.bam"),
#         re =  temp(join(output_bamdir, "preprocessing", "{samples}_recal_data.grp"))
#     params: 
#         genome = config['references']['GENOME'],
#         knowns = config['references']['KNOWNRECAL'],
#         ver_gatk = config['tools']['gatk4']['version'],
#         chrom = chroms,
#         intervals = intervals_file,
#         rname = 'recal'
#     envmodules:
#         'GATK/4.2.0.0'
#     container:
#         config['images']['wes_base']
#     threads: 24
#     shell: """
#     gatk --java-options '-Xmx48g' BaseRecalibrator \\
#         --input {input.bam} \\
#         --reference {params.genome} \\
#         {params.knowns} \\
#         --output {output.re} \\
#         --intervals {params.intervals}
#     
#     gatk --java-options '-Xmx48g' ApplyBQSR \\
#         --reference {params.genome} \\
#         --input {input.bam} \\
#         --bqsr-recal-file {output.re} \\
#         --output {output.bam} \\
#         --use-jdk-inflater \\
#         --use-jdk-deflater
#     """
