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


rule gatk_haplotypecaller:
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
        Single-sample GVCF file with raw, unfiltered SNP and indel calls. 
    """
    input: 
        bam = join(workpath, "BAM", "{name}.recal.bam"),
    output: 
        gvcf = temp(join(workpath, "haplotypecaller", "gVCFs", "{name}.{chrom}.g.vcf.gz")),
        idx = temp(join(workpath, "haplotypecaller", "gVCFs", "{name}.{chrom}.g.vcf.gz.tbi")),
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
            int(allocated("mem", "gatk_haplotypecaller", cluster).lower().rstrip('g')) * \
            int(allocated("threads", "gatk_haplotypecaller", cluster)) 
        )-1 if run_mode == "uge" \
        else allocated("mem", "gatk_haplotypecaller", cluster).lower().rstrip('g'),
    message: "Running GATK4 HaplotypeCaller on '{input.bam}' input file"
    threads: int(allocated("threads", "gatk_haplotypecaller", cluster))
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
