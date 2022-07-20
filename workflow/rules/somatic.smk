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

# Data processing rules for calling somatic variants
rule octopus:
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
        tumor = join(workpath, "BAM", "{name}.recal.bam"),
        normal = get_normal_recal_bam,
    output:
        vcf = join(workpath, "octopus", "chunks", "{region}", "{name}.vcf.gz"),
    params: 
        genome = config['references']['GENOME'],
        rname  = "octopus",
        chunk = "{region}",
        tumor = "{name}",
        wd = workpath,
        tmpdir = join(workpath, "octopus", "chunks", "{region}", "{name}_tmp"),
        model = config['references']['OCTOPUS_FOREST_MODEL'],
        error = config['references']['OCTOPUS_ERROR_MODEL'],
        # Building optional argument for paired normal 
        normal_option = lambda w: "--normal-sample {0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: 
        int(allocated("threads", "octopus", cluster))
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
        vcfs = expand(join(workpath, "octopus", "chunks", "{region}", "{{name}}.vcf.gz"), region=regions),
    output:
        vcf  = join(workpath, "octopus", "{name}.vcf"),
        lsl  = join(workpath, "octopus", "{name}.list"),
        norm = join(workpath, "octopus", "{name}.norm.vcf.gz"),
    params: 
        genome = config['references']['GENOME'],
        rname  = "octomerge",
        tumor  = "{name}",
        octopath = join(workpath, "octopus", "chunks")
    threads: 
        int(allocated("threads", "octopus_merge", cluster))
    container: 
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
        -m - \\
        -Oz \\
        --threads {threads} \\
        -f {params.genome} \\
        -o {output.norm} \\
        {output.vcf}
    """


rule gatk_mutect2:
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
        tumor = join(workpath, "BAM", "{name}.recal.bam"),
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
        int(allocated("threads", "gatk_mutect2", cluster))
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
