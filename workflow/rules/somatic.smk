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
        vcf = join(workpath, "octopus", "chunks", "{chunk}", "{name}.vcf.gz"),
    params: 
        genome = config['references']['GENOME'],
        rname  = "octopus",
        chunk = "{chunk}",
        tumor = "{name}",
        wd = workpath,
        tmpdir = join(workpath, "octopus", "chunks", "{chunk}", "{name}_tmp"),
        model = config['references']['OCTOPUS_FOREST_MODEL'],
        error = config['references']['OCTOPUS_ERROR_MODEL'],
        normal_option = lambda w: "--normal-sample {0}".format(tumor2normal[w.name]) if tumor2normal[w.name] else "",
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
    Data-processing step to merge scatter variant calls from Octopus. Octopus 
    is scattered across genomic intervals or chunks to reduce its overall runtime. 
    @Input:
        Somatic variants in VCF format (gather-per-sample-per-chunks)
    @Output:
        Per sample somatic variants in VCF format  
    """
    input:
        vcfs = expand(join(workpath, "octopus", "chunks", "{chunk}", "{{name}}.vcf.gz"), chunk=chunks),
    output:
        vcf = join(workpath, "octopus", "{name}.vcf"),
        lsl = join(workpath, "octopus", "{name}.list"),
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
        -d exact \\
        -a \\
        -f {output.lsl} \\
        -o {output.vcf} \\
        -O v
    """
