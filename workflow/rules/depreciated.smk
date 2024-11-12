# Depreciated rules that may still be useful for some projects
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


# Depreciated germline variant calling rule(s)
rule deepvariant:
    """
    Data processing step to call germline variants using deep neural 
    network. DeepVariant is a deep learning-based variant caller that 
    takes aligned reads (in BAM or CRAM format), produces pileup image 
    tensors from them, classifies each tensor using a convolutional 
    neural network, and finally reports the results in a standard VCF 
    or gVCF file. This rule runs all three steps in the deepvariant
    pipeline as a single step, i.e.: make_examples, call_variants, and
    postprocess_variants. This is not optimal for large-scale projects
    as it will consume a lot of resources inefficently (only the 2nd
    step in the dv pipeline can make use of GPU-computing). As so, it
    is better to run the 1st/3rd step on a normal compute node and run
    the 2nd step on a GPU node. This rule is depreciated. Please see 
    the deepvariant_make_examples, deepvariant_call_variants, and
    deepvariant_postprocess_variants rules for the optimal way to
    run this tool.
    @Input:
        Duplicate marked, sorted BAM file (scatter)
    @Output:
        Single-sample gVCF file with called variants
    """
    input: 
        bam = join(workpath, "BAM", "{name}.sorted.bam"),
        bai = join(workpath, "BAM", "{name}.sorted.bam.bai"),
    output:
        gvcf = join(workpath, "deepvariant", "gVCFs", "{name}.g.vcf.gz"),
        vcf  = join(workpath, "deepvariant", "VCFs", "{name}.vcf.gz"),
    params: 
        rname  = "deepvar",
        genome = config['references']['GENOME'],
        tmpdir = tmpdir,
        # Building option for glnexus config, where:
        #  @WES = --model_type=WES
        #  @WGS = --model_type=WGS
        dv_model_type = lambda _: "WES" if run_wes else "WGS",
    message: "Running DeepVariant on '{input.bam}' input file"
    threads: int(allocated("threads", "deepvariant", cluster))
    container: config['images']['deepvariant']
    envmodules: config['tools']['deepvariant']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    run_deepvariant \\
        --model_type={params.dv_model_type} \\
        --ref={params.genome} \\
        --reads={input.bam} \\
        --output_gvcf={output.gvcf} \\
        --output_vcf={output.vcf} \\
        --num_shards={threads} \\
        --intermediate_results_dir=${{tmp}}
    """

# Depreciated somatic variant calling rule(s)
rule deepsomatic:
    """
    Data processing step to call somatic variants using deep neural 
    network in tumor-normal pairs. DeepSomatic is an extension of the
    deep learning-based variant caller DeepVariant that takes aligned
    reads (in BAM or CRAM format) from tumor and normal data, produces 
    pileup image tensors from them, classifies each tensor using a CNN,
    and  finally reports somatic variants in a standard VCF or gVCF file. 
    This rule runs all three steps in the deepsomatic pipeline as a one 
    step: i.e. make_examples, call_variants, and postprocess_variants.
    This is not optimal for large-scale projects as it will consume a lot
    of resources inefficently (only the 2nd step in the dv pipeline can
    make use of GPU-computing). As so, it is better to run the 1st/3rd 
    step on a normal compute node and run the 2nd step on a GPU node.
    @Input:
        Duplicate marked, sorted Tumor-Normal BAM file (scatter)
    @Output:
        Single-sample VCF file with called somatic variants
    """
    input: 
        tumor  = join(workpath, "BAM", "{name}.sorted.bam"),
        normal = get_normal_sorted_bam
    output:
        vcf  = join(workpath, "deepsomatic", "somatic", "{name}.deepsomatic.vcf"),
    params: 
        rname  = "deepsom",
        genome = config['references']['GENOME'],
        tmpdir = tmpdir,
        # Building option for deepsomatic config, where:
        #  @WGS = --model_type=WGS
        #  @WES = --model_type=WES  (may be added in future)
        dv_model_type = "WGS",
        # Get tumor and normal sample names 
        tumor  = '{name}',
        # Building option for the paired normal sorted bam
        normal_bam_option = lambda w: "--reads_normal={0}.sorted.bam".format(
            join(workpath, "BAM", tumor2normal[w.name])
        ) if tumor2normal[w.name] else "",
        # Building option for the normal sample name
        normal_name_option = lambda w: "--sample_name_normal={0}".format(
            tumor2normal[w.name]
        ) if tumor2normal[w.name] else "",
    threads: int(allocated("threads", "deepsomatic", cluster))
    container: config['images']['deepsomatic']
    envmodules: config['tools']['deepsomatic']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'du -sh "${{tmp}}"; rm -rf "${{tmp}}"' EXIT

    # Export OpenBLAS variable to
    # control the number of threads
    # in a thread pool. By setting
    # this variable to 1, work is
    # done in the thread that ran
    # the operation, rather than
    # disbatching the work to a
    # thread pool. If this option
    # is not provided, it can lead
    # to nested parallelism.
    # See this issue for more info:
    # https://github.com/google/deepsomatic/issues/28
    export OPENBLAS_NUM_THREADS=1

    # Run deepsomatic
    run_deepsomatic \\
        --model_type={params.dv_model_type} \\
        --ref={params.genome} \\
        --reads_tumor={input.tumor} {params.normal_bam_option} \\
        --sample_name_tumor={params.tumor} {params.normal_name_option} \\
        --output_vcf={output.vcf} \\
        --num_shards={threads} \\
        --intermediate_results_dir=${{tmp}} \\
        --vcf_stats_report=false
    """
