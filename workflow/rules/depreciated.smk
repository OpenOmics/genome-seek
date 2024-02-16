# Depreciated rules that may still be useful for some projects

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