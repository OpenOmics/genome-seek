# <code>genome-seek <b>run</b></code>

## 1. About 
The `genome-seek` executable is composed of several inter-related sub commands. Please see `genome-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>genome-seek <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running genome-seek pipeline. 

Setting up the genome-seek pipeline is fast and easy! In its most basic form, <code>genome-seek <b>run</b></code> only has *two required inputs*.

## 2. Synopsis
```text
$ genome-seek run [--help] \
      [--mode {slurm,uge,local}] [--job-name JOB_NAME] [--batch-id BATCH_ID] \
      [--call-cnv] [--call-sv] [--call-hla] [--call-somatic] [--gatk-germline] \
      [--open-cravat] [--oc-annotators OC_ANNOTATORS] [--oc-modules OC_MODULES] \
      [--pairs PAIRS] [--pon PANEL_OF_NORMALS] [--wes-mode] [--wes-bed WES_BED] \
      [--skip-qc] [--tmp-dir TMP_DIR] [--silent] [--sif-cache SIF_CACHE] \
      [--singularity-cache SINGULARITY_CACHE] \
      [--resource-bundle RESOURCE_BUNDLE] \
      [--dry-run] [--threads THREADS] \
      --input INPUT [INPUT ...] \
      --output OUTPUT
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument and an output directory to store results via `--output` argument.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ or BAM file(s).**  
> *type: file(s)*  
> 
> One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should seperated by a space. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
> 
> ***Example:*** `--input .tests/*.R?.fastq.gz`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/genome-seek_out`

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. 

  `--call-cnv`            
> **Call copy number variation (CNV).**  
> *type: boolean flag*
> 
> Runs additional steps to predict copy number variation.
>
> ***Example:*** `--call-cnv`

---  
  `--call-sv`            
> **Call structural variation (SV).**  
> *type: boolean flag*
> 
> Runs additional steps to predict structural variation.
>
> ***Example:*** `--call-sv`


---  
  `--call-hla`            
> **Call HLA types.**  
> *type: boolean flag*
> 
> Runs additional steps to call HLA types.
>
> ***Example:*** `--call-hla`

---  
  `--call-somatic`            
> **Call somatic variants.**  
> *type: boolean flag*
> 
> Runs additional steps to call somatic variants. By default when this option is provided, the pipeline will perform somatic variant calling for each sample in a tumor-only mode; however, if a tumor-normal pairs file is provided via the `--pairs` option then the pipeline will call somatic variants using its matched normal. Please see the **Somatic options** section below for more information about the somatic variant calling pipeline and any additional options.
>
> ***Example:*** `--call-somatic`

---  
  `--gatk-germline`            
> **Call short germline variants using GATK4 best practices.**  
> *type: boolean flag*
> 
> Runs additional steps to call short (SNPs/INDELs) germline variants using GATK4. By default, the pipeline will call germline variants using deepvariant. If this option is provided, the pipeline will also call germline variants using GATK4's set of recommended best practices for calling short germline variants.
>
> ***Example:*** `--gatk-germline`

---  
  `--open-cravat`            
> **Run OpenCRAVAT to annotate variants.**  
> *type: boolean flag*
> 
> Runs additional steps to annotate variants with OpenCRAVAT. See the Annotation Options section below for more information about what modules are included by default and how to include more modules. 
>
> ***Example:*** `--open-cravat`

---  
  `--skip-qc`            
> **Skips over quality control steps.**  
> *type: boolean flag*
> 
> When this option is provided the pipeline will not run any of its QC related steps. This means that only data processing steps will run like trimming, alignment, and variant  calling. It is worth noting that we do not recommend skipping over QC; however, in certain scenarios it may make sense. This option is useful for testing changes to data processing steps or when evaluating the overall accuracy and precision of the pipeline.
>
> ***Example:*** `--skip-qc`

---  
  `--wes-mode`            
> **Runs the whole exome pipeline.**  
> *type: boolean flag*
> 
> By default, the whole genome sequencing (WGS) pipeline is run. This option allows a user to process and analyze whole exome sequencing data. Please note when this mode is enabled, a sub-set of the WGS rules will run. Please see the option below for more information about providing a custom exome targets/capture-kit BED file.
>
> ***Example:*** `--wes-mode`

---  
  `--wes-bed WES_BED`            
> **Path to exome targets/capture-kit BED file.**  
> *type: BED file*
>
> This file can be obtained from the manufacturer of the target capture kit that was used. By default, a set of BED files generated from GENCODE's exon annotation for protein coding gene's exon is used. Please note: This BED file should contain at least 6 columns.
>
> ***Example:*** `--wes-bed Agilent_SS_AllExons_V7_Regions.bed`

---  
  `--batch-id BATCH_ID`            
> **Unique identifer for a batch of samples.**  
> *type: string*
> *default: MD5 of the sorted file of input file names*
> 
> A batch identifer should be a string containing no spaces. If a batch identifer is not provided, one will be generated by taking the MD5 checksum of the sorted list of input file names. The batch identifer is added to some output file names to ensure files are not over-written when the pipeline is run with a different set of samples. This allows the pipeline to be run with extra samples without overriding some of results of a previous run.
>
> ***Example:*** `--batch-id WGS_2022-04-19`

### 2.3 Somatic options

Each of the following arguments are optional, and do not need to be provided. 

  `--pairs PAIRS`
> **Path to a tumor-normal pairs file.**
> *type: TSV file*
> *default: None*
> 
> The tumor-normal pairs file is used to pair a tumor sample with its match normal sample. This file should only be provided when calling somatic variants via the `--call-somatic` option. Please see the option above for more info about calling somatic mutations. By default, the `--call-somatic` option will call somatic variants for each sample in a *tumor-only* mode. 
>
> This tab delimited file contains two columns with the names of tumor and normal pairs, one per line. The header of the file needs to be `Tumor` for the tumor column and `Normal` for the normal column. The base name of each sample should be listed in the pairs file. The base name of a given sample can be determined by removing the extension from the sample's R1 FastQ file, e.g. `.R1.fastq.gz`.
>
> **Contents of example pairs file:**
> ```
> Tumor    Normal
> Sample4_CRL1622_S31  Sample10_ARK1_S37
> Sample4_CRL1622_S31  Sample11_ACI_158_S38
> ```
> 
> ***Example:*** `--pairs .tests/pairs.tsv`

---  
  `--pon PANEL_OF_NORMALS`
> **Path to a panel of normals file.**
> *type: VCF.gz file*
> *default: None*
> 
> A VCF file of containing sites observed in normal tissue. Normal in this context refers to samples derived from healthy tissue that is NOT believed to have any somatic alterations. By default, the pipeline will use a PON included with its resource bundle. You can provide your own PON with this option. The PON should be gzipped **and** there should be a tabix index for the PON in the same directory.
> 
> ***Example:*** `--pon 1000g_pon.hg38.vcf.gz`

### 2.4 Anotation options

Each of the following arguments are optional, and do not need to be provided. 

#### 2.4.1 OpenCRAVAT

  `--oc-annotators OC_ANNOTATORS`            
> **List of OpenCRAVAT annotators to use.**  
> *type: string*
> *default: config/oc_annotators.cfg*
> 
> Please note that one or more annotators can be provided. Multiple annotators should be seperate with a space. By default, the pipeline will annotate variants using  modules in the resource bundle; however, a custom module installation path can be defined using the `--oc-modules` option. Annotators listed in `config/oc_annotators.cfg` can be skipped over my commenting them out.
>
> ***Example:*** `--oc-annotators dann dann_coding`

---  
  `--oc-modules OC_MODULES`            
> **Sets the path to OpenCRAVAT's modules directory.**  
> *type: path*
> 
> This option overrides the default path any installed modules for OpenCRAVAT. Use this option if you have your own custom installation of OpenCRAVAT modules. An additional annotators installed the provided path can be add to `config/oc_annotators.cfg`. To find where
OpenCRAVAT has installed its modules, please run `oc config system`. Also, please take caution to ensure the annotators used at run time via the `--oc-annotators` option exist in the provided path.
>
> ***Example:*** `--oc-modules /data/$USER/CRAVAT/modules`

### 2.5 Orchestration options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,uge,local}`  
> **Execution Method.**  
> *type: enum*  
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local. 
> 
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running genome-seek in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***uge***    
> The ueg execution method will submit jobs to the [Univa Grid Engine](https://en.wikipedia.org/wiki/Univa_Grid_Engine). This method will submit jobs to a cluster using qsub. Please set the mode to uge when running the pipeline on LOCUS.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode. 
> 
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*
> *default: pl:genome-seek*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:genome-seek".
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `genome-seek cache` subcommand can be used to create a local SIF cache. Please see `genome-seek cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running genome-seek with this option when ever possible.
> 
> ***Example:*** `--sif-cache /data/$USER/SIFs`

---  
  `--resource-bundle RESOURCE_BUNDLE`
> **Path to a local resource bundle.**  
> *type: path*  
>
> This is a path to a local resource bundle containing all of the pipeline's reference files. Please only provide this option if you are running the pipeline outside of Biowulf. If you are running the pipeline on Biowulf, the pipeline will automatically resolve the correct path to any references files. The resource bundle contains the set of required reference files for processing any data, and it is extrememly large. Several terabytes of storage space will be needed to download it. As a result, we always recommend using the resource bundle provided by the pipeline. If you are unsure if the resource bundle is installed on your system, please contact us or your system administrator.
> 
> ***Example:*** `--resource-bundle /path/to/refs/genome-seek`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`

---  
  `--tmp-dir TMP_DIR`   
> **Max number of threads for each process.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOBID`*
> 
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
> 
> ***Example:*** `--tmp-dir /scratch/$USER/`

### 2.6 Miscellaneous options  
Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example
```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./genome-seek run --input .tests/*.R?.fastq.gz \
        --output results/ \
        --call-cnv --call-sv \
        --call-hla --open-cravat \
        --gatk-germline \
        --call-somatic \
        --mode slurm \
        --dry-run

# Step 2B.) Run the genome-seek pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./genome-seek run --input .tests/*.R?.fastq.gz \
        --output results/ \
        --call-cnv --call-sv \
        --call-hla --open-cravat \
        --gatk-germline \
        --call-somatic \
        --mode slurm
```