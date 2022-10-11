# <code>genome-seek <b>cache</b></code>

## 1. About 
The `genome-seek` executable is composed of several inter-related sub commands. Please see `genome-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>genome-seek <b>cache</b></code> sub command in more detail.

With minimal configuration, the **`cache`** sub command enables you to cache remote software containers from [Dockerhub](https://hub.docker.com/u/skchronicles). Caching remote software containers allows the pipeline to run in an offline mode where no requests are made. The cache sub command can also be used to pull our pre-built software container onto a new cluster or target system.

These containers are normally pulled onto the filesystem when the pipeline runs; however, due to network issues or DockerHub pull rate limits, it may make sense to pull the resources once so a shared cache can be created. It is worth noting that a singularity cache cannot normally be shared across users. Singularity strictly enforces that a cache is owned by the user. To get around this issue, the cache subcommand can be used to create local SIFs on the filesystem from images on DockerHub. The path of these locally cached SIFs can be passed to the run sub commands --sif-cache option.

Caching software containers is fast and easy! In its most basic form, <code>genome-seek <b>cache</b></code> only has *one required input*.

## 2. Synopsis
```text
$ ./genome-seek cache [--help] [--dry-run] \
       --sif-cache SIF_CACHE
```

The synopsis for each command shows its parameters and their usage. Optional parameters are shown in square brackets.

A user **must** provide a directory to cache remote Docker images via the `--sif-cache` argument. Once the cache has pipeline completed, the local sif cache can be passed to the `--sif-cache` option of the <code>genome-seek <b>run</b></code> subcomand. This enables the pipeline to run in an offline mode.

Use you can always use the `-h` option for information on a specific command.

### 2.1 Required Arguments

`--sif-cache SIF_CACHE` 
 
> **Path where a local cache of SIFs will be stored.**  
> *type: path*
> 
> Any images defined in *config/containers.json* will be pulled into the local filesystem. The path provided to this option can be passed to the `--sif-cache` option of the <code>genome-seek <b>run</b></code> subcomand. This allows for running the build and run pipelines in an offline mode where no requests are made to external sources. This is useful for avoiding network issues or DockerHub pull rate limits. Please see genome-seek run for more information.
> 
> ***Example:*** `--sif-cache /data/$USER/cache`

### 2.2 Options

Each of the following arguments are optional and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

---  
  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Only displays what software container will be cached locally. Does not execute anything!
>
> ***Example:*** `--dry-run`

## 3. Example
```bash 
# Step 0.) Grab an interactive node (do not run on head node)
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
module purge
module load singularity snakemake

# Step 1.) Dry run to see what will be pulled
./genome-seek cache --sif-cache /data/$USER/cache \
                 --dry-run  

# Step 2.) Cache remote resources locally.
# This command will NOT automatically submit
# a job to the cluster. As so, we recommend 
# submitting this next command to the cluster
# as a job. Download speeds will vary so it 
# is best to set the wall time a few hours. 
./genome-seek cache --sif-cache /data/$USER/cache  
```
