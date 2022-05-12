# <code>genome-seek <b>cache</b></code>

## 1. About 
The `genome-seek` executable is composed of several inter-related sub commands. Please see `genome-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>genome-seek <b>cache</b></code> sub command in more detail. With minimal configuration, the **`cache`** sub command enables you to cache remote resources for the  genome-seek pipeline. Caching remote resources allows the pipeline to run in an offline mode. The cache sub command can also be used to pull our pre-built reference bundles onto a new cluster or target system.

The cache sub command creates local cache on the filesysytem for resources hosted on DockerHub or AWS S3. These resources are normally pulled onto the filesystem when the pipeline runs; however, due to network issues or DockerHub pull rate limits, it may make sense to pull the resources once so a shared cache can be created and re-used. It is worth noting that a singularity cache cannot normally be shared across users. Singularity strictly enforces that its cache is owned by the user. To get around this issue, the cache subcommand can be used to create local SIFs on the filesystem from images on DockerHub.

## 2. Synopsis

Coming Soon!  

<!-- ```text
$ ./genome-seek cache [-h] --sif-cache SIF_CACHE \
                     [--resource-bundle RESOURCE_BUNDLE] \
                     [--dry-run] 
```

The synopsis for each command shows its parameters and their usage. Optional parameters are shown in square brackets.

A user **must** provide a directory to cache remote Docker images via the `--sif-cache` argument. Once the cache has pipeline completed, the local sif cache can be passed to the `--sif-cache` option of the <code>genome-seek <b>build</b></code> and <code>genome-seek <b>run</b></code> subcomand. This enables the build and run pipeline to run in an offline mode.

Use you can always use the `-h` option for information on a specific command.

### 2.1 Required Arguments

`--sif-cache SIF_CACHE` 
 
> **Path where a local cache of SIFs will be stored..**  
> *type: string*
> 
> Any images defined in *config/containers/images.json* will be pulled into the local filesystem. The path provided to this option can be passed to the `--sif-cache` option of the <code>genome-seek <b>build</b></code> and <code>genome-seek <b>run</b></code> subcomand. This allows for running the build and run pipelines in an offline mode where no requests are made to exteexomel sources. This is useful for avoiding network issues or DockerHub pull rate limits. Please see genome-seek build and run for more information.
> 
> ***Example:*** `--sif-cache /data/$USER/cache`

### 2.2 Options

Each of the following arguments are optional and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

---  
  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

## 3. Example
```bash 
# Step 0.) Grab an interactive node (do not run on head node)
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
module purge
module load singularity snakemake

# Step 1.) Dry run cache to see what will be pulled
./genome-seek cache --sif-cache /scratch/$USER/cache \
                 --dry-run  

# Step 2.) Cache remote resources locally 
./genome-seek cache --sif-cache /scratch/$USER/cache  
```
-->