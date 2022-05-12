# <code>genome-seek <b>unlock</b></code>

## 1. About 
The `genome-seek` executable is composed of several inter-related sub commands. Please see `genome-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>genome-seek <b>unlock</b></code> sub command in more detail. With minimal configuration, the **`unlock`** sub command enables you to unlock a pipeline output directory. 

If the pipeline fails ungracefully, it maybe required to unlock the working directory before proceeding again. Snakemake will inform a user when it maybe necessary to unlock a working directory with an error message stating: `Error: Directory cannot be locked`. 

Please verify that the pipeline is not running before running this command. If the pipeline is currently running, the workflow manager will report the working directory is locked. The is the default behavior of snakemake, and it is normal. Do NOT run this command if the pipeline is still running! Please kill the master job and it's child jobs prior to running this command.

Unlocking genome-seek pipeline output directory is fast and easy! In its most basic form, <code>genome-seek <b>unlock</b></code> only has *one required input*.

## 2. Synopsis
```text
$ ./genome-seek unlock [-h] --output OUTPUT
```

The synopsis for this command shows its parameters and their usage. Optional parameters are shown in square brackets.

A user **must** provide an output directory to unlock via `--output` argument. After running the unlock sub command, you can resume the build or run pipeline from where it left off by re-running it. 

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required Arguments  

  `--output OUTPUT` 
> **Output directory to unlock.**  
> *type: path*
> 
> Path to a previous run's output directory. This will remove a lock on the working  directory. Please verify that the pipeline is not running before running this command.  
> ***Example:*** `--output /data/$USER/genome-seek_out`

### 2.2 Options

Each of the following arguments are optional and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`


## 3. Example
```bash 
# Step 0.) Grab an interactive node (do not run on head node)
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
module purge
module load singularity snakemake

# Step 1.) Unlock a pipeline output directory
genome-seek unlock --output /data/$USER/output
```