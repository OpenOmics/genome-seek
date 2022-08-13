#!/bin/bash
# properties = {properties}

# Taking advantage of sbatch's ability to send 
# signals to jobs. This could prevent the master 
# job from running indefinitely when a rule's job 
# reaches it's walltime (the process is killed 
# ungracefully). This ungracefully exit causes the 
# rule.job.failed file to never be touched, which 
# leads to a condition that is never satified.
# To send a custom signal, like SIGUSR1, will use 
# sbatch's --signal option. We don't want to provide 
# Example sends SIGUSR1 two minutes before walltime:
# sbatch --time 00:03:00 --signal=B:USR1@120 script.sh

# Functions
function abspath() { readlink -f "$1"; }
function err() { cat <<< "$@" 1>&2; }
function fail(){
    # Ensures a process interrupted
    # with a signal exits gracefully
    # @INPUT $1 = Name of job failed file to touch
    # @INPUT $2 = Process ID to be terminated
    err "Terminating process $pid slightly before walltime to ensure graceful failure"
    err "slurmstepd error: JOB $SLURM_JOBID CANCELLED due to DUE TO TIME LIMIT ***"
    touch "$1"
    kill -15 "$2"
    exit 1
}


function run(){
    {exec_job}
}



# Get directory where script resides,
# the finished and failed file to track
# job status are in this directory
if [ -n "$SLURM_JOB_ID" ]; then
    # Running $0 from inside script will resolve to
    # /var/spool/slurm/slurmd/job${SLURM_JOBID}/jobscript.sh
    script="$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')"
else
    script="$0"
fi
snakemake_tmpdir="$(abspath "$(dirname  "$script")")"
# Create the name of the failed job file
# based on the name of this script. As an
# example, a job script  with the filename
# snakejob.muse.184.sh  will create a job
# tracking file called 184.jobfinished OR 
# 184.jobfailed depending on whether the 
# job was successful or failed.
failed_name="$(basename "$script" | awk -F '.' '{print $3".jobfailed"}')"
failed_file="${snakemake_tmpdir}/${failed_name}"
# Setup a trap to catch a user signal, SIGUSR1,
# that touches the failed job file prior to the
# job getting cancelled. If the process tracking
# job status is killed forecefully, i.e. with
# kill -9, the master job of the pipeline will
# hang indefinitely. This file-based approach to
# tracking jobs is reduces any strain related to
# querying the job scheduler.
run &
pid="$!"
trap "fail '$failed_name' '$pid'" USR1
wait $pid
