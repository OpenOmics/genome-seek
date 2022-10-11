#!/usr/bin/env bash
set -euo pipefail

function usage() { cat << EOF
cache.sh: Wrapper script for caching remote software containers.
USAGE:
  cache.sh <MODE> [OPTIONS] -s <SIF_CACHE> -i <IMAGE_URIs>

SYNOPSIS:
  This main process dictates how subsequent software containers are 
pulled onto the cluster's local filesystem. cache.sh pulls containers 
from Dockerhub locally. Docker images are converted on the fly into 
singularity image format.
  The main entry point of the pipeline calls this wrapper script. 
As so, this script can be used to manually by-pass the pipeline for 
a previously failed cache or for the purpose of debugging.

Required Positional Argument:
  [1] MODE  [Type: Str] Defines the mode of execution. More methods 
                        can be added later. Valid mode options include: 
                         a) local: uses singularity and local compute.

Required Arguments:
  -s, --sif-cache   [Type: Path]  Path to output directory to cache 
                                  software containers, i.e. SIFs.
  
  -i, --image-uris  [Type: Str]   Image(s) to pull from Dockerhub. 
                                  Multiple images are seperated by 
                                  a comma.
OPTIONS:
  -t, --tmp-dir     [Type: Path]  Path to singularity temp directory. 
                                  Singularity uses this directory when 
                                  images are pulled from DockerHub and 
                                  coverted into SIFs. If not provided, 
                                  the location to the temp dir will 
                                  default to the following location:
                                  "/tmp/$USER/SIFs/.singularity/".

  -h, --help        [Type: Bool]  Displays usage and help information.

Example:
  $ cache.sh local \\
      -s $PWD/$USER/SIFs \\
      -t $PWD/$USER/SIFs/.singularity \\
      -i 'docker://nciccbr/ccbr_arriba_2.0.0:v0.0.1,docker://nciccbr/ccbr_rna:v0.0.1'
Version:
  0.2.0
EOF
}


# Functions
function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; usage; exit 1; }
function abspath() { readlink -e "$1"; }
function parser() {
  # Adds parsed command-line args to GLOBAL $Arguments associative array
  # + KEYS = short_cli_flag ("j", "o", ...)
  # + VALUES = parsed_user_value ("MasterJobName" "/scratch/hg38", ...)
  # @INPUT "$@" = user command-line arguments
  # @CALLS check() to see if the user provided all the required arguments

  while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
      -h  | --help) usage && exit 0;;
      -s  | --sif-cache)  provided "$key" "${2:-}"; Arguments["s"]="$2"; shift; shift;;
      -i  | --image-uris) provided "$key" "${2:-}"; Arguments["i"]="$2"; shift; shift;;
      -t  | --tmp-dir)    provided "$key" "${2:-}"; Arguments["t"]="$2"; shift; shift;;
      -*  | --*) err "Error: Failed to parse unsupported argument: '${key}'."; usage && exit 1;;
      *) err "Error: Failed to parse unrecognized argument: '${key}'. Do any of your inputs have spaces?"; usage && exit 1;;
    esac
  done

  # Check for required args
  check
}


function provided() {
  # Checks to see if the argument's value exists
  # @INPUT $1 = name of user provided argument
  # @INPUT $2 = value of user provided argument
  # @CALLS fatal() if value is empty string or NULL

  if [[ -z "${2:-}" ]]; then
     fatal "Fatal: Failed to provide value to '${1}'!";
  fi
}


function check(){
  # Checks to see if user provided required arguments
  # @INPUTS $Arguments = Global Associative Array
  # @CALLS fatal() if user did NOT provide all the $required args

  # List of required arguments
  local required=("s" "i")
  #echo -e "Provided Required Inputs"
  for arg in "${required[@]}"; do
    value=${Arguments[${arg}]:-}
    if [[ -z "${value}" ]]; then
      fatal "Failed to provide all required args.. missing ${arg}"
    fi
  done
}


function retry() {
  # Tries to run a cmd 5 times before failing 
  # If a command is successful, it will break out of attempt loop
  # Failed attempts are padding with the following exponential 
  # back-off strategy {4, 16, 64, 256, 1024} in seconds
  # @INPUTS "$@"" = cmd to run 
  # @CALLS fatal() if command cannot be run in 5 attempts
  local n=1
  local max=5
  local attempt=true # flag for while loop
  while $attempt; do
    # Attempt command and break if successful 
    "$@" && attempt=false || {
      # Try again up to 5 times 
      if [[ $n -le $max ]]; then
        err "Command failed: $@"
        delay=$(( 4**$n ))
        err "Attempt: ${n}/${max}. Trying again in ${delay} seconds!\n"
        sleep $delay;
        ((n++))
      else
        fatal "Fatal: the command has failed after max attempts!"
      fi
    }
  done
}


function _pull(){
  # Caches a remote image from DockerHub
  # INPUT $1 = Snakemake Mode of execution
  # INPUT $2 = Cache output directory
  # INPUT $3 = Singularity temp directory
  # INPUT $4 = Images to pull from DockerHub

  # Check if singularity in $PATH
  # If not, try to module load singularity as a last resort
  command -V singularity &> /dev/null || { 
    command -V module &> /dev/null && 
    module purge && module load singularity
  } || fatal "Fail to find or load 'singularity', not installed on target system."

  # Execution method, currently pulls 
  # from local compute, in the future
  # options can be added to submit a 
  # to different job schedulers, like
  # PBS or SLURM, etc
  executor=${1}

  # Goto Pipeline Ouput directory
  # Create a local singularity cache in output directory
  # cache can be re-used instead of re-pulling from DockerHub everytime
  cd "$2" && export SINGULARITY_CACHEDIR="${3}"

  # unsetting XDG_RUNTIME_DIR to avoid some unsighly but harmless warnings
  unset XDG_RUNTIME_DIR

  # Run the workflow with specified executor
  case "$executor" in
    local)
          # Create directory for logfiles
          for image in ${4//,/$'\t'}; do 
            # Try to pull image from URI with 5 max attempt
            echo "Singularity pulling ${image}"
            retry singularity pull -F ${image}
          done
        ;;
      *)  echo "${executor} is not available." && \
          fatal "Failed to provide valid execution backend: ${executor}. Please use local."
        ;;
    esac
}


function main(){
  # Parses args and pulls remote resources
  # @INPUT "$@" = command-line arguments
  # @CALLS pull()

  if [ $# -eq 0 ]; then usage; exit 1; fi

  # Associative array to store parsed args
  declare -Ag Arguments

  # Positional Argument for Executor
  case $1 in
    local) Arguments["e"]="$1";;
    -h    | --help | help) usage && exit 0;;
    -*    | --*) err "Error: Failed to provide required positional argument: <local>."; usage && exit 1;;
    *) err "Error: Failed to provide valid positional argument. '${1}' is not supported. Valid option(s) are local"; usage && exit 1;;
  esac

  # Parses remaining user provided command-line arguments
  parser "${@:2}" # Remove first item of list
  mkdir -p "${Arguments[s]}"
  cache=$(abspath "${Arguments[s]}")
  dockers="${Arguments[i]}"

  # Setting defaults for non-required arguments
  tmp="${Arguments[t]:-/tmp/$USER/SIFs/.singularity/}"
  # Sanitize single quoted input
  tmp=$(echo "$tmp" | sed "s/\${SLURM_JOB_ID}/${SLURM_JOB_ID}/g" | sed "s/\$SLURM_JOB_ID/${SLURM_JOB_ID}/g")
  mkdir -p "$tmp"
  tmp=$(abspath "${tmp}")

  # Print cli options prior to running
  echo -e "cache.sh \t$(date)"
  echo -e "Running with the following parameters:"
  for key in "${!Arguments[@]}"; do echo -e "\t${key}\t${Arguments["$key"]}"; done

  # Pull software containers into SIF cache 
  # Cache remote image from DockerHub
  # INPUT $1 = Snakemake Mode of execution
  # INPUT $2 = Cache output directory
  # INPUT $3 = Singularity temp directory
  # INPUT $4 = Images to pull from DockerHub
  _pull "${Arguments[e]}" "$cache" "$tmp" "$dockers"

}


# Main: check usage, parse args, and run pipeline
main "$@"