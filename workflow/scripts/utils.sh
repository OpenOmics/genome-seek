#!/usr/bin/env bash

# Functions
function err() { cat <<< "$@" 1>&2; }
function fatal() { err "$@"; exit 1; }
function abspath() { readlink -e "$1"; }

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
        err "Attempt: ${n}/${max}. Trying again in ${delay} seconds!"
        sleep $delay;
        ((n++))
      else
        fatal "Fatal: the command has failed after max attempts!"
      fi
    }
  done
}


export -f err
export -f fatal
export -f retry
