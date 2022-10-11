#!/usr/bin/env bash

set -euo pipefail

usage="Usage: $0 <path_to_resource_bundle_directory> <output_archive.tar.gz>"

function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; err "$usage"; exit 1; }
function abspath() { readlink -f "$1"; }
function bundle() { tar -hczvf "$1" "$2"; }


function check() {
    # Checks command-line usage for required positional arguments.
    # @INPUT $1 = Path on filesytem to create a tarball
    # @INPUT $2 = Name of of output tarballl to create
    # @CALLS fatal() with incorrect usage

    die=false
    # Path to archive
    if [ -z "${1:-}" ]; then
        die=true
        err "Error: Failed to provide directory to archive."
    fi
    # Output tarball name
    if [ -z "${2:-}" ]; then
        die=true
        err "Error: Failed to output file name for archive."
    fi
    if $die; then
        fatal "Fatal: Please try again after providing the required arguments!"
    fi
}


function chunk() {
    # Splits resulting tarball into N 4GB chunks
    # @INPUT $1 = Name of of output tarballl
    # @CALLS fatal() if provided a non-supported archive

    # Strip archive file extension,
    # common tarball extensions: .tar.gz or .tgz
    prefix=''
    if [ -f $1 ] ; then
        case $1 in
            *.tar.gz) prefix="${1%.tar.gz}" ;;
            *.tgz)    prefix="${1%.tgz}"    ;;
        esac
    else
        fatal "'$1' is not supported file type"
    fi

    # Spilt file into N 4GB chunk files
    split --numeric-suffixes=1 -b 4G "$1" "${prefix}_chunk-"

    # Calculate MD5 of all the chunks
    md5sum "${prefix}_chunk-"* > "${prefix}_chunks.md5"
}


function main() {
    # Checks for required positional
    # command line arguments
    check "${1:-}" "${2:-}" 

    # Converts any relative paths to
    # absolute paths, creates uninit
    # output directories as needed,
    # runs tar command in parent dir
    # of the provided resource bundle
    # path.
    archive_dir=$(abspath "$1")
    archive=$(basename "${archive_dir%/}")
    parent_dir=$(dirname "$archive_dir")
    output_dir=$(dirname "$2")
    output_dir=$(abspath "$output_dir")

    mkdir -p "$output_dir"
    cd "$parent_dir"

    # Create archive as a tarball
    echo "Creating tarball... $2"
    bundle "$2" "$archive"

    # Splitting tarball into N
    # chunks for fast parallel
    # downloading of large
    # resource bundles
    echo "Chunking tarball... $2"
    chunk "$2"
}


main "$@"
