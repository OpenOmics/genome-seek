#!/usr/bin/env bash
set -euo pipefail

# Usage and help message
function usage() {
  cat << EOF
Usage: $(basename "$0") [-h] \\
           [-c THREADS] \\
           [-b CAPTURE_BED] \\
            -s SAMPLE_ID \\
            -t TUMOR_BAM \\
            -n NORMAL_BAM \\
            -r REF_FASTA \\
            -g GC_WIGGLE \\
            -e SPECIES \\
            -b CAPTURE

Required:
  -s SAMPLE_ID    Sample ID
  -t TUMOR_BAM    Tumor BAM file
  -n NORMAL_BAM   Normal BAM file
  -r REF_FASTA    FASTA file 
  -g GC_WIGGLE	  GC Wiggle file
  -e SPECIES      Human or Mouse

Optional:
  -c THREADS      The number of threads (default 1)
  -b CAPTURE	    WES Bed file of capture (default: none)
  -h              Display help message and exit
EOF
}

# Functions
function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; usage; exit 1; }
function timestamp() { date +"%Y-%m-%d_%H-%M-%S"; }
function require(){
  # Requires an executable is in $PATH, 
  # as a last resort it will attempt to load
  # the executable or dependency as a module
  # @INPUT $@ = List of executables to check
  for exe in "${@}"; do
    # Check if executable is in $PATH
    command -V "${exe}" &> /dev/null && continue;
    # Try to load exe as lua module
    module load "${exe}" &> /dev/null || \
      fatal "Failed to find or load '${exe}', not installed on target system."
  done
}


# Parse command-line options
# Set defaults for non-required options
num_threads=1      # default: number of threads to use
new_bait="none"    # default: none (i.e. run in WGS-mode)
while getopts s:t:n:r:c:g:e:b:h OPT; do
  case $OPT in
    s ) sample_id=$OPTARG;;
    t ) tumor_bam=$OPTARG;;
    n ) normal_bam=$OPTARG;;
    r ) reference_fasta=$OPTARG;;
    c ) num_threads=$OPTARG;;
    g ) gc_wiggle=$OPTARG;;
    e ) species=$OPTARG;;
    b ) new_bait=$OPTARG;;
    h ) usage && exit 0;;
    ? ) usage && exit 1;;
  esac
done

# Sanity check: was anything was provided?!
{ [ -z "${1:-}" ] ; } \
  && fatal "Error: Did not provide any of the required options!"

# Check if all required options 
# where provided at runtime
declare -A required_options
required_options=(
  ["s"]="${sample_id:-}" 
  ["t"]="${tumor_bam:-}" 
  ["n"]="${normal_bam:-}" 
  ["r"]="${reference_fasta:-}" 
  ["g"]="${gc_wiggle:-}" 
  ["e"]="${species:-}" 
)
echo "Running with provided required options:"
for k in "${!required_options[@]}"; do
  v="${required_options[$k]:-}"
  echo "  -$k $v"
  if [[ -z "${v}" ]]; then
    fatal "Error: Failed to provide all required options... missing -${k} OPTION"
  fi
done


# Check for software dependencies,
# as last resort try to module load 
# the specified tool or dependency
require bedtools samtools sequenza-utils


# Set sequenza-utils run options
if [ "$species" == "Human" ]; then
    echo "Human as reference so chromosome 1 to 22 & X will be analyzed";
    chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
else
    echo "Mouse as reference so chromosome 1 to 19 & X will be analyzed";
    chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"
fi


echo "[$(timestamp)] Running sequenza-utils bam2seqz with the follow bam files: ${normal_bam}  ${tumor_bam}"; 
sequenza-utils bam2seqz \
    -n "${normal_bam}" \
    -t "${tumor_bam}" \
    --fasta "${reference_fasta}" \
    -gc "${gc_wiggle}" \
    -o "sequenza_out/${sample_id}/${sample_id}.seqz.gz" \
    -C ${chromosomes} \
    --parallel ${num_threads} \
    -N 40


# Merge the output across all chromosomes,
# not sure why we are doing it this way but
# didn't want to change the original script
# too much, it gets the job done, so let's
# keep it
{
    for chr in $chromosomes; do
      if [[ "$chr" = "chr1" ]]; then
        zcat "sequenza_out/${sample_id}/${sample_id}_chr1.seqz.gz" \
        > "sequenza_out/${sample_id}/${sample_id}.combine.chr1.seqz"
      else
        zcat "sequenza_out/${sample_id}/${sample_id}_${chr}.seqz.gz" \
          | tail -n +2 \
        >> "sequenza_out/${sample_id}/${sample_id}.combine.chrall.seqz"
      fi
   done
}
# Merge and gzip, probably could add
# a clean up step here as well
cat "sequenza_out/${sample_id}/${sample_id}.combine.chr1.seqz" "sequenza_out/${sample_id}/${sample_id}.combine.chrall.seqz" \
  | gzip -c  \
> "sequenza_out/${sample_id}/${sample_id}.combine.seqz.gz"


if [ "$new_bait" != "none" ]; then
    # WES capture kit provided via
    # the -b CAPTURE option, default
    # is none or run in WGS-mode
    zgrep -v "chromosome" "sequenza_out/${sample_id}/${sample_id}.combine.seqz.gz" \
      | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' \
      | bedtools intersect -a stdin -b "$new_bait" \
      | cut -f-2,4- \
    > "sequenza_out/${sample_id}/${sample_id}.combine.seqz.body"
    gunzip -c  "sequenza_out/${sample_id}/${sample_id}.combine.seqz.gz" \
      | head -n 1 \
    > "sequenza_out/${sample_id}/head"
    cat "sequenza_out/${sample_id}/head" "sequenza_out/${sample_id}/${sample_id}.combine.seqz.body" \
      | gzip -c \
    > "sequenza_out/${sample_id}/${sample_id}.combine.ontarget.seqz.gz"

    echo "[$(timestamp)] Running sequenza-utils seqz_binning with the follow input file: sequenza_out/${sample_id}/${sample_id}.combine.ontarget.seqz.gz"; 
    sequenza-utils seqz_binning \
      --seqz  "sequenza_out/${sample_id}/${sample_id}.combine.ontarget.seqz.gz" \
      -w 1000 \
      -o "sequenza_out/${sample_id}/${sample_id}.1kb.seqz.gz"
else
    # Run in WGS-like-mode
    echo "[$(timestamp)] Running sequenza-utils seqz_binning with the follow input file: sequenza_out/${sample_id}/${sample_id}.combine.seqz.gz"; 
    sequenza-utils seqz_binning \
      --seqz "sequenza_out/${sample_id}/${sample_id}.combine.seqz.gz" \
      -w 1000 \
      -o "sequenza_out/${sample_id}/${sample_id}.1kb.seqz.gz"
fi


# Clean tmp/intermediate files
rm -f "sequenza_out/${sample_id}/${sample_id}.combine.chr1.seqz"
rm -f "sequenza_out/${sample_id}/${sample_id}.combine.chrall.seqz"