#!/usr/bin/env Rscript
library("sequenza")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("Must provide a seqz file")
} else {
    seqz_file <- args[1]
    if (! file.exists(seqz_file)) {
        stop(
            paste0(
                "Can't find this SEQZ output file: ",
                seqz_file
            )
        )
    }
}

if (length(args) > 1) {
    out_dir <- args[2]
} else {
    out_dir <- dirname(seqz_file)
}

if (length(args) > 2) {
    n_cores <- as.numeric(args[3])
} else {
    n_cores <- as.numeric(
        Sys.getenv("SLURM_CPUS_PER_TASK", 1)
    )
}
if (is.na(n_cores)) {
    n_cores = 1
}

if (length(args) > 3) {
    sampleid <- args[4]
} else {
    sampleid <- gsub(
        ".seqz.gz",
        "",
        basename(seqz_file)
    )
}

print(paste0("Using ", n_cores, " cores..."))
print("Extracting seqz data...")
date()
seqzdata <- sequenza.extract(
    seqz_file,
    min.reads = 30,
    min.reads.normal= 20
)

print("Fitting model...")
date()
# CP.example <- sequenza.fit(seqzdata, mc.cores = n_cores)
CP.example <- sequenza.fit(seqzdata)

# Sequenza.extract seems to fail if too few mutations
# num_mutations <- unlist(lapply(seqzdata$mutations, nrow))
# chrom_list <- names(num_mutations)[num_mutations >= 0]
# But it might actually be segments, idk?
# num_mutations <- unlist(lapply(seqzdata$segments, nrow))
# chrom_list <- names(num_mutations)[num_segments > 0]
chrom_list <- c(
    "chr1", "chr2", "chr3",
    "chr4", "chr5", "chr6",
    "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12",
    "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21",
    "chr22", "chrX"
)
# not_included <- setdiff(names(num_mutations), chrom_list)
# if (length(not_included) > 0) {
#    print("Excluding these chromosomes because of too few mutations...")
#    print(not_included)
# }
date()
print("Printing results...")
sequenza.results(
    sequenza.extract = seqzdata,
    cp.table = CP.example,
    sample.id = sampleid,
    out.dir = out_dir,
    chromosome.list = chrom_list
)

date()
print("Done")