#!/usr/bin/env Rscript

library("maftools")

# Parse command-line args
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILE1 <- args[2]
FILE2 <- args[3]
FILE3 <- args[4]

setwd(DIR)
mymaf <- read.maf(FILE1)

# Read list of blacklist genes
flags <- c(
    "TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B", "DNAH17",
    "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK", "HMCN1",
    "USH2A", "DNAH11", "MACF1", "MUC17", "DNAH5", "GPR98", "FAT1", "PKD1",
    "MDN1", "RNF213", "RYR1", "DNAH2", "DNAH3", "DNAH8", "DNAH1", "DNAH9",
    "ABCA13", "SRRM2", "CUBN", "SPTBN5", "PKHD1", "LRP2", "FBN3", "CDH23",
    "DNAH10", "FAT4", "RYR3", "PKHD1L1", "FAT2", "CSMD1", "PCNT", "COL6A3",
    "FRAS1", "FCGBP", "RYR2", "HYDIN", "XIRP2", "LAMA1"
)

samplesum <- getSampleSummary(mymaf)
genesum <- getGeneSummary(mymaf)

write.table(
    samplesum, file = "sample_summary.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
    genesum, file = "gene_summary.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
)

# Create plots
print('Creating TCGA comparison plot...')
tryCatch(
    {   # Try to create plot
        pdf('cohort_tcga_comparison.pdf')
        tcgaCompare(maf = mymaf, cohortName = "project_variants")
    }, error = function(e) {
        # Catch error, usually due to low mutational burden
        print("Error: Failed to create TCGA comparison plot!")
        print(e)
        plot.new()
        text(.5, .5, "Error: Failed to create TCGA comparison plot!")
    }, finally = {
        dev.off()
    }
)

print('Creating genes by Vaf plot...')
tryCatch(
    {   # Try to create plot
        pdf('cohort_genes_by_VAF.pdf')
        plotVaf(mymaf, showN = TRUE, top = 20)
    }, error = function(e) {
        # Catch error, usually due to low mutational burden
        print("Error: Failed to create Vaf plot!")
        print(e)
        plot.new()
        text(.5, .5, "Error: Failed to create Vaf plot!")
    }, finally = {
        # Close Plotting device
        dev.off()
    }
)


print('Creating MAF summary plot...')
tryCatch(
    {   # Try to create plot
        pdf(FILE2)
        plotmafSummary(mymaf)
    }, error = function(e) {
        # Catch error, usually due to low mutational burden
        print("Error: Failed to create MAF summary plot!")
        print(e)
        plot.new()
        text(.5, .5, "Error: Failed to create MAF summary plot!")
    }, finally = {
        # Close Plotting device
        dev.off()
    }
)


print('Creating Oncoplot...')
tryCatch(
    {   # Try to create plot
        pdf(FILE3)
        oncoplot(
            mymaf, writeMatrix = TRUE, showTumorSampleBarcodes = TRUE,
            genesToIgnore = flags, removeNonMutated = FALSE
        )
    }, error = function(e) {
        # Catch error, usually due to low mutational burden
        print("Error: Failed to create Oncoplot!")
        print(e)
        plot.new()
        text(.5, .5, "Error: Failed to create Oncoplot!")
    }, finally = {
        # Close Plotting device
        dev.off()
    }
)
