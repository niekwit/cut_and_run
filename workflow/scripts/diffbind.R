# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(DiffBind)

# Load Snakemake variables
bams <- snakemake@input[["ip_bam"]]
xls_peaks <- snakemake@input[["xls"]]
input_available <- snakemake@params[["input"]]
IgG_available <- snakemake@params[["igg"]]



# Load sample information
samples <- read.csv("config/samples.csv")

# Create empty data frame for DiffBind
sampleTable <- data.frame(SampleID = NA,
                          Tissue = NA,
                          Factor = NA,
                          Condition = NA,
                          Treatment = NA,
                          Replicate = NA,
                          bamReads = NA,
                          Peaks = NA,
                          PeakCaller = NA)
# add ControlID and bamControl conditionally



# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
