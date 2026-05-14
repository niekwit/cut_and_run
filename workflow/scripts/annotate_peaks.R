# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)

# Load Snakemake parameters
bed_file <- snakemake@input[["bed"]]
edb_file <- snakemake@input[["edb"]]
peak.mode <- snakemake@params[["pm"]]
txt <- snakemake@output[["txt"]]
gtf <- snakemake@params[["gtf"]]

# Load annotation database
txdb <- makeTxDbFromGFF(gtf)

# Annotate bed file
peakAnno <- annotatePeak(bed_file, tssRegion = c(-3000, 3000), TxDb = txdb)

# Tidy up annotation data
df <- as.data.frame(peakAnno@anno@elementMetadata@listData)
names(df)[1:3] <- c("peak_id", "fold_enrichment", "strand")

# Add gene names and gene biotype to annotation
load(edb_file)
df <- df %>%
  left_join(edb, by = "geneId")

# Write annotation to file
write.table(df, file = txt, quote = FALSE, sep = "\t", row.names = FALSE)
