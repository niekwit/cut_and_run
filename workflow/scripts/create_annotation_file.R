# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(rtracklayer)
library(tidyverse)
library(GenomicFeatures)

# Load Snakemake parameters
gtf <- snakemake@input[["gtf"]]

# Create EDB
# ------------------------------

# Load GTF file
db <- rtracklayer::import(gtf)

# Extract relevant information
edb <- data.frame(
  geneId = db$gene_id,
  geneName = db$gene_name,
  geneBiotype = db$gene_biotype
) %>%
  distinct()

# Save df to file as R object
save(edb, file = snakemake@output[["rdata"]])

# Create Txdb
# ------------------------------
txdb <- makeTxDbFromGFF(gtf)

# Save Txdb to file as R object
save(txdb, file = snakemake@output[["txdb"]])
