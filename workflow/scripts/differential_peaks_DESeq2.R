# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(DESeq2)
library(IHW)
library(stringr)
library(dplyr)
library(openxlsx)

# Load Snakemake variables
alpha <- snakemake@params[["alpha"]]
all.count.files <- snakemake@input[["counts"]]
lfc <- snakemake@params[["fc"]]

# Load sample information
samples <- read.csv("config/samples.csv")
names <- samples$sample
control_names <- samples$control
genotypes <- unique(samples$genotype)
treatments <- unique(samples$treatment)

if (length(treatments) > 1) {
  samples$condition <- paste0(samples$genotype, "_", samples$treatment)
} else {
  samples$condition <- samples$genotype
}

# Create sampleTable for DESeq2
references <- unique(control_names)
sampleName <- c(names, references)
reference.condition <- gsub("_[^_]+$", "", references)
condition <- c(samples$condition, reference.condition)
fileName <- vector()
for (i in seq_along(sampleName)) {
  fileName[i] <- all.count.files[grep(sampleName[i], all.count.files)]
}
sampleTable <- data.frame(sampleName = sampleName,
                        fileName = fileName,
                        condition = condition)

# Build the DESeqDataSet
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = ".",
                                       design = ~ condition)

# Save DESeqDataSet to file before filtering
save(ddsHTSeq, file = snakemake@output[["rdata"]])

# Load one htseq-count data set for annotation (gene_name and biotype are not kept by DESeq2)
annotation <- read.csv(all.count.files[1], 
                       sep = "\t", 
                       header = FALSE) %>%
  select(V1, V2, V3)
names(annotation) <- c("ensembl_gene_id", "gene_name", "biotype")

# Filter out genes with very low read counts
smallestGroupSize <- snakemake@params[["sg"]]
cumulativeFilterOut <- snakemake@params[["cfo"]]
keep <- rowSums(counts(ddsHTSeq) >= cumulativeFilterOut) >= smallestGroupSize
ddsHTSeq <- ddsHTSeq[keep, ]

# List to store data from each contrast
resList <- list()

for (r in seq_along(unique(reference.condition))) {
  # Set reference condition
  print(paste0("Setting reference condition to ", unique(reference.condition)[r], "..."))
  ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = unique(reference.condition)[r])

  # Run DESeq2
  dds <- DESeq(ddsHTSeq)
  res <- results(dds)
  
  # Get all contrasts
  contrasts <- resultsNames(dds)
  contrasts <- strsplit(contrasts," ")
  contrasts[1] <- NULL
  
  resList_contrast <- list()
  for (i in seq_along(contrasts)) {
    
    # Independent hypothesis weighting
    contrast <- contrasts[[i]]
    print(paste0("Analysing contrast ", contrast, "..."))
    resIHW <- results(dds, 
                      filterFun = ihw, 
                      alpha = alpha,
                      name = contrast,
                      lfcThreshold = lfc)
    # Print summary
    summary(resIHW)
    
    # Create df for results and add annotation
    df <- as.data.frame(resIHW) %>%
      mutate(ensembl_gene_id = res@rownames, .before = 1) %>%
      left_join(annotation, by = "ensembl_gene_id") %>%
      mutate(contrast = gsub("condition_", "", contrast)) %>%
      relocate(contrast, ensembl_gene_id, gene_name, biotype)
    
    resList[[(length(resList) + 1)]] <- df
  }
}

# Name each df in resList
names(resList) <- lapply(resList, function(df){unique(df$contrast)})

# Save each df in resList to csv
path <- dirname(snakemake@output[["xlsx"]])
dir.create(path, showWarnings = FALSE)
for (i in seq_along(resList)){
  write.csv(resList[[i]],
            paste0(path, "/", names(resList)[i], ".csv"),
            row.names = FALSE)
}

# Also write output to one xlsx file
write.xlsx(resList,
           snakemake@output[["xlsx"]],
           colNames = TRUE)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")