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
deseq2_apply_control <- snakemake@params[["control"]]
count.files <- snakemake@input[["counts"]]
print(count.files)

# Directory with htseq-count data files
count.dir <- dirname(count.files[1])

# Load sample information
samples <- read.csv("config/samples.csv")
names <- samples$sample
genotypes <- unique(samples$genotype)
treatments <- unique(samples$treatment)
controls <- unique(samples[(tolower(samples$factor) %in% c("input", "igg")), ]$sample)

if (length(treatments) > 1) {
  samples$comb <- paste0(samples$genotype, "_", samples$treatment)
} else {
  samples$comb <- paste0(samples$genotype)
}

# Get reference condition(s)
references <- unique(samples[samples$reference == "yes" ,]$comb)
if (length(references) == 0){
  stop("ERROR: No reference samples found. Please check your samples.csv file.")
}

# Check if batch column exists
if ("batch" %in% colnames(samples)) {
  batches <- unique(samples$batch)
  if (length(batches) == 1) {
    batches <- 1
    samples$batch <- 1
  }
} else {
  batches <- 1
  samples$batch <- 1
}

# Create df for DESeq2
sampleFiles <- vector()
for (i in seq_along(names)) {
  sampleFiles[i] <- count.files[grep(names[i], count.files)]
}

sampleTable <- data.frame(sampleName = names,
                          fileName = sampleFiles,
                          condition = samples$comb,
                          batch = samples$batch,
                          reference = samples$reference)

# Create DESeqDataSet
if (batches == 1) {
  ddsHTSeqCount <- DESeqDataSetFromHTSeqCount(sampleTable,
                                              design = ~condition)
} else {
  ddsHTSeqCount <- DESeqDataSetFromHTSeqCount(sampleTable,
                                              design = ~batch + condition)
}

# save DESeqDataSet to file before filtering
save(ddsHTSeqCount, file = snakemake@output[["rdata"]])

# Check if input or IgG controls are present and match to IP samples
if (deseq2_apply_control == "True") {
  if (is.null(controls)) {
    print("WARNING: no input/IgG control samples found!")
    print("Differential peaks will be calculated without control...")
    deseq2_apply_control <- "False"
  } else {
    #TO DO
    #CORRECT DATA FOR CONTROL
  }
}

# Load one htseq-count data set for annotation (gene_name and biotype are not kept by DESeq2)
annotation <- read.csv(count.files[1], sep = "\t", header = FALSE) %>%
  select(V1, V2, V3)
names(annotation) <- c("ensembl_gene_id", "gene_name", "biotype")

# Filter out genes with very low read counts
smallestGroupSize <- snakemake@params[["sg"]]
cumulativeFilterOut <- snakemake@params[["cfo"]]
keep <- rowSums(counts(ddsHTSeqCount) >= cumulativeFilterOut) >= smallestGroupSize
ddsHTSeqCount <- ddsHTSeqCount[keep, ]

# List to store data from each contrast (will be nested)
resList <- list()

for (r in seq_along(references)) {
  # Set reference condition
  ddsHTSeqCount$condition <- relevel(ddsHTSeqCount$condition, ref = references[r])

  # Run DESeq2
  dds <- DESeq(ddsHTSeqCount)
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
                      name=contrast)
    
    # Print summary
    summary(resIHW)
    
    # Create df for results and add annotation
    df <- as.data.frame(res) %>%
      mutate(ensembl_gene_id = res@rownames, .before=1) %>%
      left_join(annotation, by = "ensembl_gene_id") %>%
      relocate(ensembl_gene_id, gene_name, biotype)
    
    resList_contrast[[str_replace(contrast, "condition_", "")]] <- df
  }
  resList[[references[r]]] <- resList_contrast
}

# function to flatten nested lists (https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists/41882883#41882883)
flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1] == "list")
  out <- c(x[!morelists], unlist(x[morelists], recursive = FALSE))
  if(sum(morelists)) {
    Recall(out)
  } else {
    return(out)
  }
}

# Flatten resList
resList <- flattenlist(resList)

# Remove reference name followed by dot from list element names
names(resList) <- str_replace(names(resList), paste0("^", references, "\\."), "")

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

