# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(DiffBind)

# Load Snakemake variables
bams <- snakemake@input["bam"]
xls_peaks <- snakemake@input["xls"]
sample_sheet <- snakemake@output["sh"]
dba_file <- snakemake@output["dba"]
pca <- snakemake@output["pca"]
sc <- snakemake@output["sc"]
pp <- snakemake@output["pp"]

# Load sample information
samples <- read.csv("config/samples.csv")

# Create combined factor/condition/treatment column to create simple model
if (
  length(unique(samples$genotype)) > 1 & length(unique(samples$treatment)) > 1
) {
  samples$comb <- paste(samples$genotype, samples$treatment, sep = "_")
} else if (
  length(unique(samples$genotype)) > 1 & length(unique(samples$treatment) == 1)
) {
  samples$comb <- samples$genotype
} else if (
  length(unique(samples$genotype)) == 1 & length(unique(samples$treatment) > 1)
) {
  samples$comb <- samples$treatment
}

# Create base DiffBind sample table
db.table <- data.frame(
  SampleID = samples$sample,
  Factor = samples$factor,
  Genotype = samples$genotype,
  Treatment = samples$treatment,
  Condition = samples$comb
)
db.table$Replicate <- str_split_fixed(samples$sample, "_", 2)[, 2]

# Get sample bam files
for (i in samples$sample) {
  db.table[db.table$SampleID == i, "bamReads"] <- bams[grep(i, bams)]
}

# Get control samples and data
if ("control" %in% names(samples)) {
  for (i in samples$control) {
    db.table$ControlID <- samples$control
    db.table[db.table$ControlID == i, "bamControl"] <- bams[grep(i, bams)]
  }
}

# Add peak information
for (i in samples$sample) {
  db.table[db.table$SampleID == i, "Peaks"] <- xls_peaks[grep(i, xls_peaks)]
  db.table$"PeakCaller" <- "macs"
}

# Remove redundant columns
if (all(db.table$Condition == samples$genotype)) {
  db.table$Genotype <- NULL
} else if (all(db.table$Condition == samples$treatment)) {
  db.table$Treatment <- NULL
}

# Get reference conditions
references <- unique(samples[samples$reference == "yes", "comb"])

# Save sample sheet to file
write.table(
  db.table,
  file = sample_sheet,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Perform DiffBind analysis
dba <- dba(sampleSheet = db.table) %>%
  dba.count() %>%
  dba.normalize() %>%
  dba.contrast(minMembers = 2) %>%
  dba.analyze(design = ~Condition, bGreylist = FALSE)

# Save diffbind object to file
save(dba, file = dba_file)

# Generate correlation heatmap
pdf(sc)
plot(dba)
dev.off()

# Plot PCA
pdf(pca)
dba.plotPCA(dba, DBA_CONDITION, label = DBA_ID)
dev.off()

# Perform differential binding analysis
for (i in dba$contrasts) {
  contrast <- i$contrast
  print(contrast)
  #dba.temp <- dba.analyze(dba, contrast = i)
}


# Profile plot
pdf(pp)
dba.plotProfile(dba, samples = dba$masks$All, doPlot = TRUE)
dev.off()


# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
