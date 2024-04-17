# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(DiffBind)

# Load Snakemake variables
bams <- snakemake@input["bam"]
xls_peaks <- snakemake@input["xls"]
control_available <- snakemake@params["control"]

# Load sample information
samples <- read.csv("config/samples.csv")

# Create combined factor/condition/treatment column to create simple model
samples$comb <- paste(samples$factor, samples$condition, samples$treatment, sep = "_")

# Create DiffBind sample table
db.table <- data.frame() %>%
  mutate(SampleID = samples$sample,
         Tissue = NA,
         Factor = samples$factor,
         Condition = samples$comb,
         Treatment = samples$treatment,
         #Replicate = samples$Replicate,
         bamReads = bams,
         Peaks = xls_peaks,
         PeakCaller = "macs")

# Add ControlID and bamControl conditionally
# and load files from Snakemake
if (control_available == "True") {
  db.table <- db.table %>%
    mutate(ControlID = snakemake@wildcards["control_sample"],
           bamControl = snakemake@input["control_bam"])

  control.bams <- snakemake@input["control_bam"]
}

# Add replicate number
db.table <- db.table %>%
  group_by(Condition) %>%
  mutate(Replicate = 1:row_number(.))

print(db.table)

# Perform DiffBind analysis
dba <- tamoxifen <- dba(sampleSheet="tamoxifen.csv") %>%
  dba.count() %>%
  dba.normalize() %>%
  dba.contrast(minMembers = 2) %>%
  dba.analyze(design = "~Condition", #simple model for now, implement more complex model later
              bGreylist = FALSE)

# Save diffbind object to file 
save(dba, file = snakemake@output["dba"])

# Plot PCA
pdf(snakemake@output["pca"])
dba.plotPCA(dba, DBA_CONDITION, label = DBA_ID)
dev.off()

# Plot correlation heatmap
pdf(snakemake@output["sc"])
plot(dba)
dev.off()

# Profile plot
pdf(snakemake@output["pp"])
dba.plotProfile(dba, samples = dba$masks$All, doPlot = TRUE)
dev.off()

# Get all reference conditions
ref_conditions <- samples[samples$reference == "yes", "comb"]
if (length(ref_conditions) == 0) {
  stop("No reference samples found in samples.csv")
}

# Perform differential binding analysis for each
# pairwise comparison to each reference conditions
for (reference in ref_conditions) {
  # Set control condition
  diffbind <- dba.contrast(dba, reorderMeta = list(Condition = reference))
  
  # Get number of contrasts
  contrasts <- length(diffbind$contrasts)
}

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
