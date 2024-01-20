# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(tidyverse)
library(cowplot)
library(ggrepel)

# Load Snakemake variables
dir <- snakemake@output[["dir"]]
fc <- snakemake@params[["fc"]]
fdr <- -log10(snakemake@params[["fdr"]])

# Load csv file with DESeq2 results
csv.files <- Sys.glob(file.path(dir, "*.csv"))

# Create df for text labels in plot
df.label <- df %>%
  filter(log2FoldChange >= fc | log2FoldChange <= -fc,
         -log10(padj) > fdr) %>%
  filter(!is.na(gene_name)) %>%
  arrange(log2FoldChange) 

# Plot Volcano for each csv file
for (csv in csv.files) {
    # Load data
    df <- read.csv(csv)
  
    # Plot
    p <- ggplot(df, aes(x = log2FoldChange, 
                        y = -log10(padj),
                        colour = case_when(
                          log2FoldChange >= fc & -log10(padj) >= fdr ~ "red",
                          log2FoldChange <= -fc & -log10(padj) >= fdr ~ "navy",
                          .default = "black"))) +
        scale_colour_manual(values=c("black", "blue", "red"),
                          guide = "none") +
        geom_point(alpha = 0.5,
                   size=6) +
        geom_hline(yintercept = fdr, 
                   linetype = "dashed") +
        geom_vline(xintercept = c(-fc, fc), 
                   linetype = "dashed") +
        labs(x = "log2 Fold Change", 
             y = "-log10 adjusted p-value") +
        theme_cowplot(18) +
        theme(legend.position = "none")
    
}



# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")