# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)

# Load Snakemake parameters
xls <- snakemake@input[["xls"]]
DB <- snakemake@input[["adb"]]
peak.mode <- snakemake@params[["pm"]]
txt <- snakemake@output[["txt"]]
bed.file <- snakemake@output[["bed"]]
genome <- snakemake@params[["genome"]]

# Determine how many lines to skip in xls file (comment lines)
if (peak.mode == "narrow") {
  skip <- 21
} else if (peak.mode == "broad") {
  skip <- 22
} 

# Load xls file and convert to bed
xls <- read.table(xls, 
                  skip = skip,
                  header = TRUE)
bed <- xls %>%
  select(c("chr", "start", "end", "name", "fold_enrichment")) %>%
  mutate(strand = ".")

# Write bed to file (annotatePeak requires a file)
write.table(bed,
            file = bed.file,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

# Load annotation database
txdb <- makeTxDbFromGFF(gtf)

# Annotate bed file
peakAnno <- annotatePeak(bed.file,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb
                        )

# Plot binding relative to TSS
pdf(snakemake@output[["tss"]])
plotDistToTSS(peakAnno,
              title="Distribution of peaks loci\nrelative to TSS") +
  theme(axis.line.y = element_line(size = 0),
        axis.line.x = element_line(size = 0.5),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20))
dev.off()

# Plot annotation bar
pdf(snakemake@output[["bar"]])
plotAnnoBar(peakAnno) +
  theme(axis.line.y = element_line(size = 0),
        axis.line.x = element_line(size = 0.5),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20))
dev.off()

# Tidy up annotation data
df <- as.data.frame(peakAnno@anno@elementMetadata@listData)
names(df)[1:3] <- c("peak_id","fold_enrichment","strand")

# Add gene names and gene biotype to annotation
load(DB)
df <- df %>%
  left_join(edb, by = "geneId")

# Write annotation to file
write.table(df,
            file = txt,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")