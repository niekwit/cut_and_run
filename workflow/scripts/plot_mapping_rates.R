# redirect R output to log
slog <- file(snakemake@log[[1]], open="wt")
sink(slog, type = "output")
sink(slog, type = "message")

# Load libraries
library(tidyverse)
library(cowplot)
library(reshape2)

# Get log files
log.files <- unlist(snakemake@input["log"])
#print(log.files)
# Create empty data frames to store values
df.perc <- data.frame(sample = character(),
                      overall_alignment_rate = numeric(),
                      aligned_multiple_perc = numeric(),
                      aligned_one_perc = numeric(),
                      aligned_zero_perc = numeric())

df.number <- data.frame(sample = character(),
                        total_reads = numeric(),
                        aligned_multiple_number = numeric(),
                        aligned_one_number = numeric(),
                        aligned_zero_number = numeric())

# Read in mapping rates from log files
for (i in seq_along(log.files)) {
  # Read in log file
  log <- readLines(log.files[[i]])
  
  # Get sample name
  sample <- str_replace(basename(log.files[[i]]), ".log", "")

  # Get total number of reads
  total_reads <- as.numeric(str_replace(log[grep(" reads; of these:", log)],
                                        " reads; of these:",
                                        ""))
  
  # Get overall alignment rate
  overall_alignment_rate <- as.numeric(str_replace(log[grep("% overall alignment rate", log)], 
                                                   "% overall alignment rate",
                                                   ""))

  # Get number and percentage of reads aligned to 0, 1, or multiple locations
  read_log <- function(pattern){
    suppressWarnings({
      aligned_number <- as.numeric(str_split_1(str_trim(log[grep(pattern, log)],
                                                           side = "left"),
                                                  " "))[1]
      aligned_perc <- str_split_1(str_trim(log[grep(pattern, log)],
                                              side = "left"),
                                     " ")[2]
      aligned_perc <- as.numeric(sub("(", "", sub("%)", "", aligned_perc), fixed = TRUE))
    
      return(c(aligned_number, aligned_perc))
    })
  }
  aligned_zero_number <- read_log("aligned concordantly 0 times")[1]
  aligned_zero_perc <- read_log("aligned concordantly 0 times")[2]
  
  aligned_one_number <- read_log("aligned concordantly exactly 1 time")[1]
  aligned_one_perc <- read_log("aligned concordantly exactly 1 time")[2]
  
  aligned_multiple_number <- read_log("aligned concordantly >1 times")[1]
  aligned_multiple_perc <- read_log("aligned concordantly >1 times")[2]
  
  # Add data to data frames
  df.perc <- df.perc %>%
    add_row(sample = sample,
            overall_alignment_rate = overall_alignment_rate,
            aligned_multiple_perc = aligned_multiple_perc,
            aligned_one_perc = aligned_one_perc,
            aligned_zero_perc = aligned_zero_perc)
  
  df.number <- df.number %>%
    add_row(sample = sample,
            total_reads = total_reads,
            aligned_multiple_number = aligned_multiple_number,
            aligned_one_number = aligned_one_number,
            aligned_zero_number = aligned_zero_number)
}

# Melt data frames
df.perc <- melt(df.perc, id.vars = "sample")
df.number <- melt(df.number, id.vars = "sample")

title <- "Concordant alignment rates (Bowtie2)"

# Plot and save alignment data
plot_data <- function(df, title, y.axis, outfile){
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  p <- ggplot(df, aes(x = sample, y = value, fill = variable)) +
    geom_bar(stat = "identity",
             position = "dodge",
             colour = "black") +
    theme_cowplot(18) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 12),
          plot.margin = margin(t = 0.5, r = 1.5, b = 0.5, l = 0.5, unit = "cm")) +
    labs(title = title,
         x = NULL,
         y = y.axis, #"Number of reads aligned",
         fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c("darkblue",
                                 "grey40",
                                 "forestgreen",
                                 "red"),
                      labels = c("Total aligned",
                                 "Aligned > 1 time",
                                 "Aligned 1 time",
                                 "Did not align"))
  
  # Save to file
  ggsave(outfile, p)
}

plot_data(df.perc, "Concordant alignment rates (Bowtie2)", "Percentage of total reads", snakemake@output["rates"][[1]])
plot_data(df.number, "Concordant aligned read counts (Bowtie2)", "Number of reads", snakemake@output["counts"][[1]])

# Close redirection of output/messages
sink(slog, type = "output")
sink(slog, type = "message")
