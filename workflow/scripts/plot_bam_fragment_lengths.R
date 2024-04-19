# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(RColorBrewer)

# Load data
data <- read.delim(snakemake@input[[1]], header = TRUE, skip = 1)

# Add condition column to data
conditions <- unique(str_replace(data$Sample, "_[0-9]+$", ""))
data$condition <- str_replace(data$Sample, "_[0-9]+$", "")

# Based on value in occurrences column, copy that number of rows to new df
# (needed for plotting)
df <- data[rep(row.names(data), data$Occurrences),]

# Add column for replicate number
df$replicate <- str_extract(df$Sample, "[0-9]+$")

# Plot density plot for fragment lengths with panel for each condition
if (length(conditions) == 2) {
       ncol <- 2
       nrow <- 1
} else if (length(conditions) == 3) {
       ncol <- 3
       nrow <- 1
} else if (length(conditions) == 4) {
       ncol <- 2
       nrow <- 2
} else if (length(conditions) == 5) {
       ncol <- 3
       nrow <- 2
} else if (length(conditions) == 6) {
       ncol <- 3
       nrow <- 2
} else if (length(conditions) == 7) {
       ncol <- 3
       nrow <- 3
} else if (length(conditions) == 8) {
       ncol <- 4
       nrow <- 2
} else if (length(conditions) == 9) {
       ncol <- 3
       nrow <- 3
} else if (length(conditions) == 10) {
       ncol <- 3
       nrow <- 4
} else if (length(conditions) == 11) {
       ncol <- 3
       nrow <- 4
} else if (length(conditions) == 12) {
       ncol <- 3
       nrow <- 4
} else {
       stop("Too many conditions to plot")
}

p <- ggplot(df, aes(x = Size,
                    fill = replicate)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~condition,
             ncol = ncol,
             nrow = nrow) +
  theme_cowplot(12) +
  theme(strip.background = element_rect(color = "black",
                                        fill = "grey90",
                                        size = 1)) +
  scale_fill_manual(values = brewer.pal(length(unique(df$replicate)), "Dark2")) +
  labs(title = "BAM fragment lengths",
       x = "Fragment Length",
       y = "Density")

# Save plot
ggsave(snakemake@output[["pdf"]],
       p,
       height = length(conditions) * 0.6,
       width = 6)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")