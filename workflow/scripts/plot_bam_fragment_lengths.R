# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
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
df <- data[rep(row.names(data), data$Occurrences), ]

# Add column for replicate number
df$replicate <- str_extract(df$Sample, "[0-9]+$")

# Plot density plot for fragment lengths with panel for each condition
layouts <- list(
  `2` = c(2, 1),
  `3` = c(3, 1),
  `4` = c(2, 2),
  `5` = c(3, 2),
  `6` = c(3, 2),
  `7` = c(3, 3),
  `8` = c(4, 2),
  `9` = c(3, 3),
  `10` = c(3, 4),
  `11` = c(3, 4),
  `12` = c(3, 4)
)

n <- as.character(length(conditions))

if (!n %in% names(layouts)) {
  stop("Too many conditions to plot")
}

ncol <- layouts[[n]][1]
nrow <- layouts[[n]][2]

p <- ggplot(df, aes(x = Size, fill = replicate)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~condition, ncol = ncol, nrow = nrow) +
  theme_cowplot(12) +
  theme(
    strip.background = element_rect(
      color = "black",
      fill = "grey90",
      size = 1
    )
  ) +
  scale_fill_manual(
    values = brewer.pal(length(unique(df$replicate)), "Dark2")
  ) +
  labs(
    title = "BAM fragment lengths",
    x = "Fragment Length",
    y = "Density"
  )

# Save plot
ggsave(
  snakemake@output[["pdf"]],
  p,
  height = length(conditions) * 3,
  width = 6
)
