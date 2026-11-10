# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# load required libraries
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(scales)

#### PCA plot ####
# Load PCA data
# check.names = FALSE keeps sample names as they are (no X prepended to names
# starting with a number), so they match samples.csv exactly
data <- read.delim(
  snakemake@input[["pca"]],
  header = TRUE,
  skip = 1,
  check.names = FALSE
)

# Load sample information
sample_info <- read.csv(
  snakemake@input[["samples"]],
  header = TRUE,
  colClasses = "character"
)

# Factor of each sample: IP samples get their own factor, control samples
# (only present in the control column) are shown in the facet of every
# factor they are a control for
factors <- sample_info %>%
  dplyr::select(sample, factor)
if ("control" %in% colnames(sample_info)) {
  controls <- sample_info %>%
    dplyr::filter(control != "") %>%
    dplyr::select(sample = control, factor) %>%
    dplyr::filter(!sample %in% factors$sample) %>%
    distinct()
  factors <- bind_rows(factors, controls)
}

# Unique sample conditions
samples <- colnames(data)[2:(ncol(data) - 1)]
samples <- unique(str_replace(samples, "_[0-9]+$", ""))

# Set colours for plotting
if (length(samples) == 1) {
  colours <- "#1B9E77"
} else if (length(samples) == 2) {
  colours <- c("#1B9E77", "#D95F02")
} else if (length(samples) < 9) {
  colours <- brewer.pal(length(samples), "Dark2")
} else if (length(samples) < 13) {
  colours <- brewer.pal(length(samples), "Set3")
} else {
  colours <- hue_pal()(length(samples))
}
names(colours) <- samples

# Keep only components 1 and 2, transform and add sample information
df <- data[1:2, ] %>%
  dplyr::select(-c("Component", "Eigenvalue")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  rename(PC1 = 1, PC2 = 2) %>%
  mutate(sample_condition = str_replace(sample, "_[0-9]+$", "")) %>%
  left_join(factors, by = "sample", relationship = "many-to-many") %>%
  mutate(factor = replace_na(factor, "Other"))

# Calculate variance explained for each PC
PC1_var <- round((data$Eigenvalue[1] / sum(data$Eigenvalue)) * 100, 1)
PC2_var <- round((data$Eigenvalue[2] / sum(data$Eigenvalue)) * 100, 1)

# Create PCA plot with one facet per factor (PCA is calculated on all samples,
# so all facets share the same axes)
p <- ggplot(df, mapping = aes(x = PC1, y = PC2, colour = sample_condition)) +
  geom_point(shape = 19, size = 5) +
  geom_label_repel(
    aes(label = sample),
    size = 3.5,
    nudge_x = 0.5,
    nudge_y = 0.5
  ) +
  scale_colour_manual(values = colours) +
  facet_wrap(~factor, ncol = min(3, n_distinct(df$factor))) +
  theme_cowplot(16) +
  labs(
    x = paste0("PC1: ", PC1_var, "% variance"),
    y = paste0("PC2: ", PC2_var, "% variance")
  ) +
  theme(legend.position = "none")

# Scale plot size with the number of facets
n_facets <- n_distinct(df$factor)
n_cols <- min(3, n_facets)
n_rows <- ceiling(n_facets / n_cols)

# Save plot
ggsave(
  snakemake@output[["pca"]],
  p,
  height = 4 * n_rows,
  width = 6 * n_cols,
  limitsize = FALSE
)

#### Scree plot ####
# Scale factor for utilising whole second y-axis range
# https://stackoverflow.com/questions/65559901/add-a-second-y-axis-to-ggplot
scalefactor <- max(data$Eigenvalue) / 100

# Prepare data for scree plot
df <- data %>%
  dplyr::select(c("Component", "Eigenvalue")) %>%
  mutate(Component = paste0("PC", Component)) %>%
  mutate(
    cumulative_variance = (cumsum(Eigenvalue) /
      sum(Eigenvalue) *
      100 *
      scalefactor)
  )

# Re-level PC factors to ensure correct order
df$Component <- factor(df$Component, levels = df$Component)

# Create scree plot
s <- ggplot(df, aes(Component, cumulative_variance)) +
  geom_bar(
    aes(Component, Eigenvalue),
    stat = "identity",
    colour = "black",
    fill = "aquamarine4"
  ) +
  geom_line(
    mapping = aes(x = Component, y = cumulative_variance, group = 1),
    colour = "red",
    linewidth = 1
  ) +
  geom_point(
    mapping = aes(x = Component, y = cumulative_variance),
    colour = "red",
    fill = "white",
    shape = 21,
    size = 5,
    stroke = 1.5
  ) +
  theme_cowplot(15) +
  theme(
    axis.title.y.right = element_text(color = "red"),
    axis.text.y.right = element_text(color = "red")
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(
      transform = ~ .x / scalefactor,
      breaks = seq(0, 100, 25),
      name = "Cumulative variance explained (%)"
    ),
    expand = expansion(mult = c(0, .05))
  ) +
  labs(x = "Principal component", y = "Eigenvalue")

# Save plot
ggsave(snakemake@output[["scree"]], s, height = 4, width = 6)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
