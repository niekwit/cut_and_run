# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
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
# load PCA data 
data <- read.delim(snakemake@input[[1]], 
                   header = TRUE,
                   skip = 1)

#load sample information (replace any - with .)
sample_info <- read.csv("config/samples.csv", header = TRUE) %>%
  mutate( across(
    .cols = everything(),
    ~str_replace( ., "-", "." )
  ) )

#get unique genotypes, factors and treatments
genotypes <- unique(sample_info$genotype)
factors <- unique(sample_info$factor)
treatments <- unique(sample_info$treatment)

# set colours (genotypes)
if (length(genotypes) < 3) {
  colours <- c("#1B9E77","#D95F02")
} else {
  colours <- brewer.pal(length(genotypes), "Dark2")
}

# set shapes (treatments)
shapes <- c(21, 22, 24, 23, 25)[1:length(treatments)]

# set sizes (factors)
sizes <- c(8, 12, 16, 20, 24)[1:length(factors)]

# keep only components 1 and 2, transform and add sample information
df <- data[1:2,] %>%
  dplyr::select(-c("Component","Eigenvalue")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  rename(PC1 = 1,
         PC2 = 2) %>%
  left_join(sample_info, by = "sample")

# calculate variance explained for each PC
PC1_var <- round((data$Eigenvalue[1] / sum(data$Eigenvalue)) * 100, 1)
PC2_var <- round((data$Eigenvalue[2] / sum(data$Eigenvalue)) * 100, 1)

# create PCA plot
p <- ggplot(df, 
            mapping = aes(x = PC1, 
                          y = PC2,
                          fill = genotype, 
                          shape = treatment, 
                          size = factor)) +
  geom_point() +
  geom_label_repel(data = df , 
                   aes(label = sample,
                       fill = NULL), 
                   size = 5, 
                   nudge_x = 0.5, 
                   nudge_y = 0.5) +
  scale_fill_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  scale_size_manual(values = sizes) +
  theme_cowplot(16) +
  labs(x = paste0("PC1: ", PC1_var, "% variance"), 
       y = paste0("PC2: ", PC2_var, "% variance"), 
       Fill = "Genotype", 
       shape = "Treatment", 
       size = "Factor") 

# save plot
ggsave(snakemake@output[["pca"]], p)


#### Scree plot ####
# scale factor for utilising whole second y-axis range
# https://stackoverflow.com/questions/65559901/add-a-second-y-axis-to-ggplot
scalefactor <- max(data$Eigenvalue) / 100

# prepare data for scree plot
df <- data %>%
  dplyr::select(c("Component","Eigenvalue")) %>%
  mutate(Component = paste0("PC", Component)) %>%
  mutate(cumulative_variance = (cumsum(Eigenvalue) / sum(Eigenvalue) * 100 * scalefactor))

# create scree plot
s <- ggplot(df, aes(Component, cumulative_variance)) +
  geom_bar(aes(Component, Eigenvalue),
           stat="identity",
           colour = "black",
           fill = "aquamarine4") +
  geom_line(mapping = aes(x = Component, 
                          y = cumulative_variance,
                          group = 1),
            colour = "red",
            linewidth = 1) +
  geom_point(mapping = aes(x = Component, 
                          y = cumulative_variance),
            colour = "red",
            fill = "white",
            shape = 21,
            size = 6,
            stroke = 1.5) +
  theme_cowplot(16) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .x / scalefactor,
                                         breaks = seq(0, 100, 25),
                                         name = "Cumulative variance explained (%)"),
                     expand = expansion(mult = c(0, .05))) +
  labs(x = "Principal component", 
       y = "Eigenvalue")

# save plot
ggsave(snakemake@output[["scree"]], s)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

