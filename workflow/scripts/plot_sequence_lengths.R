# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

#args <- commandArgs(trailingOnly=TRUE)
#setwd(args[1])
#df <- read.csv(file=args[2])
#base.file <- args[3]
#outdir <- args[4]

x_low <- min(df$read_length)
x_high <- max(df$read_length)
x_range <- x_low:x_high
mean.length <- mean(df$read_length)
list <- suppressWarnings(hist(df$read_length, 
             breaks=x_range, 
             freq=TRUE,
             plot=FALSE))
y.value.label <- max(list$counts)/sum(list$count)

p <- ggplot(data=df, aes(x=read_length)) +
  geom_bar(aes(y=stat(count)/sum(count)), 
           colour="black", 
           fill="grey50") + 
  scale_x_continuous(breaks=x_range)+ 
  theme_cowplot(16) +
  theme(plot.title = element_text(size = 16)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Read length") +
  ylab("Relative frequency") +
  ggtitle(paste0("Read length distribution ",base.file)) +
  annotate("text", x=x_range[3], 
           y=y.value.label, 
           label= paste0("Mean read length: ",
                         round(mean.length,1)), 
           size=5)

ggsave(filename=paste0(outdir,"/","Read-length-frequency",base.file,".png"),
       plot=p )

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")