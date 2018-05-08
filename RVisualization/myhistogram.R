args <- commandArgs(TRUE)
library(ggplot2)
library(reshape2)
library(data.table)

inputfile <- args[1]
coltoplot <- as.numeric(args[2])
binsize <- as.numeric(args[3])

# source("/home1/04935/shaojf/myScripts/ggplot2_multiplot.R")
# binsize <- 100

inputdata <- fread(inputfile, sep = "\t")
plotdata <- inputdata[, coltoplot]
colnames(plotdata) <- "myvalues"
p0 <- ggplot(data = plotdata, aes(x = myvalues)) + 
  geom_histogram(bins = binsize) +
  # scale_x_continuous(limits = c(-1e5, 1e5)) +
  xlab("Values") + ggtitle(args[1])

pdf(file = paste(inputfile, "hist.pdf", sep = "."), width = 10, height = 6)
print(p0)
dev.off()