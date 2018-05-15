args <- commandArgs(TRUE)
# args <- c("NUP_ATAC.peak.length.tsv", "1", "2")
library(ggplot2)
library(data.table)

inputfile <- args[1]
coltogroup <- as.numeric(args[2])
coltoplot <- as.numeric(args[3])

inputdata <- fread(inputfile, sep = "\t")
plotdata <- inputdata[, .SD, .SDcols = c(coltogroup, coltoplot)]
colnames(plotdata) <- c("mygroup", "myvalues")
theme_set(theme_classic())
p0 <- ggplot(plotdata, aes(myvalues)) + 
    geom_density(aes(fill = mygroup), alpha = 0.3) +
    # geom_density(aes(colour = mygroup), alpha = 0.8) +
    xlim(0, 1000) + 
    theme(legend.position = "top", legend.box = "horizontal") + 
    labs(title = "Density plot", 
         subtitle = args[1],
         caption = date(),
         x = "Peak Length",
         fill = "# Experiment")
pdf(file = paste(inputfile, "density.pdf", sep = "."), width = 15, height = 6)
print(p0)
dev.off()
