args <- commandArgs(TRUE)
library(ggplot2)
library(reshape2)
library(data.table)

inputfile <- args[1]
coltogroup <- as.numeric(args[2])
coltoplot <- as.numeric(args[3])
minx <- as.numeric(args[4])
maxx <- as.numeric(args[5])

inputdata <- fread(inputfile, sep = "\t")
plotdata <- inputdata[, c(coltogroup, coltoplot)]
colnames(plotdata) <- c("mygroup", "myvalues")
p0 <- ggplot(plotdata, aes(x = mygroup, y = myvalues, fill = factor(mygroup))) + 
    geom_violin() + ylim(c(minx, maxx)) +
    xlab("") + ylab("myvalues") +
    guides(fill = FALSE) +
    ggtitle(args[1])
pdf(file = paste(inputfile, "violin.pdf", sep = "."), width = 10, height = 6)
print(p0)
dev.off()
