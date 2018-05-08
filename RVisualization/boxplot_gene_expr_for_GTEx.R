args <- commandArgs(TRUE)
library(ggplot2)
# args <- c("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", "33319", "GTEx_sample.tissue.txt")
datafile <- args[1]
linenumber <- as.numeric(args[2])
tissuesfile <- args[3]
x <- scan(file = datafile, skip = linenumber - 1, nlines  = 1, sep = "\t", what = matrix(""))
info <- as.matrix(read.table(tissuesfile, sep = "\t"))

data <- data.frame(info[,2], log2(as.numeric(x[-1])+1))
genename <- x[1]
# data <- melt(expr)
colnames(data) <- c("tissue", "expression")
png(filename = paste(genename, "violin.png", sep = "."), width = 1500, height = 600)
ggplot(data, aes(x = tissue, y = expression, fill = tissue)) + 
	geom_violin() + 
	stat_summary(fun.data = "mean_sdl", geom = "pointrange", color = "black") +
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()
