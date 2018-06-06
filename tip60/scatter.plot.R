library(data.table)
library(ggplot2)
args <- commandArgs(TRUE) # c("WGS.mut.KAT5.exp_seq.READ-US.tsv")
datax <- fread(args[1], sep = "\t", header = F)
colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues", "types", "MutCounts")


if(nrow(datax) > 20){
	xmin <- min(datax$MutCounts)
	xmax <- quantile(datax$MutCounts, probs = 0.95)
	myplot <- ggplot(data = datax, aes(x = MutCounts, y = NormalizedValues, colour = types, shape = genes)) + 
		geom_point() + geom_smooth(method = lm) +
		scale_x_continuous(limits = c(xmin, xmax)) +
		labs(title = paste(args[1], sep = " "), caption = date()) + 
		theme(axis.text.x = element_text(angle = 60, hjust = 1))
	png(filename = paste(args[1], "scatterplot.png", sep = "."), width = 1000, height = 600)
	print(myplot)
	dev.off()
}