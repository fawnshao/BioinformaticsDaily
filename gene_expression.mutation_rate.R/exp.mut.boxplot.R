library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
datax <- fread(args[1], sep = "\t", header = F)
colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues", "types", 
	"MutRegions", "MutCounts", "ExpressionQuantile")

targets <- c("exon", "intron", "Intergenic", "promoter-TSS")
for(i in 1:length(targets)){
	datay <- datax[MutRegions == targets[i],]
	if(nrow(datay) > 10){
		ymin <- min(datay$MutCounts)
		ymax <- quantile(datay$MutCounts, probs = 0.95)
		myplot <- ggplot(data = datay, aes(x = genes, y = MutCounts, fill = ExpressionQuantile)) + 
			geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) +
			scale_y_continuous(limits = c(ymin, ymax)) +
			labs(title = paste(args[1], targets[i], sep = " "), caption = date()) + 
			theme(axis.text.x = element_text(angle = 60, hjust = 1))
		png(filename = paste(args[1], targets[i], "boxplot.png", sep = "."), width = 400, height = 600)
		print(myplot)
		dev.off()
	}
}
