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
		ymax <- quantile(datay$MutCounts, probs = 0.9)
		myplot <- ggplot(data = datay, aes(x = genes, y = MutCounts, fill = ExpressionQuantile)) + 
			geom_violin(position = position_dodge(1)) +
			geom_dotplot(binaxis = 'y', bins = 5, stackdir = 'center', position = position_dodge(1)) + 
			scale_y_continuous(limits = c(ymin, ymax)) +
			labs(title = paste(args[1], targets[i], sep = " "), caption = date()) + 
			theme(axis.text.x = element_text(angle = 60, hjust = 1))
		png(filename = paste(args[1], targets[i], "violinplot.png", sep = "."), width = 600, height = 600)
		print(myplot)
		dev.off()
	}
}
