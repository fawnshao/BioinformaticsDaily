library(data.table)
library(ggplot2)
args <- commandArgs(TRUE) # c("WGS.mut.KAT5.exp_seq.READ-US.tsv")
datax <- fread(args[1], sep = "\t", header = F)
colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues", "types", "MutRegions", "MutCounts")

targets <- c("exon", "Intergenic", "promoter-TSS")
for(i in 1:length(targets)){
	datay <- datax[MutRegions == targets[i],]
	if(nrow(datay) > 10){
		xmax <- quantile(datay$MutCounts, probs = 0.95)
		myplot <- ggplot(data = datay, aes(x = MutCounts, y = NormalizedValues, colour = types, shape = genes)) + 
			geom_point() + geom_smooth() +
			scale_x_continuous(limits = c(0, xmax)) +
			labs(title = paste(args[1], targets[i], sep = " "), caption = date()) + 
			theme(axis.text.x = element_text(angle = 60, hjust = 1))
		png(filename = paste(args[1], targets[i], "scatterplot.png", sep = "."), width = 1000, height = 600)
		print(myplot)
		dev.off()
	}
}
