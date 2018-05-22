library(data.table)
library(ggplot2)
args <- commandArgs(TRUE) # c("WGS.mut.KAT5.exp_seq.READ-US.tsv")
datax <- fread(args[1], sep = "\t", header = F)
colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues", "types", "MutRegions", "MutCounts")

targets <- c("exon", "Intergenic", "promoter-TSS")
for(i in 1:length(targets)){
	myplot <- ggplot(data = datax[MutRegions == targets[i],], aes(x = MutCounts, y = NormalizedValues, colour = types, shape = genes)) + 
		geom_point() + 
		labs(title = paste(args[1], targets[i], sep = " "), caption = date()) + 
		theme(axis.text.x = element_text(angle = 60, hjust = 1))
	png(filename = paste(args[1], targets[i], "scatterplot.png", sep = "."), width = 600, height = 600)
	print(myplot)
	dev.off()
}