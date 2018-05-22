library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
datax <- fread(args[1], sep = "\t", header = F)
colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues")

# ymax <- quantile(datax$NormalizedValues, probs = 0.95)
myplot <- ggplot(data = datax, aes(x = genes, y = NormalizedValues, fill = genes)) + 
	geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
	# scale_y_continuous(limits = c(0, ymax)) +
	labs(title = args[1], caption = date()) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "boxplot.png", sep = "."), width = 400, height = 600)
print(myplot)
dev.off()
