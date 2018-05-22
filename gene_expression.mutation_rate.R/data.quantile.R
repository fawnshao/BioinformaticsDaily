library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
datax <- fread(args[1], sep = "\t", header = F)
colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues", "types")

ymax <- quantile(datax$NormalizedValues, probs = 0.95)
myplot <- ggplot(data = datax, aes(x = types, y = NormalizedValues, fill = genes)) + 
	geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
	scale_y_continuous(limits = c(0, ymax)) +
	# facet_wrap(. ~ genes) +
	labs(title = args[1], caption = date()) + 
	theme(axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "boxplot.png", sep = "."), width = 400, height = 600)
print(myplot)
dev.off()

tumors <- datax[!grep("Normal", types),]
tumors.upper <- tumors[ NormalizedValues > quantile(NormalizedValues, probs = 0.8),]
tumors.lower <- tumors[ NormalizedValues < quantile(NormalizedValues, probs = 0.2),]

write.table(tumors.upper, row.names = F, file = paste(args[1], "tumors.upper", "tsv", sep = "."), sep = "\t")
write.table(tumors.lower, row.names = F, file = paste(args[1], "tumors.lower", "tsv", sep = "."), sep = "\t")