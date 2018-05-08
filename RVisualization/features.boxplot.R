library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
input <- fread(args[1], sep = "\t", header = T)

for (i in 3:ncol(input)){
	png(filename = paste(args[1], colnames(input)[i], "boxplot.png", sep = "."), width = 1000, height = 600)
	datax <- input[,c(2,i),with = FALSE]
	colnames(datax) <- c("type", "value")
	ymax <- quantile(datax$value,probs = 0.95)
	myplot <- ggplot(data = datax, aes(x = type, y = value, fill = type)) + 
		geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
		scale_y_continuous(limits = c(0, ymax)) +
		ggtitle(colnames(input)[i]) + 
		theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
	print(myplot)
	dev.off()
}

