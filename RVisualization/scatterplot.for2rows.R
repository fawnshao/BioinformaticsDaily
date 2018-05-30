library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
datax <- fread(args[1], sep = "\t", header = F)
x <- as.data.frame(t(datax[,-1]))
colnames(x) <- as.matrix(datax[,1])
myplot <- ggplot(data = x, aes(x = x[,1], y = x[,2])) + 
	geom_point() + geom_smooth(method = lm) +
	labs(title = args[1], caption = date()) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "scatterplot.png", sep = "."), width = 800, height = 800)
print(myplot)
dev.off()