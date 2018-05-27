library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
datax <- fread(args[1], sep = "\t", header = F)
colnames(datax) <- c("groups", "classes", "values")

myplot <- ggplot(data = datax, aes(x = classes, y = values, fill = groups)) + 
	geom_bar(stat = "identity", position = "dodge") + 
	labs(title = args[1], caption = date()) + 
	theme(axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "barplot.png", sep = "."), width = 1200, height = 1000)
print(myplot)
dev.off()

