args <- commandArgs(TRUE)
# args <- c("NS15-Xiaoyu-MCF7-nupTT.rpkm.txt", "1.5")
library(ggplot2)

inputfile <- args[1]
foldchange.cutoff <- as.numeric(args[2])
inputdata <- read.table(inputfile, header = T, sep = "\t", row.names = 1)
dim(inputdata)
gene.max <- apply(inputdata, 1, max)
inputdata <- inputdata[gene.max > 3,]
dim(inputdata)
foldchange <- (inputdata[,2] + 1) / (inputdata[,1] + 1)
classes <- rep("Unclassfied", nrow(inputdata))
classes[foldchange > foldchange.cutoff] <- "Up"
classes[foldchange < 1 / foldchange.cutoff] <- "Down"
table(classes)
data2plot <- data.frame(inputdata, foldchange, classes)
rownames(data2plot) <- rownames(inputdata)

myplot <- ggplot(data = data2plot, aes(x = log2(NS15.Xiaoyu.tet.NUP53.2.E2 + 1), 
	y = log2(NS15.Xiaoyu.tet.NUP93.1.DOX.NUP53.E2 + 1), color = classes)) + 
	geom_point(aes(size = foldchange)) + #geom_smooth(method = lm) +
	labs(title = args[1], caption = date()) + 
pdf(file = paste(args[1], "scatterplot.log2.pdf", sep = "."), width = 12, height = 10)
print(myplot)
dev.off()

myplot <- ggplot(data = data2plot[classes!="Unclassfied",], aes(x = log2(NS15.Xiaoyu.tet.NUP53.2.E2 + 1), 
	y = log2(NS15.Xiaoyu.tet.NUP93.1.DOX.NUP53.E2 + 1), color = classes)) + 
	geom_point(aes(size = foldchange)) + #geom_smooth(method = lm) +
	labs(title = args[1], caption = date()) + 
pdf(file = paste(args[1], "scatterplot.log2.onlydiff.pdf", sep = "."), width = 12, height = 10)
print(myplot)
dev.off()

outdata <- data.frame(rownames(inputdata), data2plot)
write.table(outdata, file = paste(args[1], "ratios.xls", sep = "."), row.names = F, sep = "\t")
