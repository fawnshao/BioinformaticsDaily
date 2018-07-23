args <- commandArgs(TRUE)
# head -1 NS15-Xiaoyu-MCF7-TT.rpkm.txt | cut -f 9-12 | sed 's/.mTD FPKM//g' | awk '{print "Gene\t"$0}' > NS15-Xiaoyu-MCF7-memoryTT.rpkm.txt
# tail -n +2 NS15-Xiaoyu-MCF7-TT.rpkm.txt | cut -f 1,9-12 >> NS15-Xiaoyu-MCF7-memoryTT.rpkm.txt
# head -1 NS15-Xiaoyu-MCF7-TT.rpkm.txt | cut -f 13-14 | sed 's/.mTD FPKM//g' | awk '{print "Gene\t"$0}' > NS15-Xiaoyu-MCF7-nupTT.rpkm.txt
# tail -n +2 NS15-Xiaoyu-MCF7-TT.rpkm.txt | cut -f 1,13-14 >> NS15-Xiaoyu-MCF7-nupTT.rpkm.txt
# perl $myperl <(cut -f 1-6,8 NS15-Xiaoyu-MCF7-TT.rpkm.txt) <(sed 's/"//g' NS15-Xiaoyu-MCF7-memoryTT.rpkm.txt.ratios.xls) 0 0 | cut -f 1-9,11- > NS15-Xiaoyu-MCF7-memoryTT.rpkm.txt.ratios.annotated.xls
# args <- c("NS15-Xiaoyu-MCF7-nupTT.rpkm.txt", "1.5")
# args <- c("NS15-Xiaoyu-MCF7-memoryTT.rpkm.txt", "1.5")
library(ggplot2)

inputfile <- args[1]
foldchange.cutoff <- as.numeric(args[2])
inputdata <- read.table(inputfile, header = T, sep = "\t", row.names = 1)
dim(inputdata)
gene.max <- apply(inputdata, 1, max)
inputdata <- inputdata[gene.max > 3,]
e2.min <- apply(inputdata[,c(1,3)], 1, min)
# inputdata <- inputdata[e2.min > 2,]
dim(inputdata)
inputdata[inputdata < 0.1] <- 0.1
# d0.ratio <- (inputdata[,1] + 1) / (inputdata[,2] + 1)
# d4.ratio <- (inputdata[,3] + 1) / (inputdata[,4] + 1)
d0.ratio <- inputdata[,1] / inputdata[,2]
d4.ratio <- inputdata[,3] / inputdata[,4]
mem.ratio <- d4.ratio / d0.ratio
# basal.EtOH.ratio <- inputdata[,4] / inputdata[,2]
# basal.E2.ratio <- inputdata[,3] / inputdata[,1]
basal.ratio <- (inputdata[,3] + inputdata[,4]) / (inputdata[,1] + inputdata[,2])
classes <- rep("Unclassfied", nrow(inputdata))
# classes[d0.ratio > foldchange.cutoff & d4.ratio > foldchange.cutoff & mem.ratio > foldchange.cutoff & basal.E2.ratio > foldchange.cutoff] <- "Induced and Enhanced"
# classes[d0.ratio > foldchange.cutoff & d4.ratio > foldchange.cutoff & mem.ratio >= 1 / foldchange.cutoff & mem.ratio <= foldchange.cutoff & basal.E2.ratio >= 1 / foldchange.cutoff] <- "Induced and Memorized"
classes[d0.ratio > foldchange.cutoff & d4.ratio > foldchange.cutoff & mem.ratio >= 1 / foldchange.cutoff & e2.min > 2] <- "Induced and Memorized"
classes[d0.ratio > foldchange.cutoff & d4.ratio > foldchange.cutoff & mem.ratio < 1 / foldchange.cutoff & e2.min > 2] <- "Induced and Losing Memory"
classes[d0.ratio > foldchange.cutoff & d4.ratio < foldchange.cutoff & mem.ratio < 1 / foldchange.cutoff & e2.min > 2] <- "Induced but Forgot"
classes[d0.ratio <= foldchange.cutoff & d4.ratio > foldchange.cutoff & mem.ratio > foldchange.cutoff & basal.ratio > foldchange.cutoff & e2.min > 2] <- "Gradually Induced"
# classes[d0.ratio < 1 / foldchange.cutoff & d4.ratio < 1 / foldchange.cutoff] <- "Not Induced"
stats <- data.frame(table(classes))
data2plot <- data.frame(d0.ratio, d4.ratio, mem.ratio, classes)
rownames(data2plot) <- rownames(inputdata)

myplot <- ggplot(stats, aes(x = "", y = Freq, fill = classes)) + geom_bar(width = 1, stat = "identity") + coord_polar("y")
pdf(file = paste(args[1], "classes.pie.pdf", sep = "."), width = 12, height = 10)
print(myplot)
dev.off()

# myplot <- ggplot(data = data2plot, aes(x = d0.ratio, y = d4.ratio, color = classes)) + 
# 	geom_point() + geom_smooth(method = lm) +
# 	labs(title = args[1], caption = date()) + 
# 	coord_trans(x="log2", y="log2") + 
# 	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
# pdf(file = paste(args[1], "scatterplot.pdf", sep = "."), width = 8, height = 8)
# print(myplot)
# dev.off()

myplot <- ggplot(data = data2plot, aes(x = log2(d0.ratio), y = log2(d4.ratio), color = classes)) + 
	geom_point(aes(size = mem.ratio)) + #geom_smooth(method = lm) +
	labs(title = args[1], caption = date()) + 
	# theme(axis.text.x = element_text(angle = 60, hjust = 1))
	# theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
pdf(file = paste(args[1], "scatterplot.log2.pdf", sep = "."), width = 12, height = 10)
print(myplot)
dev.off()

myplot <- ggplot(data = data2plot[classes!="Unclassfied",], aes(x = log2(d0.ratio), y = log2(d4.ratio), color = classes)) + 
	geom_point(aes(size = mem.ratio)) + #geom_smooth(method = lm) +
	labs(title = args[1], caption = date()) + 
	# theme(axis.text.x = element_text(angle = 60, hjust = 1))
	# theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))
pdf(file = paste(args[1], "scatterplot.log2.onlydiff.pdf", sep = "."), width = 12, height = 10)
print(myplot)
dev.off()

outdata <- data.frame(rownames(inputdata), inputdata, data2plot)
write.table(outdata, file = paste(args[1], "ratios.xls", sep = "."), row.names = F, sep = "\t")
