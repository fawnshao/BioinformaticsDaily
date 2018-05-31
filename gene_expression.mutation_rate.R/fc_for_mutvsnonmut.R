library(data.table)
library(pheatmap)
args <- commandArgs(TRUE)
# args <- c("hkg.tsg.srtbyPCA.TCGA.CHD7.mut.sample")
input <- fread(args[1], sep = "\t", header = T, na.strings = "NA")
scores <- data.matrix(input[,-c(1:3)])
rownames(scores) <- as.matrix(input[,1])
count <- ncol(scores) / 2
scores.left <- scores[,1:count]
scores.right <- scores[,(count + 1):ncol(scores)]
mean.left <- apply(scores.left, 1, function(x) {mean(x, na.rm = T)})
mean.right <- apply(scores.right, 1, function(x) {mean(x, na.rm = T)})
mean.fc <- mean.left - mean.right
write.table(scores[abs(mean.fc) > 1,], file = paste(args[1],"fc1.tsv", sep = "."), sep = "\t")

colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(scores[abs(mean.fc) > 1,], scale = "none",
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()
