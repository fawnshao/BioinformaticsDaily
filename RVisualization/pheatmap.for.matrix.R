library(pheatmap)
library(data.table)
args <- commandArgs(TRUE)
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")

Plotheatmap <- function(x, y = annos, z = colors, pre = args[1]){
	png(filename = paste(pre, "pheatmap.png", sep = "."), width = 1000, height = 800)
	# pdf(file = paste(pre, "pheatmap.pdf", sep = "."), width = 12, height = 12)
	myplot <- pheatmap(x, scale = "none", annotation_row = y,
		show_rownames = F, show_colnames = T, color = z,
		cluster_cols = F, cluster_rows = F)
	dev.off()
}
Plotheatmap.cluster <- function(x, y = annos, z = colors, pre = args[1]){
	png(filename = paste(pre, "pheatmap.cluster.png", sep = "."), width = 1000, height = 800)
	# pdf(file = paste(pre, "pheatmap.pdf", sep = "."), width = 12, height = 12)
	myplot <- pheatmap(x, scale = "none", annotation_row = y,
		show_rownames = F, show_colnames = T, color = z,
		cluster_cols = T, cluster_rows = T)
	dev.off()
}
Plotheatmap.cluster.withname <- function(x, y = annos, z = colors, pre = args[1]){
	png(filename = paste(pre, "pheatmap.cluster.png", sep = "."), width = 1000, height = 1800)
	# pdf(file = paste(pre, "pheatmap.pdf", sep = "."), width = 10, height = 16)
	myplot <- pheatmap(x, scale = "none", annotation_row = y,
		show_rownames = T, show_colnames = T, color = z, fontsize_row = 1,
		cluster_cols = T, cluster_rows = T)
	dev.off()
}

scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
class[1:2647] <- "HKG"
types <- factor(class)
rownames(scores) <- paste(1:nrow(scores), genes, sep = ": ")

annos <- data.frame(class = types)
rownames(annos) <- rownames(scores)
colors <- colorRampPalette(c("white", "blue"))(10)

data.x <- scores
data.x[is.na(data.x)] <- 0
data.x[data.x > 1] <- 1
Plotheatmap(data.x, pre = args[1])
Plotheatmap.cluster(data.x, pre = args[1])

Plotheatmap.cluster(data.x[1:2647,], pre = paste(args[1], "HKG", sep = "."))
Plotheatmap.cluster(data.x[2648:3051,], pre = paste(args[1], "nonHKG", sep = "."))

Plotheatmap.cluster.withname(data.x[1:2647,], pre = paste(args[1], "HKG", sep = "."))
Plotheatmap.cluster.withname(data.x[2648:3051,], pre = paste(args[1], "nonHKG", sep = "."))

