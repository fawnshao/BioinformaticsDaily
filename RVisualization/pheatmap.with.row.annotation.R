library(data.table)
library(pheatmap)
args <- commandArgs(TRUE)
# args <- c("srt.sim.GTEx.pancancer.tsv", "srt.sim.hkg.tsg.annotation.2")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
class <- fread(args[2], sep = "\t", header = T, na.strings = "/")
scores <- data.matrix(input[,-1])
rownames(scores) <- 1:nrow(scores) # as.matrix(input[,1])
annosR <- class[,-1]
rownames(annosR) <- rownames(scores)
#####
# annosR[grep("hkg", annosR$GTEx), 1] <- "HKG"
# annosR[grep("hkg", annosR$pancancer), 2] <- "HKG"
# annosR[grep(",", annosR$pancancer), 2] <- "mixTSG"
scores[is.na(scores) | scores < 0] <- 0
#####
colors <- colorRampPalette(c("blue", "white", "red"))(100)
scores[scores > 10] <- 10
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 3000)
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 15, height = 16)
myplot <- pheatmap(scores, scale = "none", annotation_row = annosR,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()

annosR.e <- annosR[rowSums(scores)>0,]
rownames(annosR.e) <- rownames(scores[rowSums(scores)>0,])
png(filename = paste(args[1], "pheatmap.expressed.png", sep = "."), width = 2500, height = 3000)
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 15, height = 16)
myplot <- pheatmap(scores[rowSums(scores)>0,], scale = "none", annotation_row = annosR.e,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()
png(filename = paste(args[1], "pheatmap.expressed.cluster.png", sep = "."), width = 2500, height = 3000)
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 15, height = 16)
myplot <- pheatmap(scores[rowSums(scores)>0,], scale = "none", annotation_row = annosR.e,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = F)
dev.off()
write.table(data.frame(as.matrix(input[rowSums(scores)>0,1]),
	annosR.e, scores[rowSums(scores)>0,]),
	file = paste(args[1], "expressed.tsv", sep = "."), 
	quote = F, na = "/", row.names = F, sep = "\t")


annosR.g <- annosR[!is.na(annosR$GTEx) | !is.na(annosR$pancancer),]
rownames(annosR.g) <- rownames(scores[!is.na(annosR$GTEx) | !is.na(annosR$pancancer),])
png(filename = paste(args[1], "pheatmap.hkg.tsg.png", sep = "."), width = 2500, height = 3000)
# pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 15, height = 16)
myplot <- pheatmap(scores[!is.na(annosR$GTEx) | !is.na(annosR$pancancer),], 
	scale = "none", annotation_row = annosR.g,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()

if(args[3] == "Y"){
	png(filename = paste(args[1], "pheatmap.colcluster.png", sep = "."), width = 2500, height = 3000)
	myplot <- pheatmap(scores, scale = "none", annotation_row = annosR,
		show_rownames = F, show_colnames = T, color = colors, 
		cluster_cols = T, cluster_rows = F)
	dev.off()
}