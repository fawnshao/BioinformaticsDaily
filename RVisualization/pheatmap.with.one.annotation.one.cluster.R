library(data.table)
library(pheatmap)
args <- commandArgs(TRUE)
# args <- c("hkg.tsg.srtbyPCA.TCGA.NRF1.cat.sample", "hkg.tsg.srtbyPCA.TCGA.NRF1.cat.sample.annotation")
input <- fread(args[1], sep = "\t", header = T, na.strings = "NA")
class <- fread(args[2], sep = "\t", header = T)
scores <- data.matrix(input[,-c(1:3)])
rownames(scores) <- 1:nrow(scores) # as.matrix(input[,1])
colnames(scores) <- 1:ncol(scores)
annosC <- class[,2]
rownames(annosC) <- colnames(scores)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
scores[!is.na(scores) & scores > 20] <- 20
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(scores, scale = "none", annotation_col = annosC,
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = T)
dev.off()