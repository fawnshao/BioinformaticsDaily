library(data.table)
library(pheatmap)
args <- commandArgs(TRUE)
input <- fread(args[1], sep = "\t", header = T, na.strings = "NA")
class <- fread(args[2], sep = "\t", header = T, na.strings = "NA")
scores <- data.matrix(input[,-1])
rownames(scores) <- 1:nrow(scores) # as.matrix(input[,1])
annosR <- input[,-1]
rownames(annosR) <- rownames(scores)
#####
annosR[grep("hkg", GTEx), 1] <- "HKG"
annosR[grep("hkg", pancancer), 2] <- "HKG"
annosR[grep(",", pancancer), 2] <- "mixTSG"
scores[is.na(scores) | scores < 0] <- 0
#####
colors <- colorRampPalette(c("blue", "white", "red"))(100)
scores[scores > 10] <- 10
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(scores, scale = "none", annotation_row = annosR,
	show_rownames = F, show_colnames = F, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()