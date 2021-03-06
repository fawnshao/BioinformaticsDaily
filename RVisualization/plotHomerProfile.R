args <- commandArgs(TRUE)
# args <- c("Nup53blrp_E2_peaks.ATAC.profile.coverage", "Nup53blrp_E2_peaks.ATAC.profile.mat")
# args <- c("Nup53blrp_E2_peaks.ATAC.profile.coverage", "Nup53blrp_E2_peaks.srtbylogFC.ATAC.profile.mat")
# args <- c("2k.Nup53blrp_E2_peaks.ATAC.profile.coverage", "2k.Nup53blrp_E2_peaks.srtbylogFC.ATAC.profile.mat")
# args <- c("abc", "Nup53blrp_E2.Doxornot.ATAC.profile.mat")

library(ggplot2)
library(data.table)
library(reshape2)
library(pheatmap)

ggplotfile <- args[1]
pheatmapfile <- args[2]

ggplotdata <- fread(ggplotfile, sep = "\t")
colnames(ggplotdata)[1] <- "distance"
plotdata <- melt(ggplotdata, id.vars = "distance")
colnames(plotdata)[2:3] <- c("experiment", "signal")
p0 <- ggplot(plotdata, aes(x = distance, y = signal, colour = experiment)) + 
    geom_line() +
    theme(legend.position = "right", legend.box = "vertical") + 
    labs(title = "Lines plot", 
         subtitle = args[1],
         caption = date(),
         x = "Position")
pdf(file = paste(ggplotfile, "lines.pdf", sep = "."), width = 15, height = 8)
print(p0)
dev.off()

pheatmapdata <- fread(pheatmapfile, sep = "\t", header = T)
colors <- colorRampPalette(c("white", "red"))(100)

# pdf(file = paste(pheatmapfile, "pheatmap.pdf", sep = "."), width = 25, height = 12)
# myplot <- pheatmap(pheatmapdata[,-1], scale = "none",
# 	show_rownames = F, show_colnames = F, color = colors,
# 	cluster_cols = F, cluster_rows = F)
# dev.off()

x <- log2(data.matrix(pheatmapdata[,-1]) + 1)
# x <- data.matrix(pheatmapdata[,-1])
x[x < 0] <- 0
png(filename = paste(pheatmapfile, "pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(x, scale = "none",
	show_rownames = F, show_colnames = F, color = colors,
	cluster_cols = F, cluster_rows = F)
dev.off()

groups <- ncol(x)/241
for(i in 1:groups){
	png(filename = paste(pheatmapfile, i, "pheatmap.png", sep = "."), width = 1500, height = 1200)
	myplot <- pheatmap(x[,(i - 1) * 241 + c(1:241)], scale = "none",
		show_rownames = F, show_colnames = F, color = colors,
		cluster_cols = F, cluster_rows = F)
	dev.off()
}

