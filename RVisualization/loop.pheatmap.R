library(pheatmap)
mydists <- c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
myclusts <- c("ward.D", "ward.D2", "complete", "average", "mcquitty", "median")
# colors <- colorRampPalette(c("blue", "white", "red"))(10)
colors <- colorRampPalette(c("blue", "white", "red"))(10)
for(i in 1:length(mydists)){
    for(j in 1:length(myclusts)){
        pdf(file = paste("all.CFP1.HCFC1.YY1.GABPA.NRF1.otherHKGTSG.byColumn", mydists[i], myclusts[j], "pheatmap.pdf", sep = "."), width = 8, height = 6)
        myplot <- pheatmap(mydatamat, scale = "none", annotation_row = rowannos, annotation_col = colannos, 
                           show_rownames = F, show_colnames = T, color = colors, #breaks = mybreaks, 
                           cluster_cols = T, cluster_rows = T, 
                           clustering_distance_cols = mydists[i], 
                           clustering_distance_rows = mydists[i], 
                           clustering_method = myclusts[j], 
                           fontsize_col = 5, fontsize_row = 2)
        dev.off()
    }
}
