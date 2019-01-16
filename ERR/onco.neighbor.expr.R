# setwd("~/Google Drive File Stream/My Drive/wenbo.paper.check/ICGC.WGS.mut2oncogene.expr")
# 1. different oncogene have different basal expression level, how can we combine them together => zvalue!
# 2. for oncogene x, neighbor1 is highly mutated, but neighbor2 is lowly mutated, how to label?
library(ggpubr)
library(data.table)
datax <- fread("13cancers.oncogene.200kb.nononcogene.neighbors.each.pmutPerMMutPerPromoter.mat", 
               sep = "\t", na.strings = "", header = T)
datay <- fread("oncogene.200kb.nononcogene.neighbors.gene2gene", header = F)
colnames(datax)[1] <- "Neighbors" 
datax[,NeighborTo:=datay[match(Neighbors,datay$V1)]$V2]
datamat <- melt(datax, id.vars = c("Neighbors", "NeighborTo"))
datamat[is.na(value), value:=0]
datamat[,OncoNeigborsMut:=mean(value), by = c("NeighborTo", "variable")]
types <- unique(datamat$variable)

# i < -3
for(i in c(3,4,7,9,10)){
    print(types[i])
    topmutcount <- quantile(unique(datamat[variable==types[i]]$OncoNeigborsMut), probs = 0.9, na.rm = F)
    topmutneighbors <- datamat[variable==types[i]][OncoNeigborsMut > topmutcount]
    lowmutcount <- quantile(unique(datamat[variable==types[i]]$OncoNeigborsMut), probs = 0.1, na.rm = F)
    lowmutneighbors <- datamat[variable==types[i]][OncoNeigborsMut < lowmutcount]
    headtails <- rbindlist(list(topmutneighbors, lowmutneighbors))
    # sort(table(topmutneighbors$NeighborTo))
    # or just sort by promoter mutation for each neighbors
    # topmutneighbors <- datamat[variable==types[i]][order(OncoNeigborsMut,decreasing = T)][1:10]
    
    exprs <- fread(paste("exp_seq", types[i], "pMUT.txt", sep = "."), header = F, na.strings = "/")
    pmuts <- fread(paste(types[i], "WGS.pMUT", sep = "."), header = F, na.strings = "/")
    exprs[, zscore:=(V3-mean(V3))/sd(V3), by = c("V2")]
    
    oncoexprs <- exprs[V2 %in% headtails$NeighborTo][!is.na(V2)]
    toponconeighborpmuts <- pmuts[V6 %in% topmutneighbors$Neighbors]
    lowonconeighborpmuts <- pmuts[V6 %in% lowmutneighbors$Neighbors]
    toponconeighborpmuts.list <- unique(paste(
        toponconeighborpmuts$V5, 
        topmutneighbors[match(toponconeighborpmuts$V6, topmutneighbors$Neighbors)]$NeighborTo,sep = "|"))
    lowonconeighborpmuts.list <- unique(paste(
        lowonconeighborpmuts$V5, 
        lowmutneighbors[match(lowonconeighborpmuts$V6, lowmutneighbors$Neighbors)]$NeighborTo,sep = "|"))
    oncoexprs[,NeighborPromMut:=ifelse(paste(V1,V2, sep = "|") %in% toponconeighborpmuts.list, 
                                       "NeighborsHighMut", "NeighborsLowMut")]
    # oncoexprs[is.na(NeighborPromMut), NeighborPromMut:=ifelse(paste(V1,V2, sep = "|") %in% lowonconeighborpmuts.list, 
    #                                    "NeighborsLowMut", NA)]
    oncoexprs[,OncogenePromMut:=ifelse(is.na(V4), "OncoPromoterNotMut", "OncoPromoterMut")]
    # oncoexprs[,log2Var:=log2(V3 + 1)]
    # p <- ggdotplot(data = oncoexprs, x = "NeighborPromMut", y = "zscore", 
    #                color = "OncogenePromMut", fill = "OncogenePromMut", 
    #                size = 0.2, binwidth = 0.1,
    #                add = c("boxplot", "mean"))
    ymin <- quantile(oncoexprs$zscore, probs = 0.01, na.rm = T)
    ymax <- quantile(oncoexprs$zscore, probs = 0.99, na.rm = T)
    p <- ggboxplot(data = oncoexprs, x = "NeighborPromMut", y = "zscore", 
                   color = "OncogenePromMut", fill = "OncogenePromMut",
                   outlier.shape = NA, alpha = 0.2, 
                   xlab = "", ylab = "Oncogene expression",
                   size = 0.6, width = 0.6, add = c("none"))
    # p
    myplot <- ggpar(p, ylim = c(ymin, ymax))
    pdf(file = paste("top10", types[i], "pdf", sep = "."), width = 5, height = 6)
    print(myplot)
    dev.off()
    
    p <- ggboxplot(data = oncoexprs, x = "NeighborPromMut", y = "zscore", 
                   color = "NeighborPromMut", fill = "NeighborPromMut",
                   outlier.shape = NA, alpha = 0.2, 
                   xlab = "", ylab = "Oncogene expression",
                   size = 0.6, width = 0.6, add = c("none"))
    # p
    myplot <- ggpar(p, ylim = c(ymin, ymax))
    pdf(file = paste("top210", types[i], "sim.pdf", sep = "."), width = 5, height = 6)
    print(myplot)
    dev.off()
}

# # ggpar(p, ylim = c(0, 20))
# # JAK3
# toponcogenes <- sort(table(topmutneighbors$NeighborTo), decreasing = T)
# j <- 5
# p <- ggdotplot(data = oncoexprs[V2==names(toponcogenes)[j]], x = "NeighborPromMut", y = "log2Var", 
#                color = "OncogenePromMut", fill = "OncogenePromMut", 
#                size = 0.2, binwidth = 0.1, xlab = names(toponcogenes)[j],
#                add = c("violin"))
# p

