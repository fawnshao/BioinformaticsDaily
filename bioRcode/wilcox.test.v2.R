# wilcox.test(x = rnorm(100, 50, 5), y = rnorm(60, 30, 5))$p.value
library(data.table)
library(ggpubr)
args <- commandArgs(TRUE)
args <- c("onconeighbors.forPermutation.txt", "13cancers.allmut.stats", "13cancers.pMUT.stats")
totalpromoter <- 27502
datax <- fread(args[1], header = T, na.strings = "/")
datay <- fread(args[2], header = F, na.strings = "/")
dataz <- fread(args[3], header = F, na.strings = "/")
datax[,freq:=NeighborspMUTCount/NeighborspCount]
datax[,Group:=ifelse(Permutaion=="Oncogene", "Oncogene", "Permutation")]
types <- unique(datax$CancerType)
perms <- unique(datax$Permutaion)
for(i in 1:length(types)){
    xmax <- quantile(datax[CancerType==types[i]]$freq, probs = 0.99)
    p <- ggdensity(data = datax[CancerType==types[i]], 
                   color = "Group", fill = "Group", alpha = 0.3,
                   add = "mean", rug = TRUE,
                   x = "freq", y = "..density..",
                   xlab = paste(types[i], "averageMutationCountPerOncogene"))
    myplot <- ggpar(p, xlim = c(0, xmax))
    pdf(file = paste("ggdensity.onconeighbors.permutation", types[i], "pdf", sep = "."), width = 8, height = 6)
    print(myplot)
    dev.off()
}

medianvals <- matrix(nrow = length(perms), ncol = length(types))
rownames(medianvals) <- perms
colnames(medianvals) <- types
meanvals <- medianvals
for(i in 1:length(perms)){
    for(j in 1:length(types)){
        medianvals[i, j] <- median(datax[Permutaion==perms[i]][CancerType==types[j]]$freq, na.rm = T)
        meanvals[i, j] <- mean(datax[Permutaion==perms[i]][CancerType==types[j]]$freq, na.rm = T)
    }
}

mediantable <- as.data.table(melt(medianvals))
meantable <- as.data.table(melt(meanvals))
for(i in 1:length(types)){
    oncovals <- mediantable[Var2==types[i]][Var1=="Oncogene"]$value
    pvals <- nrow(mediantable[Var2==types[i]][Var1!="Oncogene"][value >= oncovals])/1000
    p <- gghistogram(mediantable[Var2==types[i]][Var1!="Oncogene"], x = "value", y = "..count..", bins = 100, 
                     color = "#00AFBB", fill = "#00AFBB", alpha = 0.6, xlab = "") + 
        geom_vline(xintercept = oncovals, 
                   color = "#FC4E07", linetype = "dashed", size = 0.5) +
        annotate("text", label = paste("p =", pvals), 
                 x = max(mediantable[Var2==types[i]][Var1!="Oncogene"]$value) * 0.9,
                 y = 50)
    pdf(file = paste("gghistogram.onconeighbors.permutation.median", types[i], "pdf", sep = "."), width = 8, height = 6)
    print(p)
    dev.off()
    oncovals <- meantable[Var2==types[i]][Var1=="Oncogene"]$value
    pvals <- nrow(meantable[Var2==types[i]][Var1!="Oncogene"][value >= oncovals])/1000
    p <- gghistogram(meantable[Var2==types[i]][Var1!="Oncogene"], x = "value", y = "..count..", bins = 100, 
                     color = "#00AFBB", fill = "#00AFBB", alpha = 0.6, xlab = "") + 
        geom_vline(xintercept = oncovals, 
                   color = "#FC4E07", linetype = "dashed", size = 0.5) +
        annotate("text", label = paste("p =", pvals), 
                 x = max(meantable[Var2==types[i]][Var1!="Oncogene"]$value) * 0.9,
                 y = 50)
    pdf(file = paste("gghistogram.onconeighbors.permutation.mean", types[i], "pdf", sep = "."), width = 8, height = 6)
    print(p)
    dev.off()
}

# 
# 
# 
# res <- data.frame(Types = types, pvals = pvals)
# write.table(res, file = args[3], append = T,  sep = "\t", row.names = F, col.names = F)
# 
# 
# summary(datax[CancerType==types[12]][Permutaion!="Oncogene"]$freq)
# summary(datax[CancerType==types[12]][Permutaion=="Oncogene"]$freq)
# summary(dataz$V2/totalpromoter)
# 
# ggdensity(data = datax[CancerType==types[12]], 
#           color = "Group", fill = "Group", alpha = 0.6,
#           x = "freq", y = "..density..")
# 
# p <- ggdensity(data = datax[CancerType!=types[12]], 
#           color = "Group", fill = "Group", alpha = 0.6,
#           x = "freq", y = "..density..")
# ggpar(facet(p, facet.by = "CancerType", ncol = 4))
