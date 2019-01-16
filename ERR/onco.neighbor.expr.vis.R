# setwd("/Volumes/GoogleDrive/My Drive/wenbo.paper.check/ICGC.WXS.WGS.CTCF.expr/final.files")
library(data.table)
library(ggpubr)
files <- list.files()
files <- files[grep(".pMUT.CTCF-mut.txt",files)]

expr.all <- data.frame()
for(i in 1:length(files)){
    a <- strsplit(files[i], split = "\\.")[[1]]
    proj <- paste(a[2], a[3], sep = ".")
    exprs <- fread(files[i], header = F, sep = "\t", na.strings = "/")
    exprs[, meanv := mean(V3), by = V2]
    exprs[, sdv := sd(V3), by = V2]
    exprs[, zscore := (V3 - meanv) / sdv]
    exprs[, project := proj]
    exprs[!is.na(V4), Type:="OncoNeighborCTCFMut"]
    expr.all <- rbind(expr.all, exprs[,c(1:2,5,8:10)])
}
expr.all[!is.na(Type) & zscore>2, labels:=paste(V5, V2, sep = " To ")]
# expr.all[!is.na(Type)][V2=="TERT"]
# V1   V2      V5    zscore         project                Type
# 1: SP201610 TERT CLPTM1L -0.148637 exp_seq.PBCA-US OncoNeighborCTCFMut

ymin <- quantile(expr.all$zscore, probs = 0.01, na.rm = T)
ymax <- quantile(expr.all$zscore, probs = 0.99, na.rm = T)
p <- ggboxplot(data = expr.all[!is.na(Type)], x = "project", y = "zscore", 
               color = "project", fill = "project",
               outlier.shape = NA, alpha = 0.2, 
               xlab = "", ylab = "Oncogene expression",
               size = 0.6, width = 0.6, add = c("none"))
myplot <- ggpar(p, legend = "none", x.text.angle = 60)
print(myplot)

p <- ggdotplot(data = expr.all[!is.na(Type)], x = "project", y = "zscore", 
               color = "project", fill = "project",
               xlab = "", ylab = "Oncogene expression",
               label = "labels", repel = T,
               binwidth = 0.1, size = 0.3, add = c("mean_sd"))
myplot <- ggpar(p, legend = "none", x.text.angle = 60)
pdf(file = "ICGC.exp_seq.onconeighbors.CTCF-mut.zscore.pdf", width = 12, height = 8)
print(myplot)
dev.off()

count <- table(expr.all[!is.na(Type)]$project)
p <- ggdotplot(data = expr.all[!is.na(Type)][project %in% names(count[count>10])], x = "project", y = "zscore", 
               color = "project", fill = "project",
               xlab = "", ylab = "Oncogene expression",
               label = "labels", repel = T,
               binwidth = 0.1, size = 0.3, add = c("mean_sd"))
myplot <- ggpar(p, legend = "none", x.text.angle = 60)
pdf(file = "ICGC.exp_seq.onconeighbors.CTCF-mut.zscore.sim.pdf", width = 12, height = 8)
print(myplot)
dev.off()

######################################
# for sv
# setwd("~/Google Drive File Stream/My Drive/wenbo.paper.check/ICGC.SV")
library(data.table)
library(ggpubr)
expr.all <- fread("oncogene.affected.deletion.tsv", header = F, sep = "\t", na.strings = "/")
expr.all[, meanv := mean(V4), by = list(V1,V3)]
expr.all[, sdv := sd(V4), by = list(V1,V3)]
expr.all[, zscore := (V4 - meanv) / sdv]
expr.all[!is.na(V5), Type:="OncoNeighborPromDel"]
expr.all[!is.na(V5) & zscore > 2, labels:=V3]
p <- ggviolin(data = expr.all, x = "V1", y = "zscore", 
               color = "Type", fill = "Type",
               xlab = "", ylab = "Oncogene expression",
               add = c("mean")) +
    geom_hline(yintercept = 2, linetype = 2)
myplot <- ggpar(p, legend = "none")
pdf(file = "ICGC.onconeighbors.all.zscore.pdf", width = 15, height = 6)
print(myplot)
dev.off()

p <- ggdotplot(data = expr.all[!is.na(Type)], x = "V1", y = "zscore", 
               color = "V1", fill = "V1",
               xlab = "", ylab = "Oncogene expression",
               label = "labels", repel = T,
               binwidth = 0.3, size = 0.6, add = c("mean")) +
    geom_hline(yintercept = 2, linetype = 2)
myplot <- ggpar(p, legend = "none")
pdf(file = "ICGC.onconeighbors.deletion.zscore.pdf", width = 8, height = 6)
print(myplot)
dev.off()

# sort(table(expr.all[!is.na(V5) & zscore > 1]$V3))
mygene <- "CCNE1"
mytype <- "OV-AU"
mytypes <- c("OV-AU", "PACA-CA")
p <- ggdotplot(expr.all[V3 == mygene & V1 %in% mytypes], x = "V1", y = "V4", 
               color = "Type", fill = "Type",
               xlab = "", ylab = mygene,
               binwidth = 0.6, size = 2, add = c("none"))
myplot <- ggpar(p, legend = "none")
myplot

mylist <- fread("sim.myselected.txt", header = F, sep = "\t", na.strings = "/")
myexpr <- expr.all[V1 %in% mylist$V1 & V3 %in% mylist$V3]
genes <- unique(myexpr$V3)
for(i in 1:length(genes)){
    ymax <- quantile(myexpr[V3==genes[i]]$V4, probs = 0.995, na.rm = T)
    p <- ggdotplot(myexpr[V3==genes[i]], x = "V1", y = "V4", 
                   color = "Type", fill = "Type",
                   xlab = "", ylab = genes[i],
                   # label = "V5", repel = T,
                   # stackdir = "center", binpositions = "all", stackgroups = TRUE, 
                   binwidth = 0.1, size = 0.2, add = c("none"))
    myplot <- ggpar(p, legend = "none", ylim = c(0, ymax))
    # myplot
    pdf(file = paste("ICGC.onconeighbors.deletion.auto", genes[i], "pdf", sep = "."), width = 4, height = 6)
    print(myplot)
    dev.off()
}

for(i in 2:length(genes)){
    ymax <- quantile(myexpr[V3==genes[i] & V1 != "BRCA-FR"]$V4, probs = 0.995, na.rm = T)
    p <- ggdotplot(myexpr[V3==genes[i] & V1 != "BRCA-FR"], x = "V1", y = "V4", 
                   color = "Type", fill = "Type",
                   xlab = "", ylab = genes[i],
                   # label = "V5", repel = T,
                   # stackdir = "center", binpositions = "all", stackgroups = TRUE, 
                   binwidth = 0.1, size = 0.2, add = c("none"))
    myplot <- ggpar(p, legend = "none", ylim = c(0, ymax))
    # myplot
    pdf(file = paste("ICGC.onconeighbors.deletion", genes[i], "pdf", sep = "."), width = 4, height = 6)
    print(myplot)
    dev.off()
}

for(i in 1){
    ymax <- quantile(myexpr[V3==genes[i] & V1 == "BRCA-FR"]$V4, probs = 0.995, na.rm = T)
    p <- ggdotplot(myexpr[V3==genes[i] & V1 == "BRCA-FR"], x = "V1", y = "V4", 
                   color = "Type", fill = "Type",
                   xlab = "", ylab = genes[i],
                   # label = "V5", repel = T,
                   # stackdir = "center", binpositions = "all", stackgroups = TRUE, 
                   binwidth = 0.2, size = 0.5, add = c("none"))
    myplot <- ggpar(p, legend = "none", ylim = c(0, ymax))
    # myplot
    pdf(file = paste("ICGC.onconeighbors.deletion", genes[i], "pdf", sep = "."), width = 3, height = 6)
    print(myplot)
    dev.off()
}

######################################
# for sv
# setwd("~/Google Drive File Stream/My Drive/wenbo.paper.check/ICGC.SV")
library(data.table)
library(ggpubr)
expr.all <- fread("oncogene.affected.deletion.ctrl.tsv", header = F, sep = "\t", na.strings = "/")
expr.all[, meanv := mean(V4), by = list(V1,V3)]
expr.all[, sdv := sd(V4), by = list(V1,V3)]
expr.all[, zscore := (V4 - meanv) / sdv]
expr.all[!is.na(V5), Type:="OncoNeighborPromDel"]
expr.all[!is.na(V6), Type:="OncoNeighborOtherDel"]
expr.all[!is.na(V5) & zscore > 2, labels:=V3]
p <- ggviolin(data = expr.all, x = "V1", y = "zscore", 
              color = "Type", fill = "Type",
              xlab = "", ylab = "Oncogene expression",
              add = c("none")) +
    geom_hline(yintercept = 2, linetype = 2)
myplot <- ggpar(p, legend = "top")
pdf(file = "ICGC.onconeighbors.ctrl.all.zscore.pdf", width = 15, height = 6)
print(myplot)
dev.off()

p <- ggdotplot(data = expr.all[!is.na(Type)], x = "V1", y = "zscore", 
               color = "Type", fill = "Type",
               xlab = "", ylab = "Oncogene expression",
               label = "labels", repel = T,
               binwidth = 0.3, size = 0.6, add = c("mean")) +
    geom_hline(yintercept = 2, linetype = 2)
myplot <- ggpar(p, legend = "top")
pdf(file = "ICGC.onconeighbors.ctrl.deletion.zscore.pdf", width = 8, height = 6)
print(myplot)
dev.off()

# sort(table(expr.all[!is.na(V5) & zscore > 1]$V3))
# x <- table(expr.all[!is.na(Type),c(3,10)])
# x[x[,1] > 0 & x[,2] > 0,]
mygene <- "KDM5A"
expr.all[V3==mygene & !is.na(Type)]
# mytype <- "PACA-CA"
mytypes <- c("PACA-CA", "LIRI-JP")
p <- ggdotplot(expr.all[V3 == mygene & V1 %in% mytypes], x = "V1", y = "V4", 
               color = "Type", fill = "Type",
               xlab = "", ylab = mygene,
               binwidth = 0.6, size = 1, add = c("none"))
myplot <- ggpar(p, legend = "top")
myplot
