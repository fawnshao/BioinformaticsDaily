# args <- commandArgs(TRUE)
library(ggpubr)
library(data.table)
library(reshape2)
library(pheatmap)
library(plyr)
args <- c("xena.mygene.txt", "xena.stage.info.txt", "143gene.txt")
exprmat <- fread(args[1], header = T)
samples <- fread(args[2], header = T)
genes <- fread(args[3], header = F)

# Code  Definition  Short Letter Code
# 01  Primary Solid Tumor TP
# 02  Recurrent Solid Tumor TR
# 03  Primary Blood Derived Cancer - Peripheral Blood TB
# 04  Recurrent Blood Derived Cancer - Bone Marrow  TRBM
# 05  Additional - New Primary  TAP
# 06  Metastatic  TM
# 07  Additional Metastatic TAM
# 08  Human Tumor Original Cells  THOC
# 09  Primary Blood Derived Cancer - Bone Marrow  TBM
# 10  Blood Derived Normal  NB
# 11  Solid Tissue Normal NT
# 12  Buccal Cell Normal  NBC
# 13  EBV Immortalized Normal NEBV
# 14  Bone Marrow Normal  NBM
# 15  sample type 15  15SH
# 16  sample type 16  16SH
# 20  Control Analyte CELLC
# 40  Recurrent Blood Derived Cancer - Peripheral Blood TRB
# 50  Cell Lines  CELL
# 60  Primary Xenograft Tissue  XP
# 61  Cell Line Derived Xenograft Tissue  XCL
# 99  sample type 99  99SH
exprmatsim <- data.matrix(exprmat[,-1])
rownames(exprmatsim) <- exprmat$sample
exprmatsim <- t(exprmatsim)
cancersample <- matrix(unlist(strsplit(rownames(exprmatsim), split = "-")), ncol = 4, byrow = T)[,4]
cancersample <- as.numeric(cancersample)
sampledesc <- data.frame(SampleTypeCodes = cancersample, 
  CancerType = samples[match(rownames(exprmatsim), samples$sample)]$`cancer type abbreviation`, 
  TumorStage = samples[match(rownames(exprmatsim), samples$sample)]$`ajcc_pathologic_tumor_stage`)
sampledesc <- data.frame(sampledesc, 
  SampleCodeTumorStage = paste(sampledesc$SampleTypeCodes, sampledesc$TumorStage))
rownames(sampledesc) <- paste(rownames(exprmatsim), 1:nrow(sampledesc), sep = ".")
rownames(exprmatsim) <- paste(rownames(exprmatsim), 1:nrow(exprmatsim), sep = ".")
allmat <- data.frame(sampledesc, exprmatsim)

ccount <- sort(table(sampledesc$CancerType), decreasing = T) 
# ct <- "BRCA"
for(i in c(1:16,18,20:26,28)){
    ct <- names(ccount)[i]
    tmpmat <- allmat[allmat$CancerType==ct & allmat$TumorStage!="" & allmat$TumorStage!="[Discrepancy]", ]
    if(nrow(tmpmat) > 50){
      tmpmat.srt <- tmpmat[order(tmpmat$SampleCodeTumorStage),]
      mysamples <- tmpmat.srt[tmpmat.srt$CancerType==ct, c(2,4)]
      mymat <- t(tmpmat.srt[tmpmat.srt$CancerType==ct, -c(1:4)])
      # nacount <- apply(mymat, 1, function(x){length(x[is.na(x)])})

      mybreaks <- c(seq(-20, -3, 1), seq(-2.9, 2.9, 0.1), seq(3, 20, 1))
      mycolors <- colorRampPalette(c("blue", "white", "red"))(length(mybreaks))
      png(filename = paste("pheatmap", ct, "png", sep = "."), width = 1500, height = 1000)
      myplot <- pheatmap(mymat, scale = "row", 
                         breaks = mybreaks, color = mycolors,
                         annotation_col = mysamples, show_colnames = F, 
                         main = paste(ct, "Left to Right: stage low to high"), 
                         cluster_cols = F, cluster_rows = T)
      dev.off()

      png(filename = paste("auto.pheatmap", ct, "png", sep = "."), width = 1500, height = 1000)
      myplot <- pheatmap(mymat, scale = "row", 
                          annotation_col = mysamples, show_colnames = F, 
                         main = paste(ct, "Left to Right: stage low to high"), 
                         cluster_cols = F, cluster_rows = T)
      dev.off()
    }
}

for(i in c(17,19,27,29:32)){
    ct <- names(ccount)[i]
    tmpmat <- allmat[allmat$CancerType==ct & allmat$TumorStage!="" & allmat$TumorStage!="[Discrepancy]", ]
    if(nrow(tmpmat) > 50){
      tmpmat.srt <- tmpmat[order(tmpmat$SampleCodeTumorStage),]
      mysamples <- tmpmat.srt[tmpmat.srt$CancerType==ct, c(2,4)]
      mymat <- t(tmpmat.srt[tmpmat.srt$CancerType==ct, -c(1:4)])

      mybreaks <- c(seq(-20, -3, 1), seq(-2.9, 2.9, 0.1), seq(3, 20, 1))
      mycolors <- colorRampPalette(c("blue", "white", "red"))(length(mybreaks))
      png(filename = paste("pheatmap", ct, "png", sep = "."), width = 1500, height = 1000)
      myplot <- pheatmap(mymat, scale = "row", 
                         breaks = mybreaks, color = mycolors,
                         annotation_col = mysamples, show_colnames = F, 
                         main = paste(ct, "Left to Right: stage low to high"), 
                         cluster_cols = F, cluster_rows = F)
      dev.off()

      png(filename = paste("auto.pheatmap", ct, "png", sep = "."), width = 1500, height = 1000)
      myplot <- pheatmap(mymat, scale = "row", 
                          annotation_col = mysamples, show_colnames = F, 
                         main = paste(ct, "Left to Right: stage low to high"), 
                         cluster_cols = F, cluster_rows = F)
      dev.off()
    }
}

for(i in 1:32){
    ct <- names(ccount)[i]
    tmpmat <- allmat[allmat$CancerType==ct & allmat$TumorStage!="" & allmat$TumorStage!="[Discrepancy]", ]
    if(nrow(tmpmat) > 50){
      tmpmat.srt <- tmpmat[order(tmpmat$SampleCodeTumorStage),]
      mysamples <- tmpmat.srt[tmpmat.srt$CancerType==ct, c(2,4)]
      mymat <- t(tmpmat.srt[tmpmat.srt$CancerType==ct, -c(1:4)])

      mybreaks <- c(seq(-20, -3, 1), seq(-2.9, 2.9, 0.1), seq(3, 20, 1))
      mycolors <- colorRampPalette(c("blue", "white", "red"))(length(mybreaks))
      pdf(file = paste("byname.pheatmap", ct, "pdf", sep = "."), width = 15, height = 15)
      myplot <- pheatmap(mymat, scale = "row", 
                         breaks = mybreaks, color = mycolors,
                         annotation_col = mysamples, show_colnames = F, 
                         main = paste(ct, "Left to Right: stage low to high"), 
                         cluster_cols = F, cluster_rows = F)
      dev.off()

      pdf(file = paste("byname.auto.pheatmap", ct, "pdf", sep = "."), width = 15, height = 15)
      myplot <- pheatmap(mymat, scale = "row", 
                          annotation_col = mysamples, show_colnames = F, 
                         main = paste(ct, "Left to Right: stage low to high"), 
                         cluster_cols = F, cluster_rows = F)
      dev.off()
    }
}

ct="PAAD"
tmpmat <- allmat[allmat$CancerType==ct & allmat$TumorStage!="" & allmat$TumorStage!="[Discrepancy]", ]
tmpmat.srt <- tmpmat[order(tmpmat$SampleCodeTumorStage),]
mysamples <- tmpmat.srt[tmpmat.srt$CancerType==ct, c(2,4)]
mymat <- t(tmpmat.srt[tmpmat.srt$CancerType==ct, -c(1:4)])

mybreaks <- c(seq(-10, -3, 1), seq(-2.9, 2.9, 0.1), seq(3, 10, 1))
mycolors <- colorRampPalette(c("blue", "white", "red"))(length(mybreaks))
pdf(file = paste("pheatmap", ct, "pdf", sep = "."), width = 15, height = 15)
myplot <- pheatmap(mymat, scale = "row", 
                   breaks = mybreaks, color = mycolors,
                   annotation_col = mysamples, show_colnames = F, 
                   main = paste(ct, "Left to Right: stage low to high"), 
                   cluster_cols = F, cluster_rows = T)
dev.off()

pdf(file = paste("auto.pheatmap", ct, "pdf", sep = "."), width = 15, height = 15)
myplot <- pheatmap(mymat, scale = "row", 
                    annotation_col = mysamples, show_colnames = F, 
                   main = paste(ct, "Left to Right: stage low to high"), 
                   cluster_cols = F, cluster_rows = T)
dev.off()
