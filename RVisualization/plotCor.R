setwd("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/all.HCT116.NUP93.PROseq/PRO.bothTreatments.bothReps/proseq.reproducibility/")
library(data.table)
library(corrplot)
library(Hmisc)
library(ggpubr)
input <- fread("HCT116.gene.featureCounts.fpkm", header = T, sep = "\t")
# input.vars <- data.matrix(input[,c(7:14)])
input.vars <- data.matrix(input[,c(7,9,11,13)])
input.vars.sum <- apply(input.vars, 1, sum)
# [input.vars.sum > 20, ]
mydata.cor <- cor(input.vars, method = c("pearson"))
# mydata.cor2 <- rcorr(input.vars)
corrplot(mydata.cor, method = "number", order = "hclust", cl.lim = range(mydata.cor))

p1 <- ggscatter(input, x = "NS25_noIAA_noTNF", y = "NS29_noIAA_noTNF",
                add = "reg.line", cor.coef = TRUE, 
                color = "tomato", fill = "tomato", size = 0.6,
                xlab = "Ctrl Rep1", ylab = "Ctrl Rep2")
p11 <- ggpar(p1, xscale = "log2", yscale = "log2")
ggsave(filename = "ggscatter.NUP93.UT.noTNF.reps.pdf", width = 8, height = 8, p11)
