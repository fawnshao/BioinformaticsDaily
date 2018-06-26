# source("https://bioconductor.org/biocLite.R")
# biocLite("RUVSeq")
library(RUVSeq)
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
args <- c("GSE98476_salmon.genes.counts.txt")
input <- read.table(args[1], header = T, sep = "\t", row.names = 1)
# filter <- apply(input, 1, function(x) length(x[x>5])>=2)
# filtered <- zfGenes[filter,]
# genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
# spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
genes <- rownames(input)[grep("^ERCC-", rownames(input), invert = TRUE)]
spikes <- rownames(input)[grep("^ERCC-", rownames(input))]
# GSM2597190	HEP10_1008_rnaseq_1
# GSM2597191	HEP14_0079_rnaseq_1
# GSM2597192	HEP14_0080_rnaseq_1
# GSM2597193	HEP14_0120_rnaseq_1
# GSM2597193	LC_YY1-mutated (p.Lys179*)
# GSM2597192	LC_YY1-mutated (p.Leu366Pro)
# GSM2597190	LC_YY1-mutated (p.Asp380Tyr)
# GSM2597191	LC_YY1 hemizygous deletion
# x <- as.factor(c(rep("ctrl", 2), c("YY1del", "YY1mut", "YY1mut", "YY1del"), rep("ctrl", 4)))
# x <- as.factor(c(rep("ctrl", 2), c("YY1mut", "YY1mut", "YY1mut", "YY1mut"), rep("ctrl", 4)))
x <- as.factor(c(rep("ctrl", 6), rep("YY1del", 2)))
input <- input[,c(1:2,7:10,3,6)]
set <- newSeqExpressionSet(as.matrix(round(input,0)), phenoData = data.frame(x, row.names = colnames(input)))
# set

set <- betweenLaneNormalization(set, which = "upper")
set1 <- RUVg(set, spikes, k = 1)
pData(set1)
# pdf(file = paste(args,"RUV.plotRLE.pdf", sep = "."))
# plotRLE(set1, outline = FALSE, ylim = c(-4, 4), col = colors[x])
# dev.off()
# pdf(file = paste(args,"RUV.plotPCA.pdf", sep = "."))
# plotPCA(set1, col = colors[x], cex = 1.2)
# dev.off()
nc <- normCounts(set1)
write.table(x = data.frame(rownames(nc), nc), file = paste(args,"RUV.normCounts.tsv", sep = "."), sep = "\t", row.names = F)


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                               colData = pData(set1),
                               design = ~ W_1 + x)
dds <- DESeq(dds)
res <- results(dds)
# res

out <- data.frame(counts(dds, normalized = F), counts(dds, normalized = T), res)
write.table(x = data.frame(rownames(out), out), file = paste(args,"RUV.tsv", sep = "."), sep = "\t", row.names = F)



# design <- model.matrix(~x + W_1, data = pData(set1))
# y <- DGEList(counts = counts(set1), group = x)
# y <- calcNormFactors(y, method = "upperquartile")
# y <- estimateGLMCommonDisp(y, design)
# y <- estimateGLMTagwiseDisp(y, design)
# fit <- glmFit(y, design)
# lrt <- glmLRT(fit, coef = 2)
# topTags(lrt)
