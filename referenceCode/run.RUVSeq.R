args <- commandArgs(TRUE)
library(RUVSeq)
# args <- c("attr2.HKGttseq.htseq.txt")
## ----readin-----------------------------------------------
myGenes <- read.table(args[1], header = TRUE, row.names = 1)
mtPCGs <- c("ENSG00000198695.2","ENSG00000198712.1","ENSG00000198727.2",
	"ENSG00000198763.3","ENSG00000198786.2","ENSG00000198804.2",
	"ENSG00000198840.2","ENSG00000198886.2","ENSG00000198888.2",
	"ENSG00000198899.2","ENSG00000198938.2","ENSG00000212907.2","ENSG00000228253.1")
## ----filter-----------------------------------------------
filter <- apply(myGenes, 1, function(x) length(x[x > 5]) >= 2)
filtered <- myGenes[filter,]
genes <- rownames(filtered)[grep("^ENSG", rownames(filtered))]
spikes <- mtPCGs
# spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

## ----store_data-------------------------------------------
# x <- as.factor(rep(c("Ctl", "Trt"), each = 2))
x <- as.factor(rep(c("Ctl", "Trt"), 2))
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names = colnames(filtered)))
# set

## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
set <- betweenLaneNormalization(set, which = "upper")

## ----ruv_spikes, fig.cap="RUVg normalization based on spike-in controls.", fig.subcap=c("RLE plot","PCA plot")----
set1 <- RUVg(set, spikes, k = 1)
pData(set1)
d <- attr(normCounts(set1),"scaled:scale")
# filtered[rownames(filtered) %in% spikes, ]
write.csv(d, file = paste(args[1], "RUVSeq.scalefactor.csv", sep = "."), quote = F)
# tst <- glm.fit(x = spikeinexpr[,1]/colSums(counts(set1))[1]*d[1], 
# 	y = spikeinexpr[,2]/colSums(counts(set1))[2]*d[2])
# coefficients(tst)
# #  1.045293
# tst <- glm.fit(x = spikeinexpr[,1]/colSums(counts(set1))[1]/d[1], 
# 	y = spikeinexpr[,2]/colSums(counts(set1))[2]/d[2])
# coefficients(tst)
# # 1.322967

# library(ggpubr)
# library(reshape2)
# library(data.table)
# # df <- data.table(CXXC1.ctl = spikeinexpr[,1]/colSums(counts(set1))[1]*d[1],
# # 	CXXC1.IAA = spikeinexpr[,2]/colSums(counts(set1))[2]*d[2])
# df <- data.table(CXXC1.ctl = spikeinexpr[,1]/colSums(counts(set1))[1]*d[1],
# 	CXXC1.IAA = spikeinexpr[,2]/colSums(counts(set1))[2]*d[2])
# p1 <- ggscatter(df, x = "CXXC1.ctl", y = "CXXC1.IAA", 
# 	color = "black", shape = 21, size = 3, 
# 	add = "reg.line",
# 	add.params = list(color = "blue", fill = "lightgray"), 
# 	conf.int = TRUE, 
# 	cor.coef = TRUE, 
# 	cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
# 	)
# ggpar(p1, xlim = c(0, max(df)), ylim = c(0, max(df)))

# tst <- glm.fit(x = spikeinexprnorm[,1]/colSums(normCounts(set1))[1]*d[1], 
# 	y = spikeinexprnorm[,2]/colSums(normCounts(set1))[2]*d[2])
# # coefficients(tst)
# 1.12353
# tst <- glm.fit(x = spikeinexprnorm[,1]/colSums(normCounts(set1))[1], 
# 	y = spikeinexprnorm[,2]/colSums(normCounts(set1))[2])
# coefficients(tst)
# # 1.26398
# tst <- glm.fit(x = spikeinexprnorm[,1], 
# 	y = spikeinexprnorm[,2])
# coefficients(tst)
# # 1.227767
# tst <- glm.fit(x = spikeinexprnorm[,1]*d[1], 
# 	y = spikeinexprnorm[,2]*d[2])
# coefficients(tst)
# # 1.091341

## ----DESeq-------------------------------------------
library(DESeq)
design <- model.matrix(~x + W_1, data = pData(set1))
y <- DGEList(counts = counts(set1), group = x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

## ----DESeq-------------------------------------------
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts(set1), colData = pData(set1), design = ~ W_1 + x)
dds <- DESeq(dds)
# resultsNames(dds)
# res <- results(dds, name = resultsNames(dds)[2])
res <- results(dds)
outputs <- data.frame(counts(set1), normCounts(set1), counts(dds, normalized = TRUE), as.data.frame(res))
write.csv(outputs, file = paste(args[1], "DESeq2_results.csv", sep = "."), quote = F)


## ----visualization-------------------------------------------
# spikeinexpr <- counts(set1)[rownames(counts(set1)) %in% spikes, ]
# spikeinexprnorm <- normCounts(set1)[rownames(normCounts(set1)) %in% spikes, ]
# spikeinexprdeg <- outputs[rownames(outputs) %in% spikes, ]

# tst <- glm.fit(x = spikeinexpr[,1], y = spikeinexpr[,2])
# coefficients(tst)
# # > coefficients(tst)
# # [1] 1.00035
# tst <- glm.fit(x = spikeinexprnorm[,1], y = spikeinexprnorm[,2])
# coefficients(tst)
# # > coefficients(tst)
# # [1] 1.227767
# tst <- glm.fit(x = spikeinexprnorm[,1]/colSums(normCounts(set1))[1], y = spikeinexprnorm[,2]/colSums(normCounts(set1))[2])
# coefficients(tst)
# # > coefficients(tst)
# # [1] 1.26398
# tst <- glm.fit(x = spikeinexpr[,1]/colSums(counts(set1))[1], y = spikeinexpr[,2]/colSums(counts(set1))[2])
# coefficients(tst)
# # > coefficients(tst)
# # [1] 1.175962
# d <- sizeFactors(dds)
# tst <- glm.fit(x = spikeinexprnorm[,1]*d[1], y = spikeinexprnorm[,2]*d[2])
# coefficients(tst)
# # > coefficients(tst)
# # [1] 1.091341
# tst <- glm.fit(x = spikeinexprdeg[,9]*d[1], y = spikeinexprdeg[,10]*d[2])
# coefficients(tst)
# # 1.00035
# tst <- glm.fit(x = spikeinexprdeg[,5]/d[1], y = spikeinexprdeg[,6]/d[2])
# coefficients(tst)
# # 1.381247
# tst <- glm.fit(x = spikeinexprdeg[,1]/d[1], y = spikeinexprdeg[,2]/d[2])
# coefficients(tst)
# # 1.125401
# tst <- glm.fit(x = spikeinexprdeg[,1]/1.15213, y = spikeinexprdeg[,2])
# coefficients(tst)
# # 1.152533
# tst <- glm.fit(x = spikeinexprdeg[,5]*1.15213, y = spikeinexprdeg[,6])
# coefficients(tst)
# # 1.06565
# dd <- attr(normCounts(set1),"scaled:scale")
# tst <- glm.fit(x = spikeinexprnorm[,1]/dd[1], y = spikeinexprnorm[,2]/dd[2])
# coefficients(tst)
# # 1.414547
# tst <- glm.fit(x = spikeinexprnorm[,1]*dd[1], y = spikeinexprnorm[,2]*dd[2])
# coefficients(tst)
# # 1.06565
# library(ggpubr)
# library(reshape2)
# library(data.table)
# mydf <- spikeinexprnorm
# mydf2 <- spikeinexprnorm * dd
# df1 <- as.data.table(mydf)
# # df <- melt(mydf2[,1:2])
# p1 <- ggscatter(df1, x = "CXXC1.ctl", y = "CXXC1.IAA", 
# 	color = "black", shape = 21, size = 3, 
# 	add = "reg.line",
# 	add.params = list(color = "blue", fill = "lightgray"), 
# 	conf.int = TRUE, 
# 	cor.coef = TRUE, 
# 	cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
# 	)
# mydf2 <- spikeinexprnorm * dd
# df2 <- as.data.table(mydf2)
# # df <- melt(mydf2[,1:2])
# p2 <- ggscatter(df2, x = "CXXC1.ctl", y = "CXXC1.IAA", 
# 	color = "black", shape = 21, size = 3, 
# 	add = "reg.line",
# 	add.params = list(color = "blue", fill = "lightgray"), 
# 	conf.int = TRUE, 
# 	cor.coef = TRUE, 
# 	cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
# 	)
# df3 <- as.data.table(spikeinexpr)
# # df <- melt(mydf2[,1:2])
# p3 <- ggscatter(df3, x = "CXXC1.ctl", y = "CXXC1.IAA", 
# 	color = "black", shape = 21, size = 3, 
# 	add = "reg.line",
# 	add.params = list(color = "blue", fill = "lightgray"), 
# 	conf.int = TRUE, 
# 	cor.coef = TRUE, 
# 	cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
# 	)
# mydf4 <- spikeinexprdeg * d
# df4 <- as.data.table(mydf4)
# # df <- melt(mydf2[,1:2])
# p4 <- ggscatter(df4, x = "CXXC1.ctl", y = "CXXC1.IAA", 
# 	color = "black", shape = 21, size = 3, 
# 	add = "reg.line",
# 	add.params = list(color = "blue", fill = "lightgray"), 
# 	conf.int = TRUE, 
# 	cor.coef = TRUE, 
# 	cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
# 	)

# # > attr(normCounts(set1),"scaled:scale")[1]/attr(normCounts(set1),"scaled:scale")[2]
# # CXXC1.ctl 
# #   1.15213 
# # > sizeFactors(dds)[1]/sizeFactors(dds)[2]
# # CXXC1.ctl 
# #  1.125008


# attr(,"scaled:scale")
# CXXC1.ctl CXXC1.IAA HCFC1.ctl HCFC1.IAA 
#  9508.251 10570.406 10557.768  9105.874 
# attr(,"scaled:center")
# CXXC1.ctl CXXC1.IAA HCFC1.ctl HCFC1.IAA 
#  584.3263  567.5855  575.8818  576.3443 

# attr(normCounts(set1),"scaled:scale")
# CXXC1.ctl CXXC1.IAA HCFC1.ctl HCFC1.IAA 
# 1.0975845 0.9526570 1.0376812 0.9120773 
# sizeFactors(dds)
# CXXC1.ctl CXXC1.IAA HCFC1.ctl HCFC1.IAA 
# 1.0922259 0.9708608 1.0526705 0.9246650 

# ############################################################
# ## from BioConductor
# ## ----data, warning=FALSE----------------------------------
# library(RUVSeq)
# library(zebrafishRNASeq)
# data(zfGenes)
# head(zfGenes)
# tail(zfGenes)

# ## ----filter-----------------------------------------------
# filter <- apply(zfGenes, 1, function(x) length(x[x>5]) >= 2)
# filtered <- zfGenes[filter,]
# genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
# spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

# ## ----store_data-------------------------------------------
# x <- as.factor(rep(c("Ctl", "Trt"), each = 3))
# set <- newSeqExpressionSet(as.matrix(filtered),
#                            phenoData = data.frame(x, row.names = colnames(filtered)))
# set

# ## ----rle, fig.cap="No normalization.",fig.subcap=c("RLE plot","PCA plot")----
# library(RColorBrewer)
# colors <- brewer.pal(3, "Set2")
# plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = colors[x])
# plotPCA(set, col = colors[x], cex = 1.2)

# ## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
# set <- betweenLaneNormalization(set, which = "upper")
# plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = colors[x])
# plotPCA(set, col = colors[x], cex = 1.2)

# ## ----ruv_spikes, fig.cap="RUVg normalization based on spike-in controls.", fig.subcap=c("RLE plot","PCA plot")----
# set1 <- RUVg(set, spikes, k = 1)
# pData(set1)
# plotRLE(set1, outline = FALSE, ylim = c(-4, 4), col = colors[x])
# plotPCA(set1, col = colors[x], cex = 1.2)
# attr(normCounts(set1),"scaled:scale")
# filtered[rownames(filtered) %in% spikes, ]
