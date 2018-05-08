# .libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(data.table)
library(ggplot2)
library(pheatmap)
# args <- commandArgs(TRUE)
args <- c("hkg.tsg.srtbyPCA.all")
input <- fread(args[1], sep = "\t", header = T)
### for old data.table
# scores <- data.matrix(input[,.SD, .SDcols = c(3:170)])
# genes <- as.matrix(input[, .SD, .SDcols = 1])
# class <- as.matrix(input[, .SD, .SDcols = 2])
### for new data.table
scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
rownames(scores) <- genes
types.sim <- as.vector(types)
types.sim[1:2551] <- "HKG"
types.sim[2552:5286] <- "singleTSG"

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
seqfeature <- scores[,1:12]
seqfeature01 <- apply(seqfeature, 2, range01)
homermotif <- scores[,13:375]
gtrdpeaks <- scores[,376:2562]
meth <- scores[,2563:2599]
histone <- scores[,2600:3485]
dnase <- scores[,3486:3524]
hic <- scores[,3525:3539]
pchic <- scores[,3540:3573]

save.image(paste(args[1], "detail.data.RData", sep = "."))

###################################################
#####         look at the following           #####
###################################################
## add hkg classification by PCA
############## FeatureSelection
args <- c("hkg.tsg.srtbyPCA.all")
load(paste(args[1], "detail.data.RData", sep = "."))
RunFS <- function(x, y = as.numeric(types), pre = args[1]){
	require(FeatureSelection)
	trainx <- x
	trainy <- y
	pre <- pre
	params_glmnet <- list(alpha = 1, family = 'gaussian', nfolds = 5, parallel = TRUE)
	params_xgboost <- list(params = list("objective" = "reg:linear", "bst:eta" = 0.001, "subsample" = 0.75, "max_depth" = 5, "colsample_bytree" = 0.75, "nthread" = 60), nrounds = 1000, print.every.n = 250, maximize = FALSE)
	params_ranger <- list(dependent.variable.name = 'y', probability = FALSE, num.trees = 1000, verbose = TRUE, mtry = 5, min.node.size = 10, num.threads = 60, classification = FALSE, importance = 'permutation')
	params_features <- list(keep_number_feat = NULL, union = TRUE)
	feat <- wrapper_feat_select(X = trainx, y = trainy, params_glmnet = params_glmnet, params_xgboost = params_xgboost, params_ranger = params_ranger, xgb_sort = 'Gain', CV_folds = 5, stratified_regr = FALSE, scale_coefs_glmnet = FALSE, cores_glmnet = 5, params_features = params_features, verbose = TRUE)
	# save.image(paste(pre, "RData", sep = ".")
	write.table(feat$all_feat$"glmnet-lasso", paste(pre, "FeatureSelection.glmnetlasso.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$xgboost, paste(pre, "FeatureSelection.xgboost.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$ranger, paste(pre, "FeatureSelection.ranger.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$union_feat, paste(pre, "FeatureSelection.union_feat.tsv", sep = "."), row.names = F, sep = "\t")
	return(feat)
}
seqfeature.feat <- RunFS(x = seqfeature, y = as.numeric(types), pre = paste(args[1],"seqfeature", sep = "."))
seqfeature01.feat <- RunFS(x = seqfeature01, y = as.numeric(types), pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.feat <- RunFS(x = homermotif, y = as.numeric(types), pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.feat <- RunFS(x = gtrdpeaks, y = as.numeric(types), pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.feat <- RunFS(x = meth, y = as.numeric(types), pre = paste(args[1],"meth", sep = "."))
histone.feat <- RunFS(x = histone, y = as.numeric(types), pre = paste(args[1],"histone", sep = "."))
dnase.feat <- RunFS(x = dnase, y = as.numeric(types), pre = paste(args[1],"dnase", sep = "."))
hic.feat <- RunFS(x = hic, y = as.numeric(types), pre = paste(args[1],"hic", sep = "."))
pchic.feat <- RunFS(x = pchic, y = as.numeric(types), pre = paste(args[1],"pchic", sep = "."))
save.image(paste(args[1], "FeatureSelection.RData", sep = "."))
# cat onlyhkg.tsg.further.srt.all.*.FeatureSelection.union_feat.tsv | awk '$2>0.5'
# cat onlyhkg.tsg.further.srt.all.dnase.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.FeatureSelection.union_feat.tsv | awk '$2>0.5' > important.features.tsv
############## FeatureSelection use simplyfy classification
args <- c("hkg.tsg.srtbyPCA.all")
load(paste(args[1], "detail.data.RData", sep = "."))
RunFS <- function(x, y = as.numeric(factor(types.sim)), pre = args[1]){
	require(FeatureSelection)
	trainx <- x
	trainy <- y
	pre <- pre
	params_glmnet <- list(alpha = 1, family = 'gaussian', nfolds = 5, parallel = TRUE)
	params_xgboost <- list( params = list("objective" = "reg:linear", "bst:eta" = 0.001, "subsample" = 0.75, "max_depth" = 5, "colsample_bytree" = 0.75, "nthread" = 60), nrounds = 1000, print.every.n = 250, maximize = FALSE)
	params_ranger <- list(dependent.variable.name = 'y', probability = FALSE, num.trees = 1000, verbose = TRUE, mtry = 5, min.node.size = 10, num.threads = 60, classification = FALSE, importance = 'permutation')
	params_features <- list(keep_number_feat = NULL, union = TRUE)
	feat <- wrapper_feat_select(X = trainx, y = trainy, params_glmnet = params_glmnet, params_xgboost = params_xgboost, params_ranger = params_ranger, xgb_sort = 'Gain', CV_folds = 5, stratified_regr = FALSE, scale_coefs_glmnet = FALSE, cores_glmnet = 5, params_features = params_features, verbose = TRUE)
	# save.image(paste(pre, "RData", sep = ".")
	write.table(feat$all_feat$"glmnet-lasso", paste(pre, "FeatureSelection.glmnetlasso.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$xgboost, paste(pre, "FeatureSelection.xgboost.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$all_feat$ranger, paste(pre, "FeatureSelection.ranger.tsv", sep = "."), row.names = F, sep = "\t")
	write.table(feat$union_feat, paste(pre, "FeatureSelection.union_feat.tsv", sep = "."), row.names = F, sep = "\t")
	return(feat)
}
seqfeature.feat <- RunFS(x = seqfeature, pre = paste(args[1],"seqfeature.sim", sep = "."))
seqfeature01.feat <- RunFS(x = seqfeature01, pre = paste(args[1],"seqfeature01.sim", sep = "."))
homermotif.feat <- RunFS(x = homermotif, pre = paste(args[1],"homermotif.sim", sep = "."))
gtrdpeaks.feat <- RunFS(x = gtrdpeaks, pre = paste(args[1],"gtrdpeaks.sim", sep = "."))
meth.feat <- RunFS(x = meth, pre = paste(args[1],"meth.sim", sep = "."))
histone.feat <- RunFS(x = histone, pre = paste(args[1],"histone.sim", sep = "."))
dnase.feat <- RunFS(x = dnase, pre = paste(args[1],"dnase.sim", sep = "."))
hic.feat <- RunFS(x = hic, pre = paste(args[1],"hic.sim", sep = "."))
pchic.feat <- RunFS(x = pchic, pre = paste(args[1],"pchic.sim", sep = "."))
save.image(paste(args[1], "FeatureSelection.sim.RData", sep = "."))

############## PCA
args <- c("hkg.tsg.srtbyPCA.all")
load(paste(args[1], "detail.data.RData", sep = "."))
RunPCA <- function(x, y = types, z = types.sim, pre = args[1]){
	require(ggplot2)
	data <- x
	mycolors <- y
	mycolors.sim <- z
	pre <- pre
	myPCA.a <- prcomp(data, scale. = F, center = F)
	myPCA.a.t <- prcomp(t(data), scale. = F, center = F)
	pca.scores <- as.data.frame(myPCA.a$x)
	pdf(file = paste(pre, "pca.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
		geom_point(aes(colour = mycolors, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
	pdf(file = paste(pre, "pca.sim.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
		geom_point(aes(colour = mycolors.sim, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()

	pca.scores <- as.data.frame(myPCA.a.t$x)
	pdf(file = paste(pre, "t.pca.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) +
		geom_point(aes(alpha = 1/10)) +
		geom_text(colour = "tomato", alpha = 0.8, size = 2) +
		ggtitle(args[1])
	print(p1)
	dev.off()
}
seqfeature.pca <- RunPCA(x = seqfeature, pre = paste(args[1],"seqfeature", sep = "."))
seqfeature01.pca <- RunPCA(x = seqfeature01, pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.pca <- RunPCA(x = homermotif, pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.pca <- RunPCA(x = gtrdpeaks, pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.pca <- RunPCA(x = meth, pre = paste(args[1],"meth", sep = "."))
histone.pca <- RunPCA(x = histone, pre = paste(args[1],"histone", sep = "."))
dnase.pca <- RunPCA(x = dnase, pre = paste(args[1],"dnase", sep = "."))
hic.pca <- RunPCA(x = hic, pre = paste(args[1],"hic", sep = "."))
pchic.pca <- RunPCA(x = pchic, pre = paste(args[1],"pchic", sep = "."))
save.image(paste(args[1], "PCA.RData", sep = "."))

############## t-SNE
# conda.R
args <- c("hkg.tsg.srtbyPCA.all")
load(paste(args[1], "detail.data.RData", sep = "."))
.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(ggplot2)
library(Rtsne)
set.seed(123)
RuntSNE <- function(x, pre = args[1]){
	data <- x
	pre <- pre
	tsne <- Rtsne(data, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
	v.tsne <- as.data.frame(tsne$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	write.table(v.tsne, paste(pre, "tSNE.tsv", sep = "."), row.names = T, sep = "\t")
	return(tsne)
}
seqfeature.tsne <- RuntSNE(x = seqfeature, pre = paste(args[1],"seqfeature", sep = "."))
seqfeature01.tsne <- RuntSNE(x = seqfeature01, pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.tsne <- RuntSNE(x = homermotif, pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.tsne <- RuntSNE(x = gtrdpeaks, pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.tsne <- RuntSNE(x = meth, pre = paste(args[1],"meth", sep = "."))
histone.tsne <- RuntSNE(x = histone, pre = paste(args[1],"histone", sep = "."))
dnase.tsne <- RuntSNE(x = dnase, pre = paste(args[1],"dnase", sep = "."))
hic.tsne <- RuntSNE(x = hic, pre = paste(args[1],"hic", sep = "."))
pchic.tsne <- RuntSNE(x = pchic, pre = paste(args[1],"pchic", sep = "."))
save.image(paste(args[1], "detail.tSNE.RData", sep = "."))
#### plot t-SNE
library(ggplot2)
args <- c("hkg.tsg.srtbyPCA.all")
load(paste(args[1], "detail.tSNE.RData", sep = "."))
PlottSNE <- function(x, y = types, z = types.sim, pre = args[1]){
	require(ggplot2)
	data <- x
	mycolors <- y
	mycolors.sim <- z
	pre <- pre
	v.tsne <- as.data.frame(data$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	pdf(file = paste(pre, "tsne.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
		geom_point(aes(colour = mycolors, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
	pdf(file = paste(pre, "tsne.sim.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
		geom_point(aes(colour = mycolors.sim, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
}
seqfeature.tsne.plot <- PlottSNE(x = seqfeature.tsne, pre = paste(args[1],"seqfeature", sep = "."))
seqfeature01.tsne.plot <- PlottSNE(x = seqfeature01.tsne, pre = paste(args[1],"seqfeature01", sep = "."))
homermotif.tsne.plot <- PlottSNE(x = homermotif.tsne, pre = paste(args[1],"homermotif", sep = "."))
gtrdpeaks.tsne.plot <- PlottSNE(x = gtrdpeaks.tsne, pre = paste(args[1],"gtrdpeaks", sep = "."))
meth.tsne.plot <- PlottSNE(x = meth.tsne, pre = paste(args[1],"meth", sep = "."))
histone.tsne.plot <- PlottSNE(x = histone.tsne, pre = paste(args[1],"histone", sep = "."))
dnase.tsne.plot <- PlottSNE(x = dnase.tsne, pre = paste(args[1],"dnase", sep = "."))
hic.tsne.plot <- PlottSNE(x = hic.tsne, pre = paste(args[1],"hic", sep = "."))
pchic.tsne.plot <- PlottSNE(x = pchic.tsne, pre = paste(args[1],"pchic", sep = "."))

# please also ref to: 
# http://r-statistics.co/Variable-Selection-and-Importance-With-R.html
###################################################
#####           look at the above             #####
###################################################

############## visualization of the important features
# cat onlyhkg.tsg.further.srt.all.dnase.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.FeatureSelection.union_feat.tsv | awk '$2>0.8' > important.features.tsv
# cat onlyhkg.tsg.further.srt.all.dnase.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.pchic.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.hic.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.seqfeature01.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.homermotif.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.meth.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.gtrdpeaks.sim.FeatureSelection.union_feat.tsv onlyhkg.tsg.further.srt.all.histone.sim.FeatureSelection.union_feat.tsv | awk '$2>0.8' > sim.important.features.tsv
#############
# refine
library(data.table)
library(ggplot2)
library(pheatmap)
args <- c("hkg.tsg.srtbyPCA.all2")
# load(paste(args[1], "detail.data.RData", sep = "."))
input <- fread(args[1], sep = "\t", header = T)
### for new data.table
scores <- data.matrix(input[,-c(1:2)])
genes <- as.matrix(input[,1])
class <- as.matrix(input[,2])
types <- factor(class)
rownames(scores) <- genes
types.sim <- as.vector(types)
types.sim[1:2551] <- "HKG"
types.sim[2552:5286] <- "singleTSG"
colnames(scores)[3574:3644] <- gsub("-","_",colnames(scores)[3574:3644])

scores.x <- scores
myscale <- function(x){
	th <- as.vector(quantile(x, probs = 0.95))
	x[x > th] <- th
	return(x)
}
scores <- apply(scores, 2, myscale)

bigdiffs <- read.table("bigdiff.features.tsv", sep = "\t", row.names = 1)
scores <- scores[,colnames(scores) %in% rownames(bigdiffs)]

feats <- read.table("important.features.tsv", sep = "\t", row.names = 1)
data <- scores[,colnames(scores) %in% rownames(feats)]
# data.1 <- scores[,match(rownames(feats),colnames(scores))]
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.x <- apply(data, 2, range01)
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "important.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "important.pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, gaps_row = c(2551, 5286),
	cluster_cols = F, cluster_rows = F)
dev.off()
#########
feats <- read.table("sim.important.features.tsv", sep = "\t", row.names = 1)
# library(pheatmap)
data <- scores[,colnames(scores) %in% rownames(feats)]
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.x <- apply(data, 2, range01)
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "sim.important.pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, gaps_row = c(2552, 5287),
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "sim.important.pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, gaps_row = c(2551, 5286),
	cluster_cols = F, cluster_rows = F)
dev.off()

#########
# manually selected
feats <- read.table("important.features.tsv", sep = "\t", row.names = 1)
data <- scores[,colnames(scores) %in% rownames(feats)]
data <- data.frame(scores[,1], data)
colnames(data)[1] <- colnames(scores)[1]
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.xx <- apply(data, 2, range01)
# data.x <- data.x[,-(10:78)]
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.xx)
# colors <- colorRampPalette(c("blue", "white", "red"))(100)
colors <- colorRampPalette(c("white","blue"))(100)
data.x <- data.xx[,c(1:9,19,20,22,33,38,45,85:88)]
png(filename = paste(args[1], "important.man.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = F, cluster_rows = F)
dev.off()


annos.sim <- data.frame(class = types.sim)
rownames(annos.sim) <- rownames(data.x)
png(filename = paste(args[1], "important.man.c.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = T)
dev.off()
png(filename = paste(args[1], "important.man.c1.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = F)
dev.off()

save.image(paste(args[1], "selected.RData", sep = "."))
####conda.newR
args <- c("hkg.tsg.srtbyPCA.all2")
load(paste(args[1], "selected.RData", sep = "."))
.libPaths("/work/04935/shaojf/stampede2/myTools/anaconda2/lib/R/library")
library(ggplot2)
library(Rtsne)
library(pheatmap)
set.seed(123)
RuntSNE <- function(x, pre = args[1]){
	data <- x
	pre <- pre
	tsne <- Rtsne(data, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
	v.tsne <- as.data.frame(tsne$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	write.table(v.tsne, paste(pre, "tSNE.tsv", sep = "."), row.names = T, sep = "\t")
	return(tsne)
}
PlottSNE <- function(x, y = types, z = types.sim, pre = args[1]){
	require(ggplot2)
	data <- x
	mycolors <- y
	mycolors.sim <- z
	pre <- pre
	v.tsne <- as.data.frame(data$Y)
	colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
	pdf(file = paste(pre, "tsne.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
		geom_point(aes(colour = mycolors, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
	pdf(file = paste(pre, "tsne.sim.pdf", sep = "."), width = 15, height = 12)
	p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
		geom_point(aes(colour = mycolors.sim, alpha = 1/10)) +
		ggtitle(pre)
	print(p1)
	dev.off()
}
myscale <- function(x){
	th <- as.vector(floor(quantile(x, probs = 0.95)))
	x[x > th] <- th
	return(x)
}
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
######################
IEratio <- scores.x[,4]/scores.x[,3]
summary(IEratio[1:2551])
summary(IEratio[2552:5286])
summary(IEratio[5287:6860])
# IEratio.scale <- myscale(IEratio)
IEratio.scale <- IEratio
IEratio.scale[IEratio.scale > 20] <- 20
IEratio.scale <- range01(IEratio)
HiChIP.Th17 <- range01(myscale(scores.x[,3538]))
HiChIP.Treg <- range01(myscale(scores.x[,3539]))
mydata <- data.frame(IEratio.scale, data.x)
mydata.tsne <- RuntSNE(x = mydata, pre = paste(args[1],"test", sep = "."))
mydata.plot <- PlottSNE(x = mydata.tsne, pre = paste(args[1],"test", sep = "."))
mydata.1 <- data.frame(IEratio.scale, data.x, HiChIP.Th17, HiChIP.Treg)
mydata.1.tsne <- RuntSNE(x = mydata.1, pre = paste(args[1],"test.1", sep = "."))
mydata.1.plot <- PlottSNE(x = mydata.1.tsne, pre = paste(args[1],"test.1", sep = "."))

colors <- colorRampPalette(c("white","blue"))(100)
png(filename = paste(args[1], "test.1.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(mydata.1, scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = F)
dev.off()
png(filename = paste(args[1], "test.1.1.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(mydata.1, scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = T)
dev.off()
######################
myselection.tsne <- RuntSNE(x = data.xx, pre = paste(args[1],"myselection", sep = "."))
myselection.plot <- PlottSNE(x = myselection.tsne, pre = paste(args[1],"myselection", sep = "."))

myselection.tsne.1 <- RuntSNE(x = data.xx[,1:88], pre = paste(args[1],"myselection.noenhancer", sep = "."))
myselection.plot.1 <- PlottSNE(x = myselection.tsne.1, pre = paste(args[1],"myselection.noenhancer", sep = "."))

myselection.tsne.2 <- RuntSNE(x = data.x, pre = paste(args[1],"myselection.noenhancer.noDHS", sep = "."))
myselection.plot.2 <- PlottSNE(x = myselection.tsne.2, pre = paste(args[1],"myselection.noenhancer.noDHS", sep = "."))

myselection.tsne.3 <- RuntSNE(x = data.x[,c(1,10:19)], pre = paste(args[1],"myselection.noenhancer.noDHS.1", sep = "."))
myselection.plot.3 <- PlottSNE(x = myselection.tsne.3, pre = paste(args[1],"myselection.noenhancer.noDHS.1", sep = "."))


png(filename = paste(args[1], "important.man.c2.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data.x[,c(1,10:19)], scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = F)
dev.off()
png(filename = paste(args[1], "important.man.c22.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(data.x[,c(1,10:19)], scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = T)
dev.off()

myscale <- function(x){
	th <- as.vector(quantile(x, probs = 0.95))
	x[x > th] <- th
	return(x)
}
seqf.loop <- apply(scores.x[,c(1:12,3535:3539,3547:3549)], 2, myscale)
seqf.loop <- apply(seqf.loop, 2, range01)

png(filename = paste(args[1], "seqf.loop.pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(seqf.loop, scale = "none", annotation_row = annos.sim,
	show_rownames = F, show_colnames = T, color = colors, 
	cluster_cols = T, cluster_rows = F)
dev.off()
myselection.tsne.4 <- RuntSNE(x = seqf.loop, pre = paste(args[1],"seqf.loop", sep = "."))
myselection.plot.4 <- PlottSNE(x = myselection.tsne.4, pre = paste(args[1],"seqf.loop", sep = "."))

IEratio <- scores.x[,4]/scores.x[,3]
summary(IEratio[1:2551])
summary(IEratio[2552:5286])
summary(IEratio[5287:6860])
IEratio.scale <- myscale(IEratio)
IEratio.scale <- range01(IEratio)
seqf.loop.1 <- data.frame(IEratio.scale,seqf.loop[,c(2,3,7:10,12,13:20)])
myselection.tsne.5 <- RuntSNE(x = seqf.loop.1, pre = paste(args[1],"seqf.loop.1", sep = "."))
myselection.plot.5 <- PlottSNE(x = myselection.tsne.5, pre = paste(args[1],"seqf.loop.1", sep = "."))


seqf.motif.loop <- apply(scores.x[,c(2,3,7:10,12,84,94,219,302,345,3535:3539,3547:3549)], 2, myscale)
seqf.motif.loop <- apply(seqf.motif.loop, 2, range01)
seqf.motif.loop.1 <- data.frame(IEratio.scale, seqf.motif.loop)
myselection.tsne.6 <- RuntSNE(x = seqf.motif.loop.1, pre = paste(args[1],"seqf.motif.loop.1", sep = "."))
myselection.plot.6 <- PlottSNE(x = myselection.tsne.6, pre = paste(args[1],"seqf.motif.loop.1", sep = "."))

myselection.tsne.7 <- RuntSNE(x = seqf.motif.loop.1[,c(1,7:21)], pre = paste(args[1],"seqf.motif.loop.2", sep = "."))
myselection.plot.7 <- PlottSNE(x = myselection.tsne.7, pre = paste(args[1],"seqf.motif.loop.2", sep = "."))

save.image(paste(args[1], "tmp.RData", sep = "."))
#########
# all
# scores.x <- scores
# myscale <- function(x){
# 	th <- as.vector(floor(quantile(x, probs = 0.95)))
# 	x[x > th] <- th
# 	return(x)
# }
data <- apply(scores.x, 2, myscale)
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
data.x <- apply(data, 2, range01)
annos <- data.frame(class = types)
rownames(annos) <- rownames(data.x)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 2500, height = 1200)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, gaps_row = c(2551, 5286),
	cluster_cols = F, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "pheatmap.pdf", sep = "."), width = 25, height = 12)
myplot <- pheatmap(data.x, scale = "none", annotation_row = annos,
	show_rownames = F, show_colnames = T, color = colors, gaps_row = c(2551, 5286),
	cluster_cols = F, cluster_rows = F)
dev.off()
#########


