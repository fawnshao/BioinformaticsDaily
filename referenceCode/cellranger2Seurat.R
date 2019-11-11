library(Seurat)
library(ggpubr)
library(dplyr)
setwd("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/NS27.NUP93.scRNA")
ut.data <- Read10X(data.dir = "X1/outs/filtered_gene_bc_matrices/hg19/")
ut <- CreateSeuratObject(counts = ut.data, project = "UT", min.cells = 3, min.features = 200)

iaa.data <- Read10X(data.dir = "X2/outs/filtered_gene_bc_matrices/hg19/")
iaa <- CreateSeuratObject(counts = iaa.data, project = "IAA", min.cells = 3, min.features = 200)

nup93 <- merge(x = ut, y = iaa, add.cell.ids = c("ut", "iaa"), project = "NUP93")
nup93[["percent.mt"]] <- PercentageFeatureSet(nup93, pattern = "^MT-")
p1 <- VlnPlot(nup93, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "Seurat.VlnPlot.pdf", p1, width = 9, height = 6)

plot1 <- FeatureScatter(nup93, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(nup93, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 <- CombinePlots(plots = list(plot1, plot2))
ggsave(filename = "Seurat.CombinePlots.pdf", p2, width = 9, height = 6)

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts
# An object of class Seurat 
# 18252 features across 1280 samples within 1 assay 
# Active assay: RNA (18252 features)
nup93.sub <- subset(nup93, subset = nFeature_RNA > 100 & percent.mt < 5)
# An object of class Seurat 
# 18252 features across 99 samples within 1 assay 
# Active assay: RNA (18252 features)

nup93 <- NormalizeData(nup93, normalization.method = "LogNormalize", scale.factor = 10000)
nup93 <- FindVariableFeatures(nup93, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
tops <- head(VariableFeatures(nup93), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nup93)
plot2 <- LabelPoints(plot = plot1, points = tops, repel = TRUE)
p3 <- CombinePlots(plots = list(plot1, plot2))
ggsave(filename = "Seurat.CombinePlots.VariableFeaturePlot.pdf", p3, width = 15, height = 6)


all.genes <- rownames(nup93)
nup93 <- ScaleData(nup93, features = all.genes)
nup93 <- RunPCA(nup93, features = VariableFeatures(object = nup93))
# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction,  DimPlot, and DimHeatmap
# Examine and visualize PCA results a few different ways
print(nup93[["pca"]], dims = 1:5, nfeatures = 5)
p4 <- VizDimLoadings(nup93, dims = 1:2, reduction = "pca")
ggsave(filename = "Seurat.VizDimLoadings.pca.pdf", p4, width = 10, height = 6)
# p4 <- VizDimLoadings(nup93, dims = 1:2, reduction = "tsne")
# ggsave(filename = "Seurat.VizDimLoadings.tsne.pdf", p4, width = 10, height = 6)
# p4 <- VizDimLoadings(nup93, dims = 1:2, reduction = "umap")
# ggsave(filename = "Seurat.VizDimLoadings.umap.pdf", p4, width = 10, height = 6)

p6 <- DimPlot(nup93, reduction = "pca", group.by = "orig.ident")
ggsave(filename = "Seurat.DimPlot.pca.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(nup93, reduction = "pca", group.by = "seurat_clusters")
ggsave(filename = "Seurat.DimPlot.pca.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(nup93, reduction = "pca", shape.by = "seurat_clusters", group.by = "orig.ident")
ggsave(filename = "Seurat.DimPlot.pca.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)


DimHeatmap(nup93, dims = 1:15, cells = 500, balanced = TRUE)
nup93 <- FindNeighbors(nup93, dims = 1:10)
nup93 <- FindClusters(nup93, resolution = 0.5)
head(Idents(nup93), 15)

nup93 <- RunUMAP(nup93, dims = 1:10)
p6 <- DimPlot(nup93, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "Seurat.DimPlot.umap.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(nup93, reduction = "umap", group.by = "seurat_clusters")
ggsave(filename = "Seurat.DimPlot.umap.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(nup93, reduction = "umap", shape.by = "seurat_clusters", group.by = "orig.ident")
ggsave(filename = "Seurat.DimPlot.umap.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)

nup93 <- RunTSNE(nup93, dims = 1:10)
p6 <- DimPlot(nup93, reduction = "tsne", group.by = "orig.ident")
ggsave(filename = "Seurat.DimPlot.tsne.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(nup93, reduction = "tsne", group.by = "seurat_clusters")
ggsave(filename = "Seurat.DimPlot.tsne.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(nup93, reduction = "tsne", shape.by = "seurat_clusters", group.by = "orig.ident")
ggsave(filename = "Seurat.DimPlot.tsne.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)


# cluster0.markers <- FindMarkers(nup93, ident.1 = 0, min.pct = 0.25)
# head(cluster0.markers, n = 5)
# cluster1.markers <- FindMarkers(nup93, ident.1 = 1, min.pct = 0.25)
# head(cluster1.markers, n = 5)
# cluster0.markers[rownames(cluster0.markers)=="KITLG",]
# nup93.markers <- FindAllMarkers(nup93, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# nup93.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# cluster0.markers <- FindMarkers(nup93, ident.1 = 0, 
#                                 logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# cluster0.markers[rownames(cluster0.markers)=="KITLG",]
# # VlnPlot(nup93, features = c("IER2", "VMP1", "IER3", "KLF10", "NEAT1", "KITLG"))
# # p7 <- VlnPlot(nup93, features = c("IER2", "VMP1", "IER3", "KLF10", "NEAT1", "KITLG"))
# cluster4.markers <- FindMarkers(nup93, ident.1 = 0, ident.2 = 4, 
#                                 logfc.threshold = 0.5, test.use = "roc", only.pos = FALSE)
# cluster4.markers[rownames(cluster4.markers)=="KITLG",]
# p7 <- VlnPlot(nup93, features = rownames(cluster4.markers)[1:20], ncol = 5, pt.size = 0.1)
# ggsave(filename = "Seurat.VlnPlot.by.cluster4vscluster0.pdf", p7, width = 20, height = 12)
# p7 <- VlnPlot(nup93, features = rownames(cluster4.markers)[1:20], ncol = 5, pt.size = 0.1, group.by = "orig.ident")
# ggsave(filename = "Seurat.VlnPlot.by.cluster4vscluster0.groupbyUTvsIAA.pdf", p7, width = 20, height = 12)
# p7 <- VlnPlot(nup93, features = rownames(cluster4.markers)[1:20], ncol = 5, pt.size = 0.1, 
#               group.by = "seurat_clusters", idents = c(0,4))
# ggsave(filename = "Seurat.VlnPlot.by.cluster4vscluster0.only04.pdf", p7, width = 20, height = 12)
# p7 <- VlnPlot(nup93, features = rownames(cluster4.markers)[1:20], ncol = 5, pt.size = 0.1, 
#               group.by = "orig.ident", idents = c(0,4))
# ggsave(filename = "Seurat.VlnPlot.by.cluster4vscluster0.only04.groupbyUTvsIAA.pdf", p7, width = 20, height = 12)


ut.markers <- FindMarkers(nup93, ident.1 = "UT", ident.2 = "IAA", group.by = "orig.ident", 
                          logfc.threshold = 0.25, test.use = "roc", only.pos = FALSE)
ut.markers[rownames(ut.markers)=="KITLG",]
which(rownames(ut.markers)=="KITLG", arr.ind = TRUE)
p8 <- VlnPlot(nup93, features = rownames(ut.markers)[1:20], ncol = 5, pt.size = 0.1, group.by = "orig.ident")
ggsave(filename = "Seurat.VlnPlot.by.UTvsIAA.pdf", p8, width = 20, height = 12)
p8 <- VlnPlot(nup93, features = rownames(ut.markers)[1:20], ncol = 5, pt.size = 0.1)
ggsave(filename = "Seurat.VlnPlot.by.UTvsIAA.plotdefault.pdf", p8, width = 20, height = 12)

# cluster4.markers.def <- FindMarkers(nup93, ident.1 = 4, ident.2 = 0, min.pct = 0.25)

cluster0.vs.cluste4 <- FindMarkers(nup93, ident.1 = 0, ident.2 = 4, 
                                   logfc.threshold = 0.5, test.use = "roc", only.pos = FALSE)
cluster0.vs.cluste4[rownames(cluster0.vs.cluste4)=="KITLG",]
which(rownames(cluster0.vs.cluste4)=="KITLG", arr.ind = TRUE)

cluster0.vs.all <- FindMarkers(nup93, ident.1 = 0, 
                               logfc.threshold = 0.5, test.use = "roc", only.pos = FALSE)
cluster0.vs.all[rownames(cluster0.vs.all)=="KITLG",]
which(rownames(cluster0.vs.all)=="KITLG", arr.ind = TRUE)
p9 <- FeaturePlot(nup93, features = c(rownames(cluster0.vs.all)[1:5], "KITLG"))
p9
nup93.markers <- FindAllMarkers(nup93, only.pos = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
nup93.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
nup93.markers[rownames(nup93.markers)=="KITLG",]
topN <- nup93.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(nup93, features = topN$gene) + NoLegend()

p10 <- DoHeatmap(nup93, features = rownames(cluster0.vs.all)[1:50])
ggsave(filename = "Seurat.DoHeatmap.by.cluster0.vs.all.pdf", p10, width = 12, height = 6)
p10 <- DoHeatmap(nup93, features = rownames(cluster0.vs.all)[1:50], group.by = "orig.ident")
ggsave(filename = "Seurat.DoHeatmap.by.cluster0.vs.all.groupbyorigident.pdf", p10, width = 12, height = 6)

p10 <- DoHeatmap(nup93, features = rownames(cluster0.vs.cluste4)[1:50])
ggsave(filename = "Seurat.DoHeatmap.by.cluster0.vs.cluster4.pdf", p10, width = 12, height = 6)
p10 <- DoHeatmap(nup93, features = rownames(cluster0.vs.cluste4)[1:50], group.by = "orig.ident")
ggsave(filename = "Seurat.DoHeatmap.by.cluster0.vs.cluster4.groupbyorigident.pdf", p10, width = 12, height = 6)

p10 <- DoHeatmap(nup93, features = rownames(ut.markers)[1:50])
ggsave(filename = "Seurat.DoHeatmap.by.ut.vs.iaa.pdf", p10, width = 12, height = 6)
p10 <- DoHeatmap(nup93, features = rownames(ut.markers)[1:50], group.by = "orig.ident")
ggsave(filename = "Seurat.DoHeatmap.by.ut.vs.iaa.groupbyorigident.pdf", p10, width = 12, height = 6)


p11 <- VlnPlot(nup93, features = rownames(cluster0.vs.cluste4)[1:20], ncol = 5, pt.size = 0.1)
ggsave(filename = "Seurat.VlnPlot.by.cluster0.vs.cluster4.pdf", p11, width = 20, height = 12)
p11 <- VlnPlot(nup93, features = rownames(cluster0.vs.cluste4)[1:20], ncol = 5, pt.size = 0.1, group.by = "orig.ident")
ggsave(filename = "Seurat.VlnPlot.by.cluster0.vs.cluster4.groupbyUTvsIAA.pdf", p11, width = 20, height = 12)
p11 <- VlnPlot(nup93, features = rownames(cluster0.vs.cluste4)[1:20], ncol = 5, pt.size = 0.1, 
              group.by = "seurat_clusters", idents = c(0,4))
ggsave(filename = "Seurat.VlnPlot.by.cluster0.vs.cluster4.only04.pdf", p11, width = 20, height = 12)
p11 <- VlnPlot(nup93, features = rownames(cluster0.vs.cluste4)[1:20], ncol = 5, pt.size = 0.1, 
              group.by = "orig.ident", idents = c(0,4))
ggsave(filename = "Seurat.VlnPlot.by.cluster0.vs.cluster4.only04.groupbyUTvsIAA.pdf", p11, width = 20, height = 12)

# saveRDS(nup93, file = "Seurat.nup93.rds")
