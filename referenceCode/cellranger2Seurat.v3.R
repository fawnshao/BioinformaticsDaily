library(Seurat)
library(ggpubr)
library(dplyr)
library(data.table)
library(gridExtra)
library(cowplot)
# setwd("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/NS27.NUP93.scRNA")
setwd("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/NS27.NUP93.scRNA/gene.eRNA/seurat/")
nup93.ut.data <- Read10X(data.dir = "X1/outs/filtered_gene_bc_matrices/hg19/")
nup93.ut <- CreateSeuratObject(counts = nup93.ut.data, project = "NUP93.UT", min.cells = 3, min.features = 200)
nup93.iaa.data <- Read10X(data.dir = "X2/outs/filtered_gene_bc_matrices/hg19/")
nup93.iaa <- CreateSeuratObject(counts = nup93.iaa.data, project = "NUP93.IAA", min.cells = 3, min.features = 200)

rad21.ut.data <- Read10X(data.dir = "R1/outs/filtered_gene_bc_matrices/hg19/")
rad21.ut <- CreateSeuratObject(counts = rad21.ut.data, project = "RAD21.UT", min.cells = 3, min.features = 200)
rad21.iaa.data <- Read10X(data.dir = "R2/outs/filtered_gene_bc_matrices/hg19/")
rad21.iaa <- CreateSeuratObject(counts = rad21.iaa.data, project = "RAD21.IAA", min.cells = 3, min.features = 200)

nup93 <- merge(x = nup93.ut, y = nup93.iaa, add.cell.ids = c("nup93.ut", "nup93.iaa"), project = "NUP93")
rad21 <- merge(x = rad21.ut, y = rad21.iaa, add.cell.ids = c("rad21.ut", "rad21.iaa"), project = "RAD21")
hct116.aid <- merge(x = nup93, y = rad21, project = "HCT116.AID")

# QC --------------------------------------------------------
hct116.aid[["percent.mt"]] <- PercentageFeatureSet(hct116.aid, pattern = "^MT-")
hct116.aid[["percent.eRNA"]] <- PercentageFeatureSet(hct116.aid, pattern = "^eRNA")
p1 <- VlnPlot(hct116.aid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.eRNA"), ncol = 4)
ggsave(filename = "hct116.aid.Seurat.VlnPlot.pdf", p1, width = 15, height = 6)
# hct116.aid$percent.mt[hct116.aid$percent.mt>50]
# head -1 R2/R2.matrix.csv | tr "," "\n" | cat -n | grep "GTCATGAAGAAGCTGC"
# 358	GTCATGAAGAAGCTGC-1
# awk -F"," '{print $1,$358}' R2/R2.matrix.csv 
# hct116.aid$nFeature_RNA[hct116.aid$nFeature_RNA < 500]
# head -1 R2/R2.matrix.csv | tr "," "\n" | cat -n | grep "GTCCCATAGCACGTCC"
# 363	GTCCCATAGCACGTCC-1
# awk -F"," '$363>0{print $1,$363}' R2/R2.matrix.csv | wc -l
# 204

basic.stats <- as.data.table(hct116.aid[[]])
setorder(basic.stats, nCount_RNA)
basic.stats[, nCount_RNA.rank := seq(1:nrow(basic.stats))]
setorder(basic.stats, nFeature_RNA)
basic.stats[, nFeature_RNA.rank := seq(1:nrow(basic.stats))]
setorder(basic.stats, percent.mt)
basic.stats[, percent.mt.rank := seq(1:nrow(basic.stats))]
setorder(basic.stats, percent.eRNA)
basic.stats[, percent.eRNA.rank := seq(1:nrow(basic.stats))]
p0 <- ggscatter(basic.stats, x = "nCount_RNA.rank", y = "nCount_RNA", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
p1 <- ggscatter(basic.stats, x = "nFeature_RNA.rank", y = "nFeature_RNA", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
p2 <- ggscatter(basic.stats, x = "percent.mt.rank", y = "percent.mt", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
p3 <- ggscatter(basic.stats, x = "percent.eRNA.rank", y = "percent.eRNA", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
ggsave(filename = "hct116.aid.basic.stats.pdf",
       grid.arrange(p0, p1, p2, p3, ncol = 1), 
       width = 10, height = 12)
# filter by MT genes --------------------------------------------------------
# hct116.aid.sub <- subset(hct116.aid, subset = nFeature_RNA > 3000 & nFeature_RNA < 10000 & nCount_RNA > 10000 & nCount_RNA < 200000 & percent.mt > 5 & percent.mt < 25)
# for gene & eRNA
# hct116.aid.sub <- subset(hct116.aid, subset = nFeature_RNA > 5000 & nFeature_RNA < 12000 & nCount_RNA > 10000 & nCount_RNA < 200000 & percent.mt > 5 & percent.mt < 25)
`%notin%` <- Negate(`%in%`)
hct116.aid.sub <- subset(hct116.aid, subset = nFeature_RNA %notin% boxplot.stats(hct116.aid$nFeature_RNA)$out & 
                             nCount_RNA %notin% boxplot.stats(hct116.aid$nCount_RNA)$out & 
                             percent.mt %notin% boxplot.stats(hct116.aid$percent.mt)$out)

plot1 <- FeatureScatter(hct116.aid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hct116.aid, feature1 = "nCount_RNA", feature2 = "percent.eRNA")
plot3 <- FeatureScatter(hct116.aid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 <- CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)
ggsave(filename = "hct116.aid.Seurat.CombinePlots.pdf", p2, width = 15, height = 6)
plot1 <- FeatureScatter(hct116.aid.sub, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hct116.aid.sub, feature1 = "nCount_RNA", feature2 = "percent.eRNA")
plot3 <- FeatureScatter(hct116.aid.sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 <- CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)
ggsave(filename = "hct116.aid.sub.Seurat.CombinePlots.pdf", p2, width = 15, height = 6)

basic.stats <- as.data.table(hct116.aid.sub[[]])
setorder(basic.stats, nCount_RNA)
basic.stats[, nCount_RNA.rank := seq(1:nrow(basic.stats))]
setorder(basic.stats, nFeature_RNA)
basic.stats[, nFeature_RNA.rank := seq(1:nrow(basic.stats))]
setorder(basic.stats, percent.mt)
basic.stats[, percent.mt.rank := seq(1:nrow(basic.stats))]
setorder(basic.stats, percent.eRNA)
basic.stats[, percent.eRNA.rank := seq(1:nrow(basic.stats))]
p0 <- ggscatter(basic.stats, x = "nCount_RNA.rank", y = "nCount_RNA", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
p1 <- ggscatter(basic.stats, x = "nFeature_RNA.rank", y = "nFeature_RNA", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
p2 <- ggscatter(basic.stats, x = "percent.mt.rank", y = "percent.mt", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
p3 <- ggscatter(basic.stats, x = "percent.eRNA.rank", y = "percent.eRNA", 
                color = "orig.ident", fill = "orig.ident", palette = "jco", size = 0.2)
ggsave(filename = "hct116.aid.sub.basic.stats.pdf",
       grid.arrange(p0, p1, p2, p3, ncol = 1), 
       width = 10, height = 12)

# before removing cell-cycle effects --------------------------------------------------------
hct116.aid.sub <- NormalizeData(hct116.aid.sub, normalization.method = "LogNormalize", scale.factor = 10000)
hct116.aid.sub <- FindVariableFeatures(hct116.aid.sub, selection.method = "vst", nfeatures = 2000)
# hct116.aid.sub <- FindVariableFeatures(hct116.aid.sub, selection.method = "vst", nfeatures = 1000)

all.genes <- rownames(hct116.aid.sub)
hct116.aid.sub <- ScaleData(hct116.aid.sub, features = all.genes)
hct116.aid.sub <- RunPCA(hct116.aid.sub, features = VariableFeatures(object = hct116.aid.sub))
p1 <- ElbowPlot(hct116.aid.sub, ndims = 50)
ggsave(filename = "hct116.aid.sub.PCA.ElbowPlot.pdf", p1, width = 10, height = 12)
hct116.aid.sub <- RunUMAP(hct116.aid.sub, dims = 1:10)
hct116.aid.sub <- RunTSNE(hct116.aid.sub, dims = 1:10)
# hct116.aid.sub <- RunUMAP(hct116.aid.sub, dims = 1:10, features = VariableFeatures(object = hct116.aid.sub))
# Error in RunUMAP.Seurat(hct116.aid.sub, dims = 1:10, features = VariableFeatures(object = hct116.aid.sub)) : 
#     Please specify only one of the following arguments: dims, features, or graph
# 

hct116.aid.sub <- FindNeighbors(hct116.aid.sub, reduction = "pca", dims = 1:10, k.param = 100)
hct116.aid.sub <- FindClusters(hct116.aid.sub, resolution = 1.5)
p1 <- DimPlot(hct116.aid.sub, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(hct116.aid.sub, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

plots <- VlnPlot(hct116.aid.sub, features = c("KITLG", "IER2", "VMP1"), 
                 split.by = "orig.ident", group.by = "seurat_clusters", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)


print(hct116.aid.sub[["pca"]], dims = 1:10, nfeatures = 5)
p4 <- VizDimLoadings(hct116.aid.sub, dims = 1:2, reduction = "pca")
ggsave(filename = "hct116.aid.sub.Seurat.VizDimLoadings.pca.pdf", p4, width = 10, height = 6)

p6 <- DimPlot(hct116.aid.sub, reduction = "pca", group.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.pca.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "pca", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.pca.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "pca", group.by = "seurat_clusters", shape.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.pca.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "pca", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.pca.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

p6 <- DimPlot(hct116.aid.sub, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.umap.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "umap", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.umap.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "umap", group.by = "seurat_clusters", shape.by= "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.umap.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.umap.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

p6 <- DimPlot(hct116.aid.sub, reduction = "tsne", group.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.tsne.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "tsne", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.tsne.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "tsne", group.by = "seurat_clusters", shape.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.tsne.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub, reduction = "tsne", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.tsne.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

aid.markers <- FindAllMarkers(hct116.aid.sub, only.pos = FALSE, test.use = "DESeq2",
                              min.pct = 0.1, logfc.threshold = 0.25)
topN <- aid.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
aid.markers[rownames(aid.markers)=="KITLG",]
aid.markers[rownames(aid.markers)=="eRNA-5690",]
aid.markers[rownames(aid.markers)=="eRNA-19175",]
aid.markers[rownames(aid.markers)=="eRNA-22479",]
# UTRN
aid.markers[grep("eRNA-", rownames(aid.markers)),]
aid.markers[rownames(aid.markers)=="CCAT1",]
# RidgePlot(hct116.aid.sub, features = c("KITLG", topN$gene), group.by = "orig.ident")
mygenes <- c("eRNA-19175","IER2", "KITLG", "DYNC1H1", "NCAPG", "PRKDC")
p1 <- RidgePlot(hct116.aid.sub, features = mygenes, ncol = 6)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.topgenes.bycluster.pdf", p1, width = 25, height = 8)
p1 <- RidgePlot(hct116.aid.sub, features = mygenes, 
                group.by = "orig.ident", ncol = 6)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.topgenes.pdf", p1, width = 25, height = 4)
p1 <- RidgePlot(hct116.aid.sub, features = mygenes, 
                idents = c(0:1,5:7), group.by = "orig.ident", ncol = 6)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.topgenes.selected.pdf", p1, width = 25, height = 4)
p1 <- FeaturePlot(hct116.aid.sub, features = mygenes, 
                  cols = c("white", "blue"), shape.by = "orig.ident", ncol = 3)
ggsave(filename = "hct116.aid.sub.Seurat.FeaturePlot.topgenes.pdf", p1, width = 20, height = 8)



# Cell cycle genes --------------------------------------------------------
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
hct116.aid.sub <- CellCycleScoring(hct116.aid.sub, 
                                   s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(hct116.aid.sub[[]])
p1 <- RidgePlot(hct116.aid.sub, features = c("PCNA", "TOP2A", "MCM6", "MKI67", "CTCF", "NCAPG"), 
                group.by = "Phase", ncol = 3)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.ccgene.Phase.pdf", p1, width = 15, height = 4)
p1 <- RidgePlot(hct116.aid.sub, features = c("PCNA", "TOP2A", "MCM6", "MKI67", "CTCF", "NCAPG"), 
                group.by = "seurat_clusters", ncol = 3)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.ccgene.cluster.pdf", p1, width = 15, height = 9)
p1 <- RidgePlot(hct116.aid.sub, features = c(cc.genes$s.genes[1:3], cc.genes$g2m.genes[1:3], mygenes),
                group.by = "Phase", ncol = 6)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.ccgenes.topgenes.Phase.pdf", p1, width = 25, height = 9)
p1 <- RidgePlot(hct116.aid.sub, features = c(cc.genes$s.genes[1:3], cc.genes$g2m.genes[1:3], mygenes),
                group.by = "seurat_clusters", ncol = 6)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.ccgenes.topgenes.seurat_clusters.pdf", p1, width = 25, height = 9)
p1 <- RidgePlot(hct116.aid.sub, features = c(cc.genes$s.genes[1:3], cc.genes$g2m.genes[1:3], mygenes),
                group.by = "orig.ident", ncol = 6)
ggsave(filename = "hct116.aid.sub.Seurat.RidgePlot.ccgenes.topgenes.orig.ident.pdf", p1, width = 25, height = 9)


hct116.aid.sub.cc <- ScaleData(hct116.aid.sub, vars.to.regress = c("S.Score", "G2M.Score"),
                               features = rownames(hct116.aid.sub))
hct116.aid.sub.cc <- RunPCA(hct116.aid.sub.cc, features = VariableFeatures(object = hct116.aid.sub.cc))
hct116.aid.sub.cc <- RunUMAP(hct116.aid.sub.cc, dims = 1:10)
hct116.aid.sub.cc <- RunTSNE(hct116.aid.sub.cc, dims = 1:10)
hct116.aid.sub.cc <- FindNeighbors(hct116.aid.sub.cc, reduction = "pca", dims = 1:10, k.param = 100)
hct116.aid.sub.cc <- FindClusters(hct116.aid.sub.cc, resolution = 1.5)


p3 <- DimPlot(hct116.aid.sub.cc, reduction = "umap", group.by = "orig.ident")
p4 <- DimPlot(hct116.aid.sub.cc, reduction = "umap", label = TRUE)
plot_grid(p1, p2, p3,p4)

plots1 <- VlnPlot(hct116.aid.sub, features = c("KITLG", "IER2", "VMP1"), 
                 split.by = "orig.ident", group.by = "seurat_clusters", 
                 pt.size = 0, combine = FALSE)
plots2 <- VlnPlot(hct116.aid.sub.cc, features = c("KITLG", "IER2", "VMP1"), 
                 split.by = "orig.ident", group.by = "seurat_clusters", 
                 pt.size = 0, combine = FALSE)
ggsave(filename = "CombinePlots.selectedgenes.pdf", 
       plot_grid(CombinePlots(plots = plots1, ncol = 1), CombinePlots(plots = plots2, ncol = 1)),
       width = 15, height = 6)


print(hct116.aid.sub.cc[["pca"]], dims = 1:5, nfeatures = 5)
p4 <- VizDimLoadings(hct116.aid.sub.cc, dims = 1:2, reduction = "pca")
ggsave(filename = "hct116.aid.sub.cc.Seurat.VizDimLoadings.pca.pdf", p4, width = 10, height = 6)

p6 <- DimPlot(hct116.aid.sub.cc, reduction = "pca", group.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.Seurat.DimPlot.pca.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "pca", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.pca.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "pca", group.by = "seurat_clusters", shape.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.pca.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "pca", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.pca.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

p6 <- DimPlot(hct116.aid.sub.cc, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.umap.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "umap", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.umap.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "umap", group.by = "seurat_clusters", shape.by= "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.umap.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.umap.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

p6 <- DimPlot(hct116.aid.sub.cc, reduction = "tsne", group.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.tsne.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "tsne", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.tsne.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "tsne", group.by = "seurat_clusters", shape.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.tsne.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.sub.cc, reduction = "tsne", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DimPlot.tsne.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

###
p1 <- RidgePlot(hct116.aid.sub.cc, features = c("PCNA", "TOP2A", "MCM6", "MKI67", "CTCF", "NCAPG"), 
                group.by = "Phase", ncol = 3)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.ccgene.Phase.pdf", p1, width = 15, height = 4)
p1 <- RidgePlot(hct116.aid.sub.cc, features = c("PCNA", "TOP2A", "MCM6", "MKI67", "CTCF", "NCAPG"), 
                group.by = "seurat_clusters", ncol = 3)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.ccgene.cluster.pdf", p1, width = 15, height = 9)
p1 <- RidgePlot(hct116.aid.sub.cc, features = c(cc.genes$s.genes[1:3], cc.genes$g2m.genes[1:3], mygenes),
                group.by = "Phase", ncol = 6)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.ccgenes.topgenes.Phase.pdf", p1, width = 25, height = 9)
p1 <- RidgePlot(hct116.aid.sub.cc, features = c(cc.genes$s.genes[1:3], cc.genes$g2m.genes[1:3], mygenes),
                group.by = "seurat_clusters", ncol = 6)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.ccgenes.topgenes.seurat_clusters.pdf", p1, width = 25, height = 9)
p1 <- RidgePlot(hct116.aid.sub.cc, features = c(cc.genes$s.genes[1:3], cc.genes$g2m.genes[1:3], mygenes),
                group.by = "orig.ident", ncol = 6)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.ccgenes.topgenes.orig.ident.pdf", p1, width = 25, height = 9)
###


cc.aid.markers <- FindAllMarkers(hct116.aid.sub.cc, only.pos = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
cc.topN <- cc.aid.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
cc.aid.markers[rownames(cc.aid.markers)=="KITLG",]
cc.aid.markers[grep("eRNA-", rownames(cc.aid.markers)),]
cc.aid.markers[rownames(cc.aid.markers)=="CCAT1",]
# mygenes <- c("eRNA-19175","IER2", "KITLG", "DYNC1H1", "NCAPG", "PRKDC")

# mygenes.cc <- c(unique(cc.aid.markers[grep("eRNA-", rownames(cc.aid.markers)),]$gene))
mygenes.cc <- unique(c(cc.aid.markers[cc.aid.markers$cluster=="3" & cc.aid.markers$avg_logFC > 0,]$gene[1:15],
                       cc.aid.markers[cc.aid.markers$cluster=="3" & cc.aid.markers$avg_logFC < 0,]$gene[1:15]))
p1 <- RidgePlot(hct116.aid.sub.cc, features = mygenes.cc, group.by = "orig.ident", ncol = 15)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.topgenes.pdf", p1, width = 40, height = 6)
p1 <- RidgePlot(hct116.aid.sub.cc, features = mygenes.cc, idents = c(1:5), group.by = "orig.ident", ncol = 15)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.topgenes.selected.pdf", p1, width = 40, height = 6)
p1 <- FeaturePlot(hct116.aid.sub.cc, features = mygenes.cc, cols = c("white", "blue"), shape.by = "orig.ident", ncol = 5)
ggsave(filename = "hct116.aid.sub.cc.Seurat.FeaturePlot.topgenes.pdf", p1, width = 30, height = 30)

p10 <- DoHeatmap(hct116.aid.sub.cc, features = rownames(cc.aid.markers))
ggsave(filename = "hct116.aid.sub.cc.Seurat.DoHeatmap.FindAllMarkers.pdf", 
       p10, width = 10, height = 16)
p10 <- DoHeatmap(hct116.aid.sub.cc, features = rownames(cc.aid.markers), group.by = "orig.ident")
ggsave(filename = "hct116.aid.sub.cc.Seurat.DoHeatmap.FindAllMarkers.groupbyorigident.pdf", 
       p10, width = 10, height = 16)
# test
RidgePlot(hct116.aid.sub.cc, features = "TCF7L2", idents = c(1:5), group.by = "orig.ident", ncol = 1)
RidgePlot(hct116.aid.sub.cc, features = c("TCF7L2","VTI1A","ACSL5"), group.by = "orig.ident", ncol = 3)
RidgePlot(hct116.aid.sub.cc, features = "IRS1", group.by = "orig.ident", ncol = 1)
RidgePlot(hct116.aid.sub.cc, features = c("COL4A4","RHBDD1","IRS1"), group.by = "orig.ident", ncol = 3)
RidgePlot(hct116.aid.sub.cc, features = c("eRNA-9800","eRNA-9801","eRNA-9802","eRNA-9803","eRNA-9804","eRNA-9805"), group.by = "orig.ident", ncol = 3)
rownames(hct116.aid.sub.cc)[grep("eRNA-980",rownames(hct116.aid.sub.cc))]


###
finalegenes <- rownames(hct116.aid.sub.cc)
finalexprmat <- hct116.aid.sub.cc[["RNA"]]@data
finalexprmat[rownames(finalexprmat)=="KITLG",]
final.markers.nup93 <- FindMarkers(hct116.aid.sub.cc, ident.1 = "1", ident.2 = "0",
                                   only.pos = FALSE, test.use = "DESeq2", 
                                   min.pct = 0.1, logfc.threshold = log2(1.5))
final.markers.nup93[rownames(final.markers.nup93)=="KITLG",]
final.markers.nup93[grep("eRNA-", rownames(final.markers.nup93)),]
final.markers.nup93[grep("eRNA-54", rownames(final.markers.nup93)),]
final.markers.nup93[rownames(final.markers.nup93)=="CCAT1",]
p10 <- DoHeatmap(hct116.aid.final, features = rownames(final.markers.nup93)[1:100], group.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DoHeatmap.onlyNUP93.groupbyorigident.pdf", 
       p10, width = 10, height = 16)

final.scup <- final.markers.nup93[final.markers.nup93$avg_logFC > 0,]
final.scdown <- final.markers.nup93[final.markers.nup93$avg_logFC < 0,]




# compared to PROseq --------------------------------------------------------
proup <- fread("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/all.HCT116.NUP93.PROseq/PRO.bothTreatments.bothReps/transcribed.gene/DE.gene/noTNF.gene.up.txt", header = F)
prodown <- fread("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/all.HCT116.NUP93.PROseq/PRO.bothTreatments.bothReps/transcribed.gene/DE.gene/noTNF.gene.down.txt", header = F)
scup <- sub.cc.aid.markers.nup93[sub.cc.aid.markers.nup93$avg_logFC > 0,]
scdown <- sub.cc.aid.markers.nup93[sub.cc.aid.markers.nup93$avg_logFC < 0,]
proup.genes <- proup$V1
prodown.genes <- prodown$V1
scup.genes <- rownames(scup)[grep("eRNA", rownames(scup), invert = T)]
scdown.genes <- rownames(scdown)[grep("eRNA", rownames(scdown), invert = T)]
library(gplots)
uplist <- venn(list(PROseq.UP = proup.genes, scRNAseq.UP = scup.genes))
downlist <- venn(list(PROseq.DOWN = prodown.genes, scRNAseq.DOWN = scdown.genes))

attr(downlist, "intersections")$PROseq.DOWN
attr(downlist, "intersections")$scRNAseq.DOWN
attr(uplist, "intersections")$PROseq.UP
attr(uplist, "intersections")$scRNAseq.UP
library(VennDiagram)
dev.off()
pdf(file = "venn.hct116.aid.sub.cc.Seurat.onlyNUP93.downgenes.pdf", width = 4, height = 4)
draw.pairwise.venn(area1 = length(attr(downlist, "intersections")$PROseq.DOWN) + length(attr(downlist, "intersections")$`PROseq.DOWN:scRNAseq.DOWN`), 
                   area2 = length(attr(downlist, "intersections")$scRNAseq.DOWN) + length(attr(downlist, "intersections")$`PROseq.DOWN:scRNAseq.DOWN`), 
                   cross.area = length(attr(downlist, "intersections")$`PROseq.DOWN:scRNAseq.DOWN`),
                   category = c("PROseq.DOWN", "scRNAseq.DOWN"), scaled = T, 
                   ext.text = F, cat.default.pos = c("text", "text"), cat.cex = 0.8, 
                   col = c("skyblue", "chocolate"), fill = c("skyblue", "chocolate"))
dev.off()
pdf(file = "venn.hct116.aid.sub.cc.Seurat.onlyNUP93.upgenes.pdf", width = 4, height = 4)
draw.pairwise.venn(area1 = length(attr(uplist, "intersections")$PROseq.UP) + length(attr(uplist, "intersections")$`PROseq.UP:scRNAseq.UP`), 
                   area2 = length(attr(uplist, "intersections")$scRNAseq.UP) + length(attr(uplist, "intersections")$`PROseq.UP:scRNAseq.UP`), 
                   cross.area = length(attr(uplist, "intersections")$`PROseq.UP:scRNAseq.UP`), 
                   category = c("PROseq.UP", "scRNAseq.UP"), scaled = T, 
                   ext.text = F, cat.default.pos = c("text", "text"), cat.cex = 0.8, 
                   col = c("skyblue", "chocolate"), fill = c("skyblue", "chocolate"))
dev.off()

# output filtered matrix --------------------------------------------------------
# table(hct116.aid.sub.cc.sub$orig.ident)
# NUP93.IAA  NUP93.UT RAD21.IAA  RAD21.UT 
# 363       186       340       306 

# table(hct116.aid.sub$orig.ident)
# NUP93.IAA  NUP93.UT RAD21.IAA  RAD21.UT 
# 607       519       344       309 
finalx1 <- hct116.aid.sub[[]][hct116.aid.sub[[]]$orig.ident=="NUP93.UT",]
finalx2 <- hct116.aid.sub.cc.sub[[]][hct116.aid.sub.cc.sub[[]]$orig.ident=="NUP93.IAA",]
finalr1 <- hct116.aid.sub[[]][hct116.aid.sub[[]]$orig.ident=="RAD21.UT",]
finalr2 <- hct116.aid.sub[[]][hct116.aid.sub[[]]$orig.ident=="RAD21.IAA",]

x1 <- fread(input = "mat2csv/X1.mat2csv.csv", header = T, sep = ",")
x1.barcodes <- gsub("-1", "", colnames(x1))
finalx1.barcodes <- gsub("nup93.ut_", "", rownames(finalx1))
finalx1.mat <- x1[, .SD, .SDcols = c(1, match(finalx1.barcodes, x1.barcodes))]
fwrite(finalx1.mat, file = "mat2csv/filtered.X1.mat2csv.csv", sep = ",")

x2 <- fread(input = "mat2csv/X2.mat2csv.csv", header = T, sep = ",")
x2.barcodes <- gsub("-1", "", colnames(x2))
finalx2.barcodes <- gsub("nup93.iaa_", "", rownames(finalx2))
finalx2.mat <- x2[, .SD, .SDcols = c(1, match(finalx2.barcodes, x2.barcodes))]
fwrite(finalx2.mat, file = "mat2csv/filtered.X2.mat2csv.csv", sep = ",")

r1 <- fread(input = "mat2csv/R1.mat2csv.csv", header = T, sep = ",")
r1.barcodes <- gsub("-1", "", colnames(r1))
finalr1.barcodes <- gsub("rad21.ut_", "", rownames(finalr1))
finalr1.mat <- r1[, .SD, .SDcols = c(1, match(finalr1.barcodes, r1.barcodes))]
fwrite(finalr1.mat, file = "mat2csv/filtered.R1.mat2csv.csv", sep = ",")

r2 <- fread(input = "mat2csv/R2.mat2csv.csv", header = T, sep = ",")
r2.barcodes <- gsub("-1", "", colnames(r2))
finalr2.barcodes <- gsub("rad21.iaa_", "", rownames(finalr2))
finalr2.mat <- r2[, .SD, .SDcols = c(1, match(finalr2.barcodes, r2.barcodes))]
fwrite(finalr2.mat, file = "mat2csv/filtered.R2.mat2csv.csv", sep = ",")

# recalculate with filtered barcodes --------------------------------------------------------
# > table(hct116.aid.sub.cc.sub$orig.ident)
# 
# NUP93.IAA  NUP93.UT RAD21.IAA  RAD21.UT 
# 363       186       340       306 
# > table(hct116.aid.sub.cc$orig.ident)
# 
# NUP93.IAA  NUP93.UT RAD21.IAA  RAD21.UT 
# 607       519       344       309 
hct116.aid.final <- subset(hct116.aid.sub.cc, subset = orig.ident != "NUP93.IAA" | seurat_clusters == 3 | seurat_clusters == 4)
# table(hct116.aid.final$orig.ident)
# 
# NUP93.IAA  NUP93.UT RAD21.IAA  RAD21.UT 
# 363       519       344       309 
finalegenes <- rownames(hct116.aid.final)
finalexprmat <- hct116.aid.final[["RNA"]]@data
finalexprmat[rownames(finalexprmat)=="KITLG",]

####
hct116.aid.final <- RunPCA(hct116.aid.final, features = VariableFeatures(object = hct116.aid.final))
hct116.aid.final <- RunUMAP(hct116.aid.final, dims = 1:10)
hct116.aid.final <- RunTSNE(hct116.aid.final, dims = 1:10)
hct116.aid.final <- FindNeighbors(hct116.aid.final, reduction = "pca", dims = 1:10, k.param = 100)
hct116.aid.final <- FindClusters(hct116.aid.final, resolution = 1.5)

print(hct116.aid.sub.cc[["pca"]], dims = 1:5, nfeatures = 5)
p4 <- VizDimLoadings(hct116.aid.final, dims = 1:2, reduction = "pca")
ggsave(filename = "hct116.aid.final.Seurat.VizDimLoadings.pca.pdf", p4, width = 10, height = 6)

p6 <- DimPlot(hct116.aid.final, reduction = "pca", group.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.pca.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "pca", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.pca.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "pca", group.by = "seurat_clusters", shape.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.pca.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "pca", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.pca.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

p6 <- DimPlot(hct116.aid.final, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.umap.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "umap", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.umap.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "umap", group.by = "seurat_clusters", shape.by= "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.umap.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.umap.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

p6 <- DimPlot(hct116.aid.final, reduction = "tsne", group.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.tsne.by.orig.ident.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "tsne", group.by = "seurat_clusters")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.tsne.by.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "tsne", group.by = "seurat_clusters", shape.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.tsne.by.orig.ident.seurat_clusters.pdf", p6, width = 6, height = 5)
p6 <- DimPlot(hct116.aid.final, reduction = "tsne", group.by = "seurat_clusters", split.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DimPlot.tsne.by.orig.ident.seurat_clusters.split.pdf", p6, width = 12, height = 5)

####
final.markers <- FindAllMarkers(hct116.aid.final, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
final.cc.topN <- final.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
final.markers[rownames(final.markers)=="KITLG",]
final.markers[grep("eRNA-", rownames(final.markers)),]
final.markers[grep("eRNA-54", rownames(final.markers)),]
final.markers[rownames(final.markers)=="CCAT1",]

# > hct116.aid.final[[]][hct116.aid.final[[]][,1]=="NUP93.UT",7]
# [1] 0 0 0 0 3 0 0 0 0 0 0 0 0 0 3 3 0 3 0 3 3 3 0 0 0 3 0 3 3 3 0 3 0 0 3 0 0 0 0 0 0 0 3 3 0 0 3 0 3 0 3 0 3 0 3 0
# [57] 0 3 0 0 0 0 0 0 3 0 0 0 0 0 3 0 0 3 0 0 0 0 3 0 0 0 3 0 3 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 0 0 0 0 3 3 3 0 0 0 3 3
# [113] 0 3 0 3 0 0 0 3 0 3 0 0 3 0 3 0 3 3 3 0 0 0 0 0 0 3 0 0 3 0 0 0 0 3 0 3 0 3 3 0 0 3 0 0 0 0 0 3 3 0 0 0 0 3 0 0
# [169] 0 0 0 0 3 3 0 0 0 0 0 0 0 3 3 3 3 0 3 3 3 0 3 3 0 0 3 3 0 0 3 3 0 0 0 0 0 3 0 3 3 0 0 0 0 3 3 0 0 0 0 3 3 3 3 0
# [225] 0 3 0 0 0 0 0 0 0 0 3 0 3 0 0 3 3 0 3 0 0 0 0 0 3 3 0 3 0 3 3 0 3 3 0 3 0 3 0 0 0 0 3 0 0 3 0 3 3 0 3 3 0 0 0 3
# [281] 0 3 0 0 0 3 3 3 0 3 3 0 0 0 3 0 3 0 3 0 0 0 0 0 0 0 3 3 3 0 0 0 0 0 3 0 0 0 0 3 3 0 3 3 0 3 0 0 3 3 3 0 3 0 3 0
# [337] 0 0 0 3 0 0 0 3 0 0 0 0 0 0 0 3 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 3 0 0 0 4 0 0 3 3 0 3 3 3 3 0 0 3 0 0
# [393] 0 0 0 0 0 3 3 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 3 0 0 3 3 0 0 3 0 3 0 3 0 0 0 3 0 0 3 0 0 3 3 0 0 0
# [449] 3 0 3 3 0 3 0 3 0 3 3 0 0 0 3 0 3 0 3 0 0 3 0 3 0 0 0 0 0 0 3 3 0 0 0 0 0 3 0 3 0 0 0 3 0 0 3 0 3 3 0 0 3 0 3 0
# [505] 3 0 3 3 0 0 3 3 0 0 3 0 0 0 0
# Levels: 0 1 2 3 4 5
# > hct116.aid.final[[]][hct116.aid.final[[]][,1]=="NUP93.IAA",7]
# [1] 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 1 3 1 1 1 1 1 1 1 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 3 1 1 1 1 3 1 1 1
# [57] 3 1 1 1 0 3 1 3 1 1 1 1 1 1 1 1 1 1 1 3 1 1 1 3 1 1 1 1 1 1 1 1 3 3 1 1 3 1 1 1 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1
# [113] 3 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 1 3 1 1 1 1 1 1 3 1 1 1 1 1 3 1 1 3 1 1 1 3 1 1 3
# [169] 3 1 0 1 1 3 1 1 1 3 1 1 3 1 1 1 1 1 1 1 3 1 1 1 1 1 3 3 1 3 1 1 1 1 1 1 1 1 1 1 1 3 3 1 1 1 1 1 3 3 1 1 1 3 1 1
# [225] 0 1 1 3 1 1 3 1 1 1 1 1 3 1 3 1 1 3 1 1 3 1 1 1 1 1 1 0 3 1 1 1 1 1 1 1 1 1 1 1 3 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1
# [281] 1 1 1 1 1 3 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 1 3 1 1 1 1 1 3 1 1 1 3 1 1 1 1 1 3 3 1 1 1 3 1 1 1 1
# [337] 3 1 1 1 1 1 1 1 3 3 1 3 1 1 1 1 1 1 1 3 1 1 3 3 1 3 3
# Levels: 0 1 2 3 4 5
final.markers.nup93 <- FindMarkers(hct116.aid.final, ident.1 = "1", ident.2 = "0",
                                   only.pos = FALSE, test.use = "poisson", 
                                   min.pct = 0.5, logfc.threshold = log2(1.5))
final.markers.nup93[rownames(final.markers.nup93)=="KITLG",]
final.markers.nup93[grep("eRNA-", rownames(final.markers.nup93)),]
final.markers.nup93[grep("eRNA-54", rownames(final.markers.nup93)),]
final.markers.nup93[rownames(final.markers.nup93)=="CCAT1",]
p10 <- DoHeatmap(hct116.aid.final, features = rownames(final.markers.nup93)[1:100], group.by = "orig.ident")
ggsave(filename = "hct116.aid.final.Seurat.DoHeatmap.onlyNUP93.groupbyorigident.pdf", 
       p10, width = 10, height = 16)

final.scup <- final.markers.nup93[final.markers.nup93$avg_logFC > 0,]
final.scdown <- final.markers.nup93[final.markers.nup93$avg_logFC < 0,]
# proup.genes <- proup$V1
# prodown.genes <- prodown$V1
final.scup.genes <- rownames(final.scup)[grep("eRNA", rownames(final.scup), invert = T)]
final.scdown.genes <- rownames(final.scdown)[grep("eRNA", rownames(final.scdown), invert = T)]
# library(gplots)
final.uplist <- venn(list(PROseq.UP = proup.genes, scRNAseq.UP = final.scup.genes))
final.downlist <- venn(list(PROseq.DOWN = prodown.genes, scRNAseq.DOWN = final.scdown.genes))

attr(final.downlist, "intersections")$PROseq.DOWN
attr(final.downlist, "intersections")$scRNAseq.DOWN
attr(final.uplist, "intersections")$PROseq.UP
attr(final.uplist, "intersections")$scRNAseq.UP
# library(VennDiagram)
# dev.off()
pdf(file = "venn.hct116.aid.final.Seurat.onlyNUP93.downgenes.pdf", width = 4, height = 4)
draw.pairwise.venn(area1 = length(attr(final.downlist, "intersections")$PROseq.DOWN) + length(attr(final.downlist, "intersections")$`PROseq.DOWN:scRNAseq.DOWN`), 
                   area2 = length(attr(final.downlist, "intersections")$scRNAseq.DOWN) + length(attr(final.downlist, "intersections")$`PROseq.DOWN:scRNAseq.DOWN`), 
                   cross.area = length(attr(final.downlist, "intersections")$`PROseq.DOWN:scRNAseq.DOWN`),
                   category = c("PROseq.DOWN", "scRNAseq.DOWN"), scaled = T, 
                   ext.text = F, cat.default.pos = c("text", "text"), cat.cex = 0.8, 
                   col = c("skyblue", "chocolate"), fill = c("skyblue", "chocolate"))
dev.off()
pdf(file = "venn.hct116.aid.final.Seurat.onlyNUP93.upgenes.pdf", width = 4, height = 4)
draw.pairwise.venn(area1 = length(attr(final.uplist, "intersections")$PROseq.UP) + length(attr(final.uplist, "intersections")$`PROseq.UP:scRNAseq.UP`), 
                   area2 = length(attr(final.uplist, "intersections")$scRNAseq.UP) + length(attr(final.uplist, "intersections")$`PROseq.UP:scRNAseq.UP`), 
                   cross.area = length(attr(final.uplist, "intersections")$`PROseq.UP:scRNAseq.UP`), 
                   category = c("PROseq.UP", "scRNAseq.UP"), scaled = T, 
                   ext.text = F, cat.default.pos = c("text", "text"), cat.cex = 0.8, 
                   col = c("skyblue", "chocolate"), fill = c("skyblue", "chocolate"))
dev.off()



# calculate cell-cell variation --------------------------------------------------------
# slotNames(hct116.aid.sub.cc.sub)
# [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"       "neighbors"    "reductions"   "project.name"
# [9] "misc"         "version"      "commands"     "tools"    
names(hct116.aid.sub.cc.sub)
# [1] "RNA"  "pca"  "umap" "tsne"
slot(hct116.aid.sub.cc.sub, name = "meta.data")
show(hct116.aid.sub.cc.sub)
class(hct116.aid.sub.cc.sub)

tstmat <- hct116.aid.sub.cc.sub[["RNA"]]@meta.features
tstmat[rownames(tstmat)=="KITLG",]
tstmat <- hct116.aid.sub.cc.sub[["RNA"]]@var.features
length(tstmat)


# Slots
# 
# raw.data
# The raw project data
# 
# data
# The normalized expression matrix (log-scale)
# 
# scale.data
# scaled (default is z-scoring each gene) expression matrix; used for dimmensional reduction and heatmap visualization
# 
# var.genes
# Vector of genes exhibiting high variance across single cells
genes <- rownames(hct116.aid.sub.cc.sub)
exprmat <- hct116.aid.sub.cc.sub[["RNA"]]@data
exprmat[rownames(exprmat)=="KITLG",]
i=10811
mydata <- data.table(Group = c(hct116.aid.sub.cc.sub$orig.ident, rep("All", length(hct116.aid.sub.cc.sub$orig.ident))), 
                     NormExpr = c(exprmat[rownames(exprmat)==genes[i],], exprmat[rownames(exprmat)==genes[i],]))
ggdensity(mydata, x = "NormExpr", fill = "Group", color = "Group")
boxplot.stats(mydata[Group == "NUP93.IAA"]$NormExpr)
boxplot.stats(mydata[Group == "NUP93.UT"]$NormExpr)
ggecdf(mydata, x = "NormExpr", linetype = "Group", color = "Group", palette = "aaas")
ggecdf(mydata[Group!="All"], x = "NormExpr", color = "Group", palette = "npg")
ggecdf(mydata, x = "NormExpr")


i=which(rownames(exprmat)=="MALAT1")
mydata <- data.table(Group = hct116.aid.sub.cc.sub$orig.ident, 
                     NormExpr = exprmat[rownames(exprmat)==genes[i],])
ggdensity(mydata[grep("NUP93",Group)], x = "NormExpr", fill = "Group", color = "Group")

# save.image(file = "seurat.image")
##############################################
# back from txburst
# myvargenes <- c("HMGCL","ZBTB20","NT5DC1","AIMP2","AAAS","ELP2","MT-ATP6","eRNA-21734")
myvargenes <- c("AIMP2", "AAAS", "ELP2")
RidgePlot(hct116.aid.sub.cc, features = myvargenes, idents = c(1:5), group.by = "orig.ident", ncol = 3)
myvargenes <- c("MT-ATP6","AAAS","ELP2")
RidgePlot(hct116.aid.sub.cc, features = myvargenes, idents = c(1:5), group.by = "orig.ident", ncol = 3)

p1 <- RidgePlot(hct116.aid.sub.cc, features = myvargenes, group.by = "orig.ident", ncol = 4)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.topvargenes.pdf", p1, width = 10, height = 6)
p1 <- RidgePlot(hct116.aid.sub.cc, features = myvargenes, idents = c(1:5), group.by = "orig.ident", ncol = 4)
ggsave(filename = "hct116.aid.sub.cc.Seurat.RidgePlot.topvargenes.selected.pdf", p1, width = 10, height = 6)
p1 <- FeaturePlot(hct116.aid.sub.cc, features = myvargenes, cols = c("white", "blue"), shape.by = "orig.ident", ncol = 4)
ggsave(filename = "hct116.aid.sub.cc.Seurat.FeaturePlot.topvargenes.pdf", p1, width = 10, height = 30)

###
myburstgenes <- c("HMGCL","RNF2","SMYD2","ZBTB20","MARVELD2","NT5DC1","AIMP2","PATL1","MALAT1","AAAS","RP11-532F12.5","SCRN2","COIL","ELP2","CDK5RAP1","CEP250","MISP","CYTH2","ZNF580","THOC5","MT-ATP6","eRNA-21734")
# i=which(rownames(exprmat)=="MALAT1")
for(j in 1:length(myburstgenes)){
    i <- which(rownames(exprmat)==myburstgenes[j])
    mydata <- data.table(Group = hct116.aid.sub.cc.sub$orig.ident, 
                         NormExpr = exprmat[rownames(exprmat)==genes[i],])
    p0 <- ggdensity(mydata[grep("NUP93",Group)], x = "NormExpr", fill = "Group", color = "Group")
    ggsave(filename = paste("ggdensity.onlyNUP93",myburstgenes[j], "pdf", sep = "."), 
           p0, width = 6, height = 4)
    p0 <- ggdensity(mydata, x = "NormExpr", fill = "Group", color = "Group")
    ggsave(filename = paste("ggdensity.bothNUP93.RAD21",myburstgenes[j], "pdf", sep = "."), 
           p0, width = 6, height = 4)
}

myburstgenes <- c("HMGCL","RNF2","SMYD2","ZBTB20","MARVELD2","NT5DC1","AIMP2","PATL1","MALAT1","AAAS","RP11-532F12.5","SCRN2","COIL","ELP2","CDK5RAP1","CEP250","MISP","CYTH2","ZNF580","THOC5","MT-ATP6","eRNA-21734")
# i=which(rownames(finalexprmat)=="MALAT1")
for(j in 1:length(myburstgenes)){
    i <- which(rownames(finalexprmat)==myburstgenes[j])
    mydata <- data.table(Group = hct116.aid.final$orig.ident, 
                         NormExpr = finalexprmat[rownames(finalexprmat)==finalegenes[i],])
    p0 <- ggdensity(mydata[grep("NUP93",Group)], x = "NormExpr", fill = "Group", color = "Group")
    ggsave(filename = paste("ggdensity.onlyNUP93",myburstgenes[j], "pdf", sep = "."), 
           p0, width = 6, height = 4)
    p0 <- ggdensity(mydata, x = "NormExpr", fill = "Group", color = "Group")
    ggsave(filename = paste("ggdensity.bothNUP93.RAD21",myburstgenes[j], "pdf", sep = "."), 
           p0, width = 6, height = 4)
}

