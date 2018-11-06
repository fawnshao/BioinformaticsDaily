library(igraph)
library(ggplot2)
library(reshape2)
library(data.table)
library(pheatmap)
args <- commandArgs(TRUE)
# add CpG island lenth for node feature
# setwd("~/Google Drive File Stream/My Drive/housekeeping_genes/interaction.community/igraph")
# args <- c("interaction.SRR5831489", "gene.EBV.sim.CGI.txt")
# args <- c("interaction.SRR5831490", "gene.EBV.sim.CGI.txt")
# args <- c("interaction.SRR5831511", "gene.EBV.sim.CGI.txt")
# args <- c("interaction.SRR5831512", "gene.EBV.sim.CGI.txt")
## args <- c("interaction", "nodes")
# 作者：何燕杰
# 链接：https://www.zhihu.com/question/22610633/answer/143644471
# 来源：知乎
# 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。
# 
# 点度中心性（degree）
# 设想一下，你在微信上有个账号，那么是不是意味着微信好友数量越多，那么你的社交圈子越广？
# （假设都是真实好友，不考虑微商神马的奇葩情况）比如我有20个好友，那么意味着20个结点与我相连。
# 如果你有50个好友，那么意味着你的点度中心度比我高，社交圈子比我广。
# 这个就是点度中心性的概念。
# 当然，刚才这个情况是无向图的情形，如果是有向图，需要考虑的出度和入度的问题。
# 在刚才的基础上拓展一下，假如我们要比较你在微博和微信上的点度中心度，刚才的方法是否适用？
# 如果说使用微信与微博的人数差不多，那么的确可以。
# 但是如果说用户数量不一样呢？那么我们需要考虑到去规模化的问题，这就是标准化的点度中心性的理念。
# 
# 接近中心性（closeness）
# 对于了解图论的朋友而言，最短路这个概念一定不陌生。
# 我们设想一个实际生活中的场景，比如你要建一个大型的娱乐商场，
# 你可能会希望周围的顾客到达这个商场的距离都可以尽可能地短。
# 这个就涉及到接近中心性的概念，接近中心性的值为路径长度的倒数。
# 接近中心性需要考量每个结点到其它结点的最短路的平均长度。
# 也就是说，对于一个结点而言，它距离其它结点越近，那么它的中心度越高。
# 一般来说，那种需要让尽可能多的人使用的设施，它的接近中心度一般是比较高的。
# 
# 中介中心性（betweenness）
# 这个度量很有意思。这个有点像是我们身边那种社交达人，
# 我们认识的不少朋友可能都是通过他/她认识的，这个人起到了中介的作用。
# 中介中心性指的是一个结点担任其它两个结点之间最短路的桥梁的次数。
# 一个结点充当“中介”的次数越高，它的中介中心度就越大。
# 如果要考虑标准化的问题，可以用一个结点承担最短路桥梁的次数除以所有的路径数量。

# generate graph ------------------------------------------------------------
links <- read.table(args[1], header = T, sep = "\t")
nodes <- read.table(args[2], header = T, sep = "\t", na.strings = "/")
network <- graph_from_data_frame(d = links, directed = F)
node.name <- names(V(graph = network))
node.degree <- degree(graph = network)
node.closeness <- closeness(graph = network)
node.betweenness <- betweenness(graph = network)

# basic node features ------------------------------------------------------------
node.feature <- data.table(node.name, node.degree, node.closeness, node.betweenness)
loopnodes.feature <- data.table(node.feature, nodes[match(node.name, nodes$Gene),])
loopnodes.feature[is.na(Gene), Type:="Enh"]
loopnodes.feature[is.na(Type), Type:="Other"]
loopnodes.feature[,Lengthquantile := cut(Length, breaks = quantile(Length, probs = seq(0, 1, 0.25), na.rm = T), 
                                         labels = paste("LenQ", 1:4, sep = ""), include.lowest = TRUE)]
loopnodes.feature[,CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25), na.rm = T),
                                      labels = paste("CpGQ", 1:4, sep = ""), include.lowest = TRUE)]
loopnodes.feature[,CpGRquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25), na.rm = T), 
                                       labels = paste("CpGRQ", 1:4, sep = ""), include.lowest = TRUE)]
loopnodes.feature[EBV.medianexpression > 0,
                  Exprquantile := cut(EBV.medianexpression, 
                                      breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T),
                                      labels = paste("ExprQ", 1:4, sep = ""), include.lowest = TRUE)]
loopnodes.feature[EBV.medianexpression == 0, Exprquantile := "NULL"]
loopnodes.feature[is.na(EBV.medianexpression), Exprquantile := "NULL"]
loopnodes.feature[,CGIquantile := cut(cpgIslandLenth, breaks = quantile(cpgIslandLenth, probs = seq(0, 1, 0.25), na.rm = T), 
                                      labels = paste("cgiQ", 1:4, sep = ""), include.lowest = TRUE)]
loopnodes.feature[is.na(Lengthquantile), Lengthquantile := "NULL"]
loopnodes.feature[is.na(CpGquantile), CpGquantile := "NULL"]
loopnodes.feature[is.na(CpGRquantile), CpGRquantile := "NULL"]
loopnodes.feature[is.na(CGIquantile), CGIquantile := "NULL"]
loopnodes.feature[is.na(Gene), Lengthquantile:="Enh"]
loopnodes.feature[is.na(Gene), CpGquantile:="Enh"]
loopnodes.feature[is.na(Gene), CpGRquantile:="Enh"]
loopnodes.feature[is.na(Gene), Exprquantile:="Enh"]
loopnodes.feature[is.na(Gene), CGIquantile:="Enh"]
loopnodes.feature.stats.1 <- table(loopnodes.feature[,c(9,13)])
loopnodes.feature.stats.2 <- table(loopnodes.feature[,c(9,16)])
loopnodes.feature.stats.3 <- table(loopnodes.feature[,c(9,12)])

dataMat <- loopnodes.feature.stats.1 / rowSums(loopnodes.feature.stats.1)
plotMat <- melt(dataMat)
plotMat$Type <- factor(plotMat$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
myplot <- ggplot(plotMat, aes(x = Type, y = value, fill = CpGquantile)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste(args[1], "loopnodes.feature.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

dataMat <- loopnodes.feature.stats.2 / rowSums(loopnodes.feature.stats.2)
plotMat <- melt(dataMat)
plotMat$Type <- factor(plotMat$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
myplot <- ggplot(plotMat, aes(x = Type, y = value, fill = CGIquantile)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste(args[1], "loopnodes.feature.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

dataMat <- loopnodes.feature.stats.3 / rowSums(loopnodes.feature.stats.3)
plotMat <- melt(dataMat)
plotMat$Type <- factor(plotMat$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
myplot <- ggplot(plotMat, aes(x = Type, y = value, fill = Lengthquantile)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste(args[1], "loopnodes.feature.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(loopnodes.feature$node.degree), 
                    Type = loopnodes.feature$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Type, y = node.degree, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.degree.NAasNA.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(loopnodes.feature$node.closeness), 
                    Type = loopnodes.feature$Type)
datay$Type <- factor(datay$Type, c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Type, y = node.closeness, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.closeness.NAasNA.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(loopnodes.feature$node.betweenness), 
                    Type = loopnodes.feature$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
myplot <- ggplot(datay, aes(x = Type, y = log2(node.betweenness), fill = Type)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "loopnodes.feature.node.betweenness.NAasNA.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(loopnodes.feature$node.degree), 
                    CpGquantile = loopnodes.feature$CpGquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CpGquantile, y = node.degree, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.degree.NAasNA.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(loopnodes.feature$node.closeness), 
                    CpGquantile = loopnodes.feature$CpGquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CpGquantile, y = node.closeness, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.closeness.NAasNA.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(loopnodes.feature$node.betweenness), 
                    CpGquantile = loopnodes.feature$CpGquantile)
myplot <- ggplot(datay, aes(x = CpGquantile, y = log2(node.betweenness), fill = CpGquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "loopnodes.feature.node.betweenness.NAasNA.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(loopnodes.feature$node.degree), 
                    CGIquantile = loopnodes.feature$CGIquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CGIquantile, y = node.degree, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.degree.NAasNA.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(loopnodes.feature$node.closeness), 
                    CGIquantile = loopnodes.feature$CGIquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CGIquantile, y = node.closeness, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.closeness.NAasNA.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(loopnodes.feature$node.betweenness), 
                    CGIquantile = loopnodes.feature$CGIquantile)
myplot <- ggplot(datay, aes(x = CGIquantile, y = log2(node.betweenness), fill = CGIquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "loopnodes.feature.node.betweenness.NAasNA.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(loopnodes.feature$node.degree), 
                    Lengthquantile = loopnodes.feature$Lengthquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = node.degree, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.degree.NAasNA.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(loopnodes.feature$node.closeness), 
                    Lengthquantile = loopnodes.feature$Lengthquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = node.closeness, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.closeness.NAasNA.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(loopnodes.feature$node.betweenness), 
                    Lengthquantile = loopnodes.feature$Lengthquantile)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = log2(node.betweenness), fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "loopnodes.feature.node.betweenness.NAasNA.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(loopnodes.feature$node.degree), 
                    Exprquantile = loopnodes.feature$Exprquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Exprquantile, y = node.degree, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.degree.NAasNA.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(loopnodes.feature$node.closeness), 
                    Exprquantile = loopnodes.feature$Exprquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Exprquantile, y = node.closeness, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "loopnodes.feature.node.closeness.NAasNA.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(loopnodes.feature$node.betweenness), 
                    Exprquantile = loopnodes.feature$Exprquantile)
myplot <- ggplot(datay, aes(x = Exprquantile, y = log2(node.betweenness), fill = Exprquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "loopnodes.feature.node.betweenness.NAasNA.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

# basic gene features ------------------------------------------------------------
gene.feature <- data.table(nodes, node.feature[match(nodes$Gene, node.feature$node.name),])
gene.feature[is.na(Type), Type:="Other"]
gene.feature[,Lengthquantile := cut(Length, breaks = quantile(Length, probs = seq(0, 1, 0.25), na.rm = T), 
                                    labels = paste("LenQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature[,CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25), na.rm = T),
                                 labels = paste("CpGQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature[,CpGRquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25), na.rm = T), 
                                  labels = paste("CpGRQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature[EBV.medianexpression > 0,
             Exprquantile := cut(EBV.medianexpression, 
                                 breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T),
                                 labels = paste("ExprQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature[EBV.medianexpression == 0, Exprquantile := "NULL"]
gene.feature[is.na(EBV.medianexpression), Exprquantile := "NULL"]
gene.feature[,CGIquantile := cut(cpgIslandLenth, breaks = quantile(cpgIslandLenth, probs = seq(0, 1, 0.25), na.rm = T), 
                                 labels = paste("cgiQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature[is.na(Lengthquantile), Lengthquantile := "NULL"]
gene.feature[is.na(CpGquantile), CpGquantile := "NULL"]
gene.feature[is.na(CpGRquantile), CpGRquantile := "NULL"]
gene.feature[is.na(CGIquantile), CGIquantile := "NULL"]
gene.feature.stats.1 <- table(gene.feature[,c(5,13)])
gene.feature.stats.2 <- table(gene.feature[,c(5,16)])
gene.feature.stats.3 <- table(gene.feature[,c(5,12)])

dataMat <- gene.feature.stats.1 / rowSums(gene.feature.stats.1)
plotMat <- melt(dataMat)
plotMat$Type <- factor(plotMat$Type, levels = c("HKG.1","HKG.2","Other","EBV.TSG", "other.TSG"))
myplot <- ggplot(plotMat, aes(x = Type, y = value, fill = CpGquantile)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste(args[1], "gene.feature.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

dataMat <- gene.feature.stats.2 / rowSums(gene.feature.stats.2)
plotMat <- melt(dataMat)
plotMat$Type <- factor(plotMat$Type, levels = c("HKG.1","HKG.2","Other","EBV.TSG", "other.TSG"))
myplot <- ggplot(plotMat, aes(x = Type, y = value, fill = CGIquantile)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste(args[1], "gene.feature.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

dataMat <- gene.feature.stats.3 / rowSums(gene.feature.stats.3)
plotMat <- melt(dataMat)
plotMat$Type <- factor(plotMat$Type, levels = c("HKG.1","HKG.2","Other","EBV.TSG", "other.TSG"))
myplot <- ggplot(plotMat, aes(x = Type, y = value, fill = Lengthquantile)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste(args[1], "gene.feature.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(gene.feature$node.degree), 
                    Type = gene.feature$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Type, y = node.degree, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "gene.feature.node.degree.NAasNA.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(gene.feature$node.closeness), 
                    Type = gene.feature$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Type, y = node.closeness, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "gene.feature.node.closeness.NAasNA.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(gene.feature$node.betweenness), 
                    Type = gene.feature$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
myplot <- ggplot(datay, aes(x = Type, y = log2(node.betweenness), fill = Type)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "gene.feature.node.betweenness.NAasNA.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(gene.feature$node.degree), 
                    CpGquantile = gene.feature$CpGquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CpGquantile, y = node.degree, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "gene.feature.node.degree.NAasNA.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(gene.feature$node.closeness), 
                    CpGquantile = gene.feature$CpGquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CpGquantile, y = node.closeness, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "gene.feature.node.closeness.NAasNA.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(gene.feature$node.betweenness), 
                    CpGquantile = gene.feature$CpGquantile)
myplot <- ggplot(datay, aes(x = CpGquantile, y = log2(node.betweenness), fill = CpGquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "gene.feature.node.betweenness.NAasNA.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(gene.feature$node.degree), 
                    CGIquantile = gene.feature$CGIquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CGIquantile, y = node.degree, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "gene.feature.node.degree.NAasNA.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(gene.feature$node.closeness), 
                    CGIquantile = gene.feature$CGIquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = CGIquantile, y = node.closeness, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "gene.feature.node.closeness.NAasNA.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(gene.feature$node.betweenness), 
                    CGIquantile = gene.feature$CGIquantile)
myplot <- ggplot(datay, aes(x = CGIquantile, y = log2(node.betweenness), fill = CGIquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "gene.feature.node.betweenness.NAasNA.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(gene.feature$node.degree), 
                    Lengthquantile = gene.feature$Lengthquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = node.degree, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "gene.feature.node.degree.NAasNA.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(gene.feature$node.closeness), 
                    Lengthquantile = gene.feature$Lengthquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = node.closeness, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "gene.feature.node.closeness.NAasNA.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(gene.feature$node.betweenness), 
                    Lengthquantile = gene.feature$Lengthquantile)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = log2(node.betweenness), fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "gene.feature.node.betweenness.NAasNA.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(gene.feature$node.degree), 
                    Exprquantile = gene.feature$Exprquantile)
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Exprquantile, y = node.degree, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "gene.feature.node.degree.NAasNA.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(gene.feature$node.closeness), 
                    Exprquantile = gene.feature$Exprquantile)
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Exprquantile, y = node.closeness, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "gene.feature.node.closeness.NAasNA.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(gene.feature$node.betweenness), 
                    Exprquantile = gene.feature$Exprquantile)
myplot <- ggplot(datay, aes(x = Exprquantile, y = log2(node.betweenness), fill = Exprquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "gene.feature.node.betweenness.NAasNA.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.degree = as.numeric(gene.feature$node.degree), 
                    Exprquantile = gene.feature$Exprquantile)
datay[is.na(node.degree),node.degree:=0]
ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Exprquantile, y = node.degree, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "gene.feature.node.degree.NAas0.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.closeness = as.numeric(gene.feature$node.closeness), 
                    Exprquantile = gene.feature$Exprquantile)
datay[is.na(node.closeness),node.closeness:=0]
ymin <- quantile(datay$node.closeness, probs = 0.01, na.rm = T)
ymax <- quantile(datay$node.closeness, probs = 0.99, na.rm = T)
myplot <- ggplot(datay, aes(x = Exprquantile, y = node.closeness, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "gene.feature.node.closeness.NAas0.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(node.betweenness = as.numeric(gene.feature$node.betweenness), 
                    Exprquantile = gene.feature$Exprquantile)
datay[is.na(node.betweenness),node.betweenness:=0]
myplot <- ggplot(datay, aes(x = Exprquantile, y = log2(node.betweenness), fill = Exprquantile)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "gene.feature.node.betweenness.NAas0.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

# community features ------------------------------------------------------------
network.long <- graph_from_data_frame(d = links[links[,4]=="Long",], directed = F)
long.node.name <- names(V(graph = network.long))
long.node.degree <- degree(graph = network.long)
long.node.closeness <- closeness(graph = network.long)
long.node.betweenness <- betweenness(graph = network.long)

network.thick <- graph_from_data_frame(d = links[links[,3]=="Thick",], directed = F)
thick.node.name <- names(V(graph = network.thick))
thick.node.degree <- degree(graph = network.thick)
thick.node.closeness <- closeness(graph = network.thick)
thick.node.betweenness <- betweenness(graph = network.thick)

dg <- decompose.graph(network)
allnodes <- length(V(graph = network))
node2community <- data.frame(Node = rep(NA, allnodes), 
                             Community = rep(NA, allnodes), 
                             CommunitySize = rep(NA, allnodes))
k <- 1
for(i in 1:length(dg)){
    a <- names(V(graph = dg[[i]]))
    size <- length(a)
    for(j in 1:size){
        node2community[k, ] <- c(a[j], paste("comp", i, sep = "."), size)
        k <- k + 1
    }
}

dg.long <- decompose.graph(network.long)
allnodes <- length(V(graph = network.long))
node2community.long <- data.frame(Node = rep(NA, allnodes), 
                             Community = rep(NA, allnodes), 
                             CommunitySize = rep(NA, allnodes))
k <- 1
for(i in 1:length(dg.long)){
    a <- names(V(graph = dg.long[[i]]))
    size <- length(a)
    for(j in 1:size){
        node2community.long[k, ] <- c(a[j], paste("long.comp", i, sep = "."), size)
        k <- k + 1
    }
}

dg.thick <- decompose.graph(network.thick)
allnodes <- length(V(graph = network.thick))
node2community.thick <- data.frame(Node = rep(NA, allnodes), 
                                  Community = rep(NA, allnodes), 
                                  CommunitySize = rep(NA, allnodes))
k <- 1
for(i in 1:length(dg.thick)){
    a <- names(V(graph = dg.thick[[i]]))
    size <- length(a)
    for(j in 1:size){
        node2community.thick[k, ] <- c(a[j], paste("thick.comp", i, sep = "."), size)
        k <- k + 1
    }
}

node2allcommunity <- data.frame(node2community, 
                                node2community.long[match(x = node2community[,1], table = node2community.long[,1]),], 
                                node2community.thick[match(x = node2community[,1], table = node2community.thick[,1]),])
# node2allcommunity[node2allcommunity[,5]=="long.comp.1" & !is.na(node2allcommunity[,5]),c(1:3,5:6,8:9)]
node2allcommunity.out <- data.table(node2allcommunity, nodes[match(node2allcommunity[,1], nodes$Gene),])
# write.table(node2allcommunity.out, file = paste(args[1], "igraph.node2allcommunity.out.tsv", sep = "."), 
            # quote = F, sep = "\t", row.names = F, col.names = T)
gene2allcommunity <- data.table(nodes, node2allcommunity[match(nodes$Gene, node2allcommunity[,1]),])
gene2allcommunity[is.na(Type), Type:="Other"]
gene2allcommunity[, FullCommnunitySize:= paste("Total", CommunitySize, sep = ".")]
gene2allcommunity[, LongCommnunitySize:= paste("Long", CommunitySize.1, sep = ".")]
gene2allcommunity[, ThickCommnunitySize:= paste("Thick", CommunitySize.2, sep = ".")]
gene2allcommunity[, FullCommnunityFlag:= ifelse(is.na(Community), "Lonely", "Total")]
gene2allcommunity[, LongCommnunityFlag:= ifelse(is.na(Community.1), "Lonely", "Long")]
gene2allcommunity[, ThickCommnunityFlag:= ifelse(is.na(Community.2), "Lonely", "Thick")]
gene2allcommunity.stats.1 <- table(gene2allcommunity[,c(5,20)])
gene2allcommunity.stats.2 <- table(gene2allcommunity[,c(5,21)])
gene2allcommunity.stats.3 <- table(gene2allcommunity[,c(5,22)])
gene2allcommunity.stats <- rbind(gene2allcommunity.stats.1, gene2allcommunity.stats.2, gene2allcommunity.stats.3)
rownames(gene2allcommunity.stats) <- paste(c(rep("Total",5), rep("Long",5), rep("Thick",5)), 
                                           rownames(gene2allcommunity.stats), sep = ":")
colnames(gene2allcommunity.stats) <- c("Lonely", "InCommnunity")
dataMat <- gene2allcommunity.stats / rowSums(gene2allcommunity.stats)
pdf(file = paste(args[1], "pheatmap.network.compnent.Total.Long.Thick.pdf", sep = "."), width = 6, height = 8)
pheatmap(dataMat, show_rownames = T, scale = "none")
dev.off()
plotMat <- melt(dataMat)
myplot <- ggplot(plotMat, aes(x = Var1, y = value, fill = Var2)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste(args[1], "barplot.network.compnent.Total.Long.Thick.pdf", sep = "."), width = 10, height = 8)
print(myplot)
dev.off()

# node edge attritubes ------------------------------------------------------------
# get all the edges connected to the specific nodes
# adjacent_vertices(network, c(1, 34))
# edge.attributes(graph = dg.long[[1]])
# g <-dg.long[[1]]
# idx <- match(names(V(g))[1], V(g)$name)
# E(g) [ from(idx) ]$Distance

# dg <- decompose.graph(network)
allnodes <- length(V(graph = network))
node2edges <- data.frame(Node = rep(NA, allnodes), 
                         ThickCount = rep(NA, allnodes), 
                         ThinCount = rep(NA, allnodes),
                         LongCount = rep(NA, allnodes), 
                         ShortCount = rep(NA, allnodes),
                         medPETs = rep(NA, allnodes), 
                         medDisBin = rep(NA, allnodes),
                         PPCount = rep(NA, allnodes), 
                         PECount = rep(NA, allnodes))
k <- 1
for(i in 1:length(dg)){
    a <- names(V(graph = dg[[i]]))
    size <- length(a)
    for(j in 1:size){
        e <- E(graph = dg[[i]]) [from(j)]
        m <- e$Thickness
        n <- e$Distance
        o <- e$PETs
        p <- e$DisBin
        q <- e$LoopType
        node2edges[k, ] <- c(a[j], length(m[m=="Thick"]), length(m[m=="Thin"]), 
                             length(n[n=="Long"]), length(n[n=="Short"]), median(o), median(p),
                             length(q[q=="PP"]), length(q[q=="PE"]))
        k <- k + 1
    }
}
node2edges.attr <- data.table(node2edges, nodes[match(node2edges[,1],nodes$Gene),])
node2edges.attr[is.na(Gene), Type:="Enh"]
node2edges.attr[is.na(Type), Type:="Other"]
node2edges.attr[,Lengthquantile := cut(Length, breaks = quantile(Length, probs = seq(0, 1, 0.25), na.rm = T), 
                                    labels = paste("lenQ", 1:4, sep = ""), include.lowest = TRUE)]
node2edges.attr[,CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25), na.rm = T), 
                          labels = paste("cpgQ", 1:4, sep = ""), include.lowest = TRUE)]
node2edges.attr[,CpGRquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25), na.rm = T), 
                                    labels = paste("cpgrQ", 1:4, sep = ""), include.lowest = TRUE)]
node2edges.attr[,Exprquantile := cut(EBV.medianexpression, breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T), 
                                     labels = paste("expQ", 1:4, sep = ""), include.lowest = TRUE)]
node2edges.attr[,CGIquantile := cut(cpgIslandLenth, breaks = quantile(cpgIslandLenth, probs = seq(0, 1, 0.25), na.rm = T), 
                                     labels = paste("cgiQ", 1:4, sep = ""), include.lowest = TRUE)]
node2edges.attr[is.na(Gene),Lengthquantile:="Enh"]
node2edges.attr[is.na(Gene),CpGquantile:="Enh"]
node2edges.attr[is.na(Gene),CpGRquantile:="Enh"]
node2edges.attr[is.na(Gene),Exprquantile:="Enh"]
node2edges.attr[is.na(Gene),CGIquantile:="Enh"]

# by Type ------------------------------------------------------------
datay <- data.table(ThickCount = as.numeric(node2edges.attr$ThickCount), 
                    Type = node2edges.attr$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$ThickCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = ThickCount, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "ThickCount.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(LongCount = as.numeric(node2edges.attr$LongCount), 
                    Type = node2edges.attr$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$LongCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = LongCount, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "LongCount.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medPETs = as.numeric(node2edges.attr$medPETs), 
                    Type = node2edges.attr$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$medPETs, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = medPETs, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medPETs.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medDisBin = as.numeric(node2edges.attr$medDisBin), 
                    Type = node2edges.attr$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$medDisBin, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = medDisBin, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medDisBin.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PPCount = as.numeric(node2edges.attr$PPCount), 
                    Type = node2edges.attr$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$PPCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = PPCount, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PPCount.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PECount = as.numeric(node2edges.attr$PECount), 
                    Type = node2edges.attr$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$PECount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = PECount, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PECount.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()
# by CpGquantile ------------------------------------------------------------
datay <- data.table(ThickCount = as.numeric(node2edges.attr$ThickCount), 
                    CpGquantile = node2edges.attr$CpGquantile)
ymax <- quantile(datay$ThickCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = ThickCount, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "ThickCount.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(LongCount = as.numeric(node2edges.attr$LongCount), 
                    CpGquantile = node2edges.attr$CpGquantile)
ymax <- quantile(datay$LongCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = LongCount, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "LongCount.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medPETs = as.numeric(node2edges.attr$medPETs), 
                    CpGquantile = node2edges.attr$CpGquantile)
ymax <- quantile(datay$medPETs, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = medPETs, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medPETs.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medDisBin = as.numeric(node2edges.attr$medDisBin), 
                    CpGquantile = node2edges.attr$CpGquantile)
ymax <- quantile(datay$medDisBin, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = medDisBin, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medDisBin.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PPCount = as.numeric(node2edges.attr$PPCount), 
                    CpGquantile = node2edges.attr$CpGquantile)
ymax <- quantile(datay$PPCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = PPCount, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PPCount.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PECount = as.numeric(node2edges.attr$PECount), 
                    CpGquantile = node2edges.attr$CpGquantile)
ymax <- quantile(datay$PECount, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = PECount, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PECount.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()
# by Lengthquantile ------------------------------------------------------------
datay <- data.table(ThickCount = as.numeric(node2edges.attr$ThickCount), 
                    Lengthquantile = node2edges.attr$Lengthquantile)
ymax <- quantile(datay$ThickCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = ThickCount, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "ThickCount.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(LongCount = as.numeric(node2edges.attr$LongCount), 
                    Lengthquantile = node2edges.attr$Lengthquantile)
ymax <- quantile(datay$LongCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = LongCount, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "LongCount.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medPETs = as.numeric(node2edges.attr$medPETs), 
                    Lengthquantile = node2edges.attr$Lengthquantile)
ymax <- quantile(datay$medPETs, probs = 0.99)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = medPETs, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medPETs.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medDisBin = as.numeric(node2edges.attr$medDisBin), 
                    Lengthquantile = node2edges.attr$Lengthquantile)
ymax <- quantile(datay$medDisBin, probs = 0.99)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = medDisBin, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medDisBin.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PPCount = as.numeric(node2edges.attr$PPCount), 
                    Lengthquantile = node2edges.attr$Lengthquantile)
ymax <- quantile(datay$PPCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = PPCount, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PPCount.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PECount = as.numeric(node2edges.attr$PECount), 
                    Lengthquantile = node2edges.attr$Lengthquantile)
ymax <- quantile(datay$PECount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Lengthquantile, y = PECount, fill = Lengthquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PECount.by.Lengthquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()
# by Exprquantile ------------------------------------------------------------
datay <- data.table(ThickCount = as.numeric(node2edges.attr$ThickCount), 
                    Exprquantile = node2edges.attr$Exprquantile)
ymax <- quantile(datay$ThickCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Exprquantile, y = ThickCount, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "ThickCount.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(LongCount = as.numeric(node2edges.attr$LongCount), 
                    Exprquantile = node2edges.attr$Exprquantile)
ymax <- quantile(datay$LongCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Exprquantile, y = LongCount, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "LongCount.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medPETs = as.numeric(node2edges.attr$medPETs), 
                    Exprquantile = node2edges.attr$Exprquantile)
ymax <- quantile(datay$medPETs, probs = 0.99)
myplot <- ggplot(datay, aes(x = Exprquantile, y = medPETs, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medPETs.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medDisBin = as.numeric(node2edges.attr$medDisBin), 
                    Exprquantile = node2edges.attr$Exprquantile)
ymax <- quantile(datay$medDisBin, probs = 0.99)
myplot <- ggplot(datay, aes(x = Exprquantile, y = medDisBin, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medDisBin.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PPCount = as.numeric(node2edges.attr$PPCount), 
                    Exprquantile = node2edges.attr$Exprquantile)
ymax <- quantile(datay$PPCount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Exprquantile, y = PPCount, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PPCount.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PECount = as.numeric(node2edges.attr$PECount), 
                    Exprquantile = node2edges.attr$Exprquantile)
ymax <- quantile(datay$PECount, probs = 0.99)
myplot <- ggplot(datay, aes(x = Exprquantile, y = PECount, fill = Exprquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PECount.by.Exprquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()
# by CGIquantile ------------------------------------------------------------
datay <- data.table(ThickCount = as.numeric(node2edges.attr$ThickCount), 
                    CGIquantile = node2edges.attr$CGIquantile)
ymax <- quantile(datay$ThickCount, probs = 0.95)
myplot <- ggplot(datay, aes(x = CGIquantile, y = ThickCount, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "ThickCount.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(LongCount = as.numeric(node2edges.attr$LongCount), 
                    CGIquantile = node2edges.attr$CGIquantile)
ymax <- quantile(datay$LongCount, probs = 0.95)
myplot <- ggplot(datay, aes(x = CGIquantile, y = LongCount, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "LongCount.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medPETs = as.numeric(node2edges.attr$medPETs), 
                    CGIquantile = node2edges.attr$CGIquantile)
ymax <- quantile(datay$medPETs, probs = 0.95)
myplot <- ggplot(datay, aes(x = CGIquantile, y = medPETs, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medPETs.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(medDisBin = as.numeric(node2edges.attr$medDisBin), 
                    CGIquantile = node2edges.attr$CGIquantile)
ymax <- quantile(datay$medDisBin, probs = 0.99)
myplot <- ggplot(datay, aes(x = CGIquantile, y = medDisBin, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "medDisBin.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PPCount = as.numeric(node2edges.attr$PPCount), 
                    CGIquantile = node2edges.attr$CGIquantile)
ymax <- quantile(datay$PPCount, probs = 0.95)
myplot <- ggplot(datay, aes(x = CGIquantile, y = PPCount, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PPCount.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(PECount = as.numeric(node2edges.attr$PECount), 
                    CGIquantile = node2edges.attr$CGIquantile)
ymax <- quantile(datay$PECount, probs = 0.95)
myplot <- ggplot(datay, aes(x = CGIquantile, y = PECount, fill = CGIquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "PECount.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

table(node2edges.attr[PECount > 0]$Type) / table(node2edges.attr$Type)
table(node2edges.attr[PPCount > 0]$Type) / table(node2edges.attr$Type)
table(node2edges.attr[LongCount > 0]$Type) / table(node2edges.attr$Type)
table(node2edges.attr[ThickCount > 0]$Type) / table(node2edges.attr$Type)

# Feature Selection Based on looping data ------------------------------------------------------------
node2edges.attr[,Thick.percetage:=as.numeric(ThickCount)/(as.numeric(ThickCount)+as.numeric(ThinCount))]
node2edges.attr[,Long.percetage:=as.numeric(LongCount)/(as.numeric(LongCount)+as.numeric(ShortCount))]
node2edges.attr[,PP.percetage:=as.numeric(PPCount)/(as.numeric(PPCount)+as.numeric(PECount))]

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
RunFS <- function(x, y = as.numeric(types), pre = args[1]){
    require(FeatureSelection)
    trainx <- x
    trainy <- y
    pre <- pre
    params_glmnet <- list(alpha = 1, family = 'gaussian', nfolds = 5, parallel = TRUE)
    params_xgboost <- list(params = list("objective" = "reg:linear", "bst:eta" = 0.001, 
                                         "subsample" = 0.75, "max_depth" = 5, 
                                         "colsample_bytree" = 0.75, "nthread" = 60), 
                           nrounds = 1000, print.every.n = 250, maximize = FALSE)
    params_ranger <- list(dependent.variable.name = 'y', probability = FALSE, 
                          num.trees = 1000, verbose = TRUE, mtry = 5, 
                          min.node.size = 10, num.threads = 60, 
                          classification = FALSE, importance = 'permutation')
    params_features <- list(keep_number_feat = NULL, union = TRUE)
    feat <- wrapper_feat_select(X = trainx, y = trainy, 
                                params_glmnet = params_glmnet, 
                                params_xgboost = params_xgboost, 
                                params_ranger = params_ranger, 
                                xgb_sort = 'Gain', CV_folds = 5, 
                                stratified_regr = FALSE, scale_coefs_glmnet = FALSE, 
                                cores_glmnet = 5, params_features = params_features, verbose = TRUE)
    write.table(feat$all_feat$"glmnet-lasso", paste(pre, "FeatureSelection.glmnetlasso.tsv", sep = "."), row.names = F, sep = "\t")
    write.table(feat$all_feat$xgboost, paste(pre, "FeatureSelection.xgboost.tsv", sep = "."), row.names = F, sep = "\t")
    write.table(feat$all_feat$ranger, paste(pre, "FeatureSelection.ranger.tsv", sep = "."), row.names = F, sep = "\t")
    write.table(feat$union_feat, paste(pre, "FeatureSelection.union_feat.tsv", sep = "."), row.names = F, sep = "\t")
    return(feat)
}
traindata <- node2edges.attr[Type!="Enh",c(2:9,11:13,15:16)]
types <- node2edges.attr[Type!="Enh"]$Type
traindata <- data.matrix(traindata)
traindata[is.na(traindata)] <- 0
traindata <- apply(traindata, 2, range01)
traindata <- data.frame(traindata, data.matrix(node2edges.attr[Type!="Enh",c(22:24)]))


node2edges.feat <- RunFS(x = traindata, y = as.numeric(types), 
                         pre = paste(args[1],"FeatureSelection.importance.tsv", sep = "."))

RunPCA <- function(x, y = types, pre = args[1]){
    require(ggplot2)
    data <- x
    mycolors <- y
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
    
    pca.scores <- as.data.frame(myPCA.a.t$x)
    pdf(file = paste(pre, "t.pca.pdf", sep = "."), width = 15, height = 12)
    p1 <- ggplot(data = pca.scores, aes(x = PC1, y = PC2, label = rownames(pca.scores))) +
        geom_point(aes(alpha = 1/10)) +
        geom_text(colour = "tomato", alpha = 0.8, size = 2) +
        ggtitle(args[1])
    print(p1)
    dev.off()
}
node2edges.pca <- RunPCA(x = traindata, y = types, pre = paste(args[1],"node2edges", sep = "."))


set.seed(123)
RuntSNE <- function(x, pre = args[1]){
    # require(tsne)
    require(Rtsne)
    data <- x
    pre <- pre
    # tsne <- tsne(data, initial_dims = 2, perplexity = 30, max_iter = 500)
    tsne <- Rtsne(data, check_duplicates = FALSE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
    v.tsne <- as.data.frame(tsne$Y)
    colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
    write.table(v.tsne, paste(pre, "tSNE.tsv", sep = "."), row.names = T, sep = "\t")
    return(tsne)
}
node2edges.tsne <- RuntSNE(x = traindata, pre = paste(args[1], "node2edges", sep = "."))

PlottSNE <- function(x, y = types, pre = args[1]){
    require(ggplot2)
    data <- x
    mycolors <- y
    pre <- pre
    v.tsne <- as.data.frame(data$Y)
    colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
    pdf(file = paste(pre, "tsne.pdf", sep = "."), width = 15, height = 12)
    p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
        geom_point(aes(colour = mycolors, alpha = 1/10)) +
        ggtitle(pre)
    print(p1)
    dev.off()
}
node2edges.tsne.plot <- PlottSNE(x = node2edges.tsne, pre = paste(args[1],"node2edges", sep = "."))

node2edges.tsne.1 <- RuntSNE(x = traindata[types!="Other",], pre = paste(args[1], "node2edges.1", sep = "."))
PlottSNE(x = node2edges.tsne.1, y = types[types!="Other"], pre = paste(args[1], "node2edges.1", sep = "."))
PlottSNE(x = node2edges.tsne, y = node2edges.attr[Type!="Enh"]$CpGquantile, pre = paste(args[1],"node2edges.2", sep = "."))
PlottSNE(x = node2edges.tsne, y = node2edges.attr[Type!="Enh"]$Lengthquantile, pre = paste(args[1],"node2edges.3", sep = "."))
PlottSNE(x = node2edges.tsne, y = node2edges.attr[Type!="Enh"]$Exprquantile, pre = paste(args[1],"node2edges.4", sep = "."))
PlottSNE(x = node2edges.tsne, y = node2edges.attr[Type!="Enh"]$CGIquantile, pre = paste(args[1],"node2edges.5", sep = "."))

node2edges.tsne.2 <- RuntSNE(x = traindata[,5:16], pre = paste(args[1], "node2edges.2", sep = "."))
PlottSNE(x = node2edges.tsne.2, y = node2edges.attr[Type!="Enh"]$Type, pre = paste(args[1],"node2edges.2.1", sep = "."))
PlottSNE(x = node2edges.tsne.2, y = node2edges.attr[Type!="Enh"]$CpGquantile, pre = paste(args[1],"node2edges.2.2", sep = "."))
PlottSNE(x = node2edges.tsne.2, y = node2edges.attr[Type!="Enh"]$Lengthquantile, pre = paste(args[1],"node2edges.2.3", sep = "."))
PlottSNE(x = node2edges.tsne.2, y = node2edges.attr[Type!="Enh"]$Exprquantile, pre = paste(args[1],"node2edges.2.4", sep = "."))
PlottSNE(x = node2edges.tsne.2, y = node2edges.attr[Type!="Enh"]$CGIquantile, pre = paste(args[1],"node2edges.2.5", sep = "."))


# Fature Selection 2 ------------------------------------------------------------
node.feature <- data.table(node.name, node.degree, node.closeness, node.betweenness)
gene2communityfeature <- data.table(nodes, node2community[match(nodes$Gene, node2community[,1]),],
                                    node.feature[match(nodes$Gene, node.feature$node.name),],
                                    node2edges[match(nodes$Gene, node2edges$Node),])
gene2communityfeature[is.na(Type), Type:="Other"]
gene2communityfeature[,Lengthquantile := cut(Length, breaks = quantile(Length, probs = seq(0, 1, 0.25), na.rm = T), 
                                       labels = paste("LenQ", 1:4, sep = ""), include.lowest = TRUE)]
gene2communityfeature[,CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25), na.rm = T), 
                                    labels = paste("CpGQ", 1:4, sep = ""), include.lowest = TRUE)]
gene2communityfeature[,CpGRquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25), na.rm = T), 
                                     labels = paste("CpGRQ", 1:4, sep = ""), include.lowest = TRUE)]
gene2communityfeature[EBV.medianexpression > 0,
                      Exprquantile := cut(EBV.medianexpression, 
                                          breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T),
                                          labels = paste("ExprQ", 1:4, sep = ""), include.lowest = TRUE)]
gene2communityfeature[,CGIquantile := cut(cpgIslandLenth, breaks = quantile(cpgIslandLenth, probs = seq(0, 1, 0.25), na.rm = T), 
                                    labels = paste("cgiQ", 1:4, sep = ""), include.lowest = TRUE)]
gene2communityfeature[EBV.medianexpression == 0, Exprquantile := "NULL"]
gene2communityfeature[is.na(Lengthquantile), Lengthquantile := "NULL"]
gene2communityfeature[is.na(CpGquantile), CpGquantile := "NULL"]
gene2communityfeature[is.na(CpGRquantile), CpGRquantile := "NULL"]
gene2communityfeature[is.na(CGIquantile), CGIquantile := "NULL"]
gene2communityfeature[,edgeFlag := cut(as.numeric(gene2communityfeature$PPCount), 
                                       breaks = quantile(as.numeric(gene2communityfeature$PPCount), 
                                                         probs = seq(0, 1, 0.25), na.rm = T), 
                                       labels = paste("edgeQ", 1:4, sep = ""), include.lowest = TRUE)]
gene2communityfeature[is.na(edgeFlag), edgeFlag := "NULL"]
gene2communityfeature[,Thick.percetage:=as.numeric(ThickCount)/(as.numeric(ThickCount)+as.numeric(ThinCount))]
gene2communityfeature[,Long.percetage:=as.numeric(LongCount)/(as.numeric(LongCount)+as.numeric(ShortCount))]
gene2communityfeature[,PP.percetage:=as.numeric(PPCount)/(as.numeric(PPCount)+as.numeric(PECount))]

trainD <- data.matrix(gene2communityfeature[,c(2:4,6:7,10,12:14,16:23)])
types <- gene2communityfeature$Type
trainD[is.na(trainD)] <- 0
trainD <- apply(trainD, 2, range01)
trainD <- data.frame(trainD, data.matrix(gene2communityfeature[,c(30:32)]))
trainD[is.na(trainD)] <- 0

gene2communityfeature.feat <- RunFS(x = trainD, y = as.numeric(types), 
                         pre = paste(args[1],"allgenes.importance.tsv", sep = "."))
gene2communityfeature.tsne <- RuntSNE(x = trainD, pre = paste(args[1], "allgenes.communityfeature", sep = "."))
PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$Type, pre = paste(args[1],"allgenes.communityfeature.Type", sep = "."))
PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$CpGquantile, pre = paste(args[1],"allgenes.communityfeature.CpG", sep = "."))
PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$Lengthquantile, pre = paste(args[1],"allgenes.communityfeature.Length", sep = "."))
PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$Exprquantile, pre = paste(args[1],"allgenes.communityfeature.Expr", sep = "."))
PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$CGIquantile, pre = paste(args[1],"allgenes.communityfeature.CGI", sep = "."))
v.tsne <- as.data.frame(gene2communityfeature.tsne$Y)[types!="Other",]
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
pdf(file = paste(args[1],"allgenes.communityfeature.Type.1", "tsne.pdf", sep = "."), width = 15, height = 12)
p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
    geom_point(aes(colour = gene2communityfeature[Type!="Other"]$Type, alpha = 1/10))
print(p1)
dev.off()
v.tsne <- as.data.frame(gene2communityfeature.tsne$Y)[types!="Other" & types!="other.TSG",]
colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
pdf(file = paste(args[1],"allgenes.communityfeature.Type.2", "tsne.pdf", sep = "."), width = 15, height = 12)
p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
    geom_point(aes(colour = gene2communityfeature[Type!="Other"][Type!="other.TSG"]$Type, alpha = 1/10))
print(p1)
dev.off()
# v.tsne <- as.data.frame(gene2communityfeature.tsne$Y)[types!="Other",]
# colnames(v.tsne) <- c("tSNE.1", "tSNE.2")
# gene2communityfeature[,edgeFlag := cut(as.numeric(gene2communityfeature$PECount), 
#                                        breaks = quantile(as.numeric(gene2communityfeature$PECount), 
#                                                          probs = seq(0, 1, 0.25), na.rm = T), 
#                                        labels = paste("edgeQ", 1:4, sep = ""), include.lowest = TRUE)]
# gene2communityfeature[is.na(edgeFlag), edgeFlag := "NULL"]
# # pdf(file = paste(args[1],"allgenes.communityfeature.Type.1", "tsne.pdf", sep = "."), width = 15, height = 12)
# p1 <- ggplot(data = v.tsne, aes(x = tSNE.1, y = tSNE.2)) + 
#     geom_point(aes(colour = gene2communityfeature[Type!="Other"]$edgeFlag, alpha = 1/10))
# print(p1)
# # dev.off()


gene2communityfeature.feat <- RunFS(x = trainD[types!="Other",], y = as.numeric(types[types!="Other"]), 
                                    pre = paste(args[1],"allgenes.importance.1.tsv", sep = "."))
gene2communityfeature.tsne.1 <- RuntSNE(x = trainD[types!="Other",], pre = paste(args[1], "allgenes.communityfeature.1", sep = "."))
PlottSNE(x = gene2communityfeature.tsne.1, y = gene2communityfeature[Type!="Other"]$Type, pre = paste(args[1],"allgenes.communityfeature.1.Type", sep = "."))
PlottSNE(x = gene2communityfeature.tsne.1, y = gene2communityfeature[Type!="Other"]$CpGquantile, pre = paste(args[1],"allgenes.communityfeature.1.CpG", sep = "."))
PlottSNE(x = gene2communityfeature.tsne.1, y = gene2communityfeature[Type!="Other"]$Lengthquantile, pre = paste(args[1],"allgenes.communityfeature.1.Length", sep = "."))
PlottSNE(x = gene2communityfeature.tsne.1, y = gene2communityfeature[Type!="Other"]$Exprquantile, pre = paste(args[1],"allgenes.communityfeature.1.Expr", sep = "."))
PlottSNE(x = gene2communityfeature.tsne.1, y = gene2communityfeature[Type!="Other"]$CGIquantile, pre = paste(args[1],"allgenes.communityfeature.1.CGI", sep = "."))

# # Fature Selection 3 ------------------------------------------------------------
# trainDS <- data.matrix(gene2communityfeature[,c(30:32)])
# types <- gene2communityfeature$Type
# trainDS[is.na(trainDS)] <- 0
# 
# # gene2edgep.feat <- RunFS(x = trainDS, y = as.numeric(types), pre = paste(args[1],"allgenes.importance.2.tsv", sep = "."))
# gene2edgep.tsne <- RuntSNE(x = trainDS, pre = paste(args[1], "allgenes.edgePfeature.2", sep = "."))
# PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$Type, pre = paste(args[1],"allgenes.edgePfeature.Type", sep = "."))
# PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$CpGquantile, pre = paste(args[1],"allgenes.edgePfeature.CpG", sep = "."))
# PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$Lengthquantile, pre = paste(args[1],"allgenes.edgePfeature.Length", sep = "."))
# PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$Exprquantile, pre = paste(args[1],"allgenes.edgePfeature.Expr", sep = "."))
# PlottSNE(x = gene2communityfeature.tsne, y = gene2communityfeature$CGIquantile, pre = paste(args[1],"allgenes.edgePfeature.CGI", sep = "."))

# PP : PE ------------------------------------------------------------
datatest <- node2edges.attr[,c(8:9,14,18,21)]
datatest[,PP.percetage:=as.numeric(PPCount)/(as.numeric(PPCount)+as.numeric(PECount))]
datatest[,PE.percetage:=as.numeric(PECount)/(as.numeric(PPCount)+as.numeric(PECount))]
zz <- melt(datatest[Type!="Enh",c(3,6:7)], id.vars = "Type")
zz$Type <- factor(zz$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
myplot <- ggplot(data = zz, aes(x = Type, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.PP.PE.ratio.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = Type)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.PP.PE.ratio.by.Type.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

zz <- melt(datatest[Type!="Enh",c(4,6:7)], id.vars = "CpGquantile")
myplot <- ggplot(data = zz, aes(x = CpGquantile, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.PP.PE.ratio.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = CpGquantile)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.PP.PE.ratio.by.CpGquantile.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

zz <- melt(datatest[Type!="Enh",c(5,6:7)], id.vars = "CGIquantile")
myplot <- ggplot(data = zz, aes(x = CGIquantile, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.PP.PE.ratio.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = CGIquantile)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.PP.PE.ratio.by.CGIquantile.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

# Thick : Thin ------------------------------------------------------------
datatest <- node2edges.attr[,c(2:3,14,18,21)]
datatest[,Thick.percetage:=as.numeric(ThickCount)/(as.numeric(ThickCount)+as.numeric(ThinCount))]
datatest[,Thin.percetage:=as.numeric(ThinCount)/(as.numeric(ThickCount)+as.numeric(ThinCount))]
zz <- melt(datatest[Type!="Enh",c(3,6:7)], id.vars = "Type")
zz$Type <- factor(zz$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
myplot <- ggplot(data = zz, aes(x = Type, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Thickness.ratio.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = Type)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Thickness.ratio.by.Type.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

zz <- melt(datatest[Type!="Enh",c(4,6:7)], id.vars = "CpGquantile")
myplot <- ggplot(data = zz, aes(x = CpGquantile, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Thickness.ratio.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = CpGquantile)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Thickness.ratio.by.CpGquantile.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

zz <- melt(datatest[Type!="Enh",c(5,6:7)], id.vars = "CGIquantile")
myplot <- ggplot(data = zz, aes(x = CGIquantile, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Thickness.ratio.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = CGIquantile)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Thickness.ratio.by.CGIquantile.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

# Long : Short ------------------------------------------------------------
datatest <- node2edges.attr[,c(4:5,14,18,21)]
datatest[,Long.percetage:=as.numeric(LongCount)/(as.numeric(LongCount)+as.numeric(ShortCount))]
datatest[,Short.percetage:=as.numeric(ShortCount)/(as.numeric(LongCount)+as.numeric(ShortCount))]
zz <- melt(datatest[Type!="Enh",c(3,6:7)], id.vars = "Type")
zz$Type <- factor(zz$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
myplot <- ggplot(data = zz, aes(x = Type, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Distance.ratio.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = Type)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Distance.ratio.by.Type.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

zz <- melt(datatest[Type!="Enh",c(4,6:7)], id.vars = "CpGquantile")
myplot <- ggplot(data = zz, aes(x = CpGquantile, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Distance.ratio.by.CpGquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = CpGquantile)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Distance.ratio.by.CpGquantile.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

zz <- melt(datatest[Type!="Enh",c(5,6:7)], id.vars = "CGIquantile")
myplot <- ggplot(data = zz, aes(x = CGIquantile, y = value, fill = variable)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Distance.ratio.by.CGIquantile.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(data = zz, aes(x = variable, y = value, fill = CGIquantile)) + geom_boxplot()
pdf(file = paste(args[1], "LoopedNode.Distance.ratio.by.CGIquantile.1.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

datatest <- node2edges.attr[Type!="Enh",c(2:5,8:9,14,18)]
datatest[,Thick.percetage:=as.numeric(ThickCount)/(as.numeric(ThickCount)+as.numeric(ThinCount))]
datatest[,Long.percetage:=as.numeric(LongCount)/(as.numeric(LongCount)+as.numeric(ShortCount))]
datatest[,PP.percetage:=as.numeric(PPCount)/(as.numeric(PPCount)+as.numeric(PECount))]
datatest[Type=="HKG.1",Type:="HKG"]
datatest[Type=="HKG.2",Type:="HKG"]
datam <- data.matrix(datatest[Type!="Other",9:11])
rownames(datam) <- 1:nrow(datam)
annos_row <- data.frame(Type=datatest[Type!="Other"]$Type, CpGquantile=datatest[Type!="Other"]$CpGquantile)
rownames(annos_row) <- rownames(datam)
datam <- datam[order(datam[,3],datam[,2],datam[,1],decreasing = T),]
annos_row <- annos_row[order(datam[,3],datam[,2],datam[,1],decreasing = T),]
pdf(file = paste(args[1], "LoopedNode.percentage.all.pdf", sep = "."), width = 6, height = 10)
pheatmap(datam, show_rownames = F, scale = "none", annotation_row = annos_row, cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "LoopedNode.percentage.HKG.pdf", sep = "."), width = 6, height = 10)
pheatmap(datam[annos_row$Type=="HKG",], show_rownames = F, scale = "none", cluster_rows = F)
dev.off()
pdf(file = paste(args[1], "LoopedNode.percentage.nonHKG.pdf", sep = "."), width = 6, height = 10)
pheatmap(datam[annos_row$Type!="HKG",], show_rownames = F, scale = "none", cluster_rows = F)
dev.off()
# ggplot(data = datatest[Type!="Enh"][Type!="Other"], aes(x = PP.percetage, y = as.numeric(medPETs), colour = CpGquantile)) + geom_point()
save.image(file = paste(args[1], "RData", sep = "."))
