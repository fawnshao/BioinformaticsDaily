library(igraph)
library(ggplot2)
library(reshape2)
library(data.table)
library(pheatmap)
args <- commandArgs(TRUE)
# awk '$NF!="PE"' interaction.SRR5831489 > interaction.SRR5831489.pp
# awk '$NF!="PP"' interaction.SRR5831489 > interaction.SRR5831489.pe
# args <- c("interaction.SRR5831489.pe", "gene.EBV.sim.txt")
# setwd("~/Google Drive File Stream/My Drive/housekeeping_genes/interaction.community/igraph")
# args <- c("interaction.SRR5831489", "gene.EBV.sim.txt")
# args <- c("interaction.SRR5831490", "gene.EBV.sim.txt")
# args <- c("interaction.SRR5831511", "gene.EBV.sim.txt")
# args <- c("interaction.SRR5831512", "gene.EBV.sim.txt")
## args <- c("interaction", "nodes")
## args <- c("interaction.SRR5831489.HKG.links", "interaction.SRR5831489.HKG.links.nodes")
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

links <- read.table(args[1], header = T, sep = "\t")
nodes <- read.table(args[2], header = T, sep = "\t", na.strings = "/")
network <- graph_from_data_frame(d = links, directed = F)
node.name <- names(V(graph = network))
node.degree <- degree(graph = network)
node.closeness <- closeness(graph = network)
# node.closeness.1 <- closeness(graph = network, weights = links$PETs)
node.betweenness <- betweenness(graph = network)
nodestats <- data.frame(NodeDegree = node.degree, 
                        NodeCloseness = node.closeness,
                        NodeBetweenness = node.betweenness)
# betweeness centrality for each node for grouping
# vertices$group <- edge.betweenness.community(network)$membership
# returns a list of three graphs
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
# node2community

output <- data.frame(node2community, nodestats[match(node2community[,1], node.name),],
                     nodes[match(node2community[,1], nodes$Gene),])
write.table(output, file = paste(args[1], "igraph.out.tsv", sep = "."), 
            quote = F, sep = "\t", row.names = F, col.names = T)


summary(as.numeric(output$CommunitySize))
datastats <- as.data.table(output)
datastats[is.na(Gene), Type:="Enh"]
datastats[is.na(Type), Type:="other"]

datay <- data.table(NodeDegree = as.numeric(datastats$NodeDegree), 
                    Type = datastats$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$NodeDegree, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = NodeDegree, fill = Type)) +
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "NodeDegree.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(NodeDegree = as.numeric(datastats$NodeDegree), 
                    CpGCount = datastats$CpG.Count)
datay[,CpGquantile := cut(CpGCount, breaks = quantile(CpGCount, probs = seq(0, 1, 0.25), na.rm = T), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
datay[is.na(CpGquantile), CpGquantile:="Enh"]
ymax <- quantile(datay$NodeDegree, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = NodeDegree, fill = CpGquantile)) +
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste(args[1], "NodeDegree.by.CpGCount.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

datay <- data.table(NodeCloseness = as.numeric(datastats$NodeCloseness), 
                    Type = datastats$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$NodeCloseness, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = NodeCloseness, fill = Type)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "NodeCloseness.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(NodeCloseness = as.numeric(datastats$NodeCloseness), 
                    CpGCount = datastats$CpG.Count)
datay[,CpGquantile := cut(CpGCount, breaks = quantile(CpGCount, probs = seq(0, 1, 0.25), na.rm = T), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
datay[is.na(CpGquantile), CpGquantile:="Enh"]
ymin <- quantile(datay$NodeCloseness, probs = 0.01)
ymax <- quantile(datay$NodeCloseness, probs = 0.99)
myplot <- ggplot(datay, aes(x = CpGquantile, y = NodeCloseness, fill = CpGquantile)) +
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "NodeCloseness.by.CpGCount.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

datay <- data.table(NodeBetweenness = as.numeric(datastats$NodeBetweenness), 
                    Type = datastats$Type)
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$NodeBetweenness, probs = 0.99)
myplot <- ggplot(datay, aes(x = Type, y = log2(NodeBetweenness), fill = Type)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "NodeBetweenness.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)
dev.off()

datay <- data.table(NodeBetweenness = as.numeric(datastats$NodeBetweenness), 
                    CpGCount = datastats$CpG.Count)
datay[,CpGquantile := cut(CpGCount, breaks = quantile(CpGCount, probs = seq(0, 1, 0.25), na.rm = T), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
datay[is.na(CpGquantile), CpGquantile:="Enh"]
ymin <- quantile(datay$NodeBetweenness, probs = 0.1)
ymax <- quantile(datay$NodeBetweenness, probs = 0.7)
myplot <- ggplot(datay, aes(x = CpGquantile, y = log2(NodeBetweenness), fill = CpGquantile)) +
    geom_boxplot(width = 0.8) #+ scale_y_continuous(limits = c(ymin, ymax))
pdf(file = paste(args[1], "NodeBetweenness.by.CpGCount.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()


datay <- datastats[,c(2,3,9,11)]
datay[,CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25), na.rm = T), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
datay[is.na(CpGquantile), CpGquantile:="Enh"]
datay[TypeCpG=="Enh.Enh",TypeCpG:="Enh"]
CountMat <- table(datay[CommunitySize >= 3, c(1,6)])
pdf(file = paste(args[1], "pheatmap.Community.detail.pdf", sep = "."), width = 6, height = 8)
pheatmap(CountMat, show_rownames = F, scale = "row") #, color = colorRampPalette(c("blue", "white", "red"))(11))
dev.off()

CountMat <- table(datastats[CommunitySize >= 3, c(2,11)])
# CountMat[CountMat==0] <- NA
pdf(file = paste(args[1], "pheatmap.Community.pdf", sep = "."), width = 6, height = 8)
pheatmap(CountMat, show_rownames = F, scale = "row") #, color = colorRampPalette(c("blue", "white", "red"))(11))
dev.off()

datastats[is.na(Gene),CpGquantile := "Enh"]
# datastats[!is.na(Gene),CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25)), 
#                           labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
datastats[!is.na(Gene),CpGquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25)), 
                                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
CountMat <- table(datastats[CommunitySize >= 3, c(2,13)])
ppcount <- apply(CountMat[,1:4], 1, sum)
counts_filtered_df <- CountMat[apply(CountMat, MARGIN = 1, FUN = function(x) sd(x) != 0),]
pdf(file = paste(args[1], "pheatmap.CommunityCpG.pdf", sep = "."), width = 6, height = 8)
pheatmap(counts_filtered_df, show_rownames = F, scale = "row") #, color = colorRampPalette(c("blue", "white", "red"))(11))
# pheatmap(CountMat[,1:4], show_rownames = F, scale = "row")
# pheatmap(CountMat[ppcount > 1,1:4], show_rownames = F, scale = "row", cluster_cols = T)
dev.off()
combination <- paste(round(CountMat[,1]/ppcount,1), round(CountMat[,2]/ppcount,1), 
                     round(CountMat[,3]/ppcount,1), round(CountMat[,4]/ppcount,1), sep = "|")
sort(table(combination), decreasing = T)[1:10]

datastats[is.na(Gene),Exprquantile := "Enh"]
datastats[!is.na(EBV.medianexpression),Exprquantile := cut(EBV.medianexpression, 
                                           breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T), 
                                           labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
CountMat <- table(datastats[CommunitySize >= 3, c(2,14)])
counts_filtered_df <- CountMat[apply(CountMat, MARGIN = 1, FUN = function(x) sd(x) != 0),]
pdf(file = paste(args[1], "pheatmap.CommunityExpr.pdf", sep = "."), width = 6, height = 8)
pheatmap(counts_filtered_df, show_rownames = F, scale = "row") 
#, clustering_distance_rows = "euclidean", clustering_method = "mcquitty") #, color = colorRampPalette(c("blue", "white", "red"))(11))
# heatmap(CountMat, scale = "row", labRow = NA, ke)
dev.off()

CommunityEdges <- data.frame(Community = NA, CommunitySize = NA, CommunityAllEdges = NA, 
                             ThinCount = NA, ThickCount = NA, ShortCount = NA, LongCount = NA)
k <- 1
for(i in 1:length(dg)){
    a <- as.data.frame(edge.attributes(graph = dg[[i]]))
    size <- length(a)
    for(j in 1:size){
        node2community[k, ] <- c(a[j], paste("comp", i, sep = "."), size)
        k <- k + 1
    }
}
# as.data.frame(edge.attributes(graph = dg[[1]]))
# table(edge.attributes(graph = dg[[1]])$Thickness)

##################
# for ggplot2
datax <- output[!is.na(output$Gene),]
datax$Type <- factor(datax$Type, levels = c("HKG.1", "HKG.2", "other", "EBV.TSG", "other.TSG"))
datax[is.na(datax$Type),11] <- "other"

datay <- data.table(CommunitySize = as.numeric(datax$CommunitySize), Type = datax$Type)
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "other", "EBV.TSG", "other.TSG"))
myplot <- ggplot(datay, aes(x = Type, y = CommunitySize, fill = Type)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "CommunitySize.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()  
    
datay <- data.table(CommunitySize = as.numeric(datax$CommunitySize), 
                    CpGRatio = datax$CpG.Ratio)
                    # CpGRatio = as.numeric(as.character(datax$CpG.Ratio)))
datay[,CpGquantile := cut(CpGRatio, breaks = quantile(CpGRatio, probs = seq(0, 1, 0.25)), 
                           labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
myplot <- ggplot(datay, aes(x = CpGquantile, y = CommunitySize, fill = CpGquantile)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "CommunitySize.by.CpGRatio.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

datay <- data.table(CommunitySize = as.numeric(datax$CommunitySize), 
                    CpGCount = datax$CpG.Count)
datay[,CpGquantile := cut(CpGCount, breaks = quantile(CpGCount, probs = seq(0, 1, 0.25)), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
myplot <- ggplot(datay, aes(x = CpGquantile, y = CommunitySize, fill = CpGquantile)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "CommunitySize.by.CpGCount.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

datay <- data.table(NodeCloseness = as.numeric(datax$NodeCloseness), 
                    CpGCount = datax$CpG.Count)
datay[,CpGquantile := cut(CpGCount, breaks = quantile(CpGCount, probs = seq(0, 1, 0.25)), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
myplot <- ggplot(datay, aes(x = CpGquantile, y = NodeCloseness, fill = CpGquantile)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "NodeCloseness.by.CpGCount.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

datay <- data.table(CommunitySize = as.numeric(datax$CommunitySize), 
                    Expr = datax$EBV.medianexpression)
datay[,Exprquantile := cut(Expr, breaks = quantile(Expr, probs = seq(0, 1, 0.25), na.rm = T), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
datay[,ComSizequantile := cut(CommunitySize, breaks = quantile(unique(CommunitySize), probs = seq(0, 1, 0.25)), 
                           labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
myplot <- ggplot(datay, aes(x = Exprquantile, y = CommunitySize, fill = Exprquantile)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "CommunitySize.by.Expr.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()
myplot <- ggplot(datay, aes(x = ComSizequantile, y = Expr, fill = ComSizequantile)) +
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "Expr.by.CommunitySize.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()

datay <- data.table(NodeBetweenness = as.numeric(datax$NodeBetweenness), 
                    Expression = datax$EBV.medianexpression)
myplot <- ggplot(datay, aes(x = log2(NodeBetweenness + 1), y = Expression)) + geom_point()
pdf(file = paste(args[1], "NodeBetweenness.EBVmedianexpression.pdf", sep = "."), width = 8, height = 8)
print(myplot)  
dev.off()


# median/mean CpG for a community??????
datay <- data.table(Community = datax$Community, 
                    CommunitySize = as.numeric(datax$CommunitySize), 
                    CpGCount = as.numeric(datax$CpG.Count),
                    CpGRatio = as.numeric(datax$CpG.Ratio),
                    Expression = datax$EBV.medianexpression)
datay[, MedianCpG:=median(CpGCount, na.rm = T), by = Community]
datay[, MedianCpGR:=median(CpGRatio, na.rm = T), by = Community]
datay[, MedianExpr:=median(Expression, na.rm = T), by = Community]
dataz <- unique(datay[,c(1,2,6)])
myplot <- ggplot(dataz[CommunitySize > 5], aes(x = log2(CommunitySize), y = MedianCpG)) + geom_point()
pdf(file = paste(args[1], "CommunitySize.medianCpGCount.pdf", sep = "."), width = 8, height = 8)
print(myplot)  
dev.off()
dataz <- unique(datay[,c(1,2,7)])
myplot <- ggplot(dataz[CommunitySize > 5], aes(x = log2(CommunitySize), y = MedianCpGR)) + geom_point()
pdf(file = paste(args[1], "CommunitySize.medianCpGRatio.pdf", sep = "."), width = 8, height = 8)
print(myplot)  
dev.off()
dataz <- unique(datay[,c(1,2,8)])
myplot <- ggplot(dataz[CommunitySize > 5], aes(x = log2(CommunitySize), y = MedianExpr)) + geom_point()
pdf(file = paste(args[1], "CommunitySize.MedianExpr.pdf", sep = "."), width = 8, height = 8)
print(myplot)  
dev.off()

# 
# cor(as.numeric(datax$CommunitySize), datax$CpG.Count, method = "spearman")
# cor(as.numeric(datax$CommunitySize), datax$CpG.Ratio, method = "spearman")
# cor(as.numeric(datax$CommunitySize), datax$Length, method = "spearman")
# cor(datax$NodeDegree, datax$CpG.Count, method = "spearman")
# cor(datax$NodeDegree, datax$CpG.Ratio, method = "spearman")
# cor(datax$NodeDegree, datax$Length, method = "spearman")
# cor(datax$NodeCloseness, datax$CpG.Count, method = "spearman")
# cor(datax$NodeCloseness, datax$CpG.Ratio, method = "spearman")
# cor(datax$NodeCloseness, datax$Length, method = "spearman")
# cor(datax$NodeBetweenness, datax$CpG.Count, method = "spearman")
# cor(datax$NodeBetweenness, datax$CpG.Ratio, method = "spearman")
# cor(datax$NodeBetweenness, datax$Length, method = "spearman")
# cor(datax$NodeBetweenness, datax$EBV.medianexpression, method = "spearman")

# cliques: complete subgraph
# cls <- clusters(network)
# g <- sample_gnp(20, 0.1)
# plot(g)
# clique_num(g)
# cliques(g, min = 3)
# largest_cliques(g)
# max_cliques(g)

# cliq <- cliques(network, min = 5)
# plot(dg[[1]])
# bigclique <- largest_cliques(dg[[2]])
# maxclique <- max_cliques(dg[[2]])
# cliq <- cliques(dg[[2]], min = 3)

# node2clique <- data.frame(Clique = NA, Node = NA, Length = NA, 
#                           CpG.Count = NA, CpG.Ratio = NA, 
#                           Type = NA, EBV.medianexpression = NA)
node2clique <- data.frame(Clique = NA, CliqueSize = NA, Node = NA)
p <- 0
for(i in 1:length(dg)){
    # cliq <- cliques(dg[[i]], min = 3)
    cliq <- largest_cliques(dg[[i]])
    if(length(cliq) > 0 && length(cliq) < 1000){
        print(paste(i, "/", length(dg), ": ",  length(cliq), sep = ""))
        for(j in 1:length(cliq)){
            vx <- names(cliq[[j]])
            if(length(vx[grep("ENSG", vx)]) >= 3){
                for(k in 1:length(vx)){
                    p <- p + 1
                    cname <- paste("component", i, "clique", j, sep = ".")
                    node2clique[p,] <- c(cname, length(vx), vx[k])
                    # node2clique[p,3:7] <- nodes[nodes[,1] %in% vx[k],2:6]
                }
            }
        }
    }
}
# node2clique
# x <- nodes[match(node2clique$Node, nodes[,1]),]
# dim(x)
# dim(node2clique)
node2clique.out <- data.frame(node2clique, nodes[match(node2clique$Node, nodes[,1]),])
write.table(node2clique.out, file = paste(args[1], "igraph.node2clique.out.tsv", sep = "."), 
            quote = F, sep = "\t", row.names = F, col.names = T)


##################
# for ggplot2
r3.rown <- sample(1:nrow(node2clique.out), size = 3000, replace = TRUE)
rands1 <- data.table(Clique = paste("r3", rep(1:1000, 3), sep = "."), 
                    CliqueSize = rep("r3", 3000), 
                    Gene = node2clique.out$Gene[r3.rown],
                    CpGCount = as.numeric(node2clique.out$CpG.Count[r3.rown]),
                    CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r3.rown]),
                    Expression = node2clique.out$EBV.medianexpression[r3.rown])
r4.rown <- sample(1:nrow(node2clique.out), size = 4000, replace = TRUE)
rands2 <- data.table(Clique = paste("r4", rep(1:1000, 4), sep = "."), 
                     CliqueSize = rep("r4", 4000), 
                     Gene = node2clique.out$Gene[r4.rown],
                     CpGCount = as.numeric(node2clique.out$CpG.Count[r4.rown]),
                     CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r4.rown]),
                     Expression = node2clique.out$EBV.medianexpression[r4.rown])
r5.rown <- sample(1:nrow(node2clique.out), size = 5000, replace = TRUE)
rands3 <- data.table(Clique = paste("r5", rep(1:1000, 5), sep = "."), 
                     CliqueSize = rep("r5", 5000), 
                     Gene = node2clique.out$Gene[r5.rown],
                     CpGCount = as.numeric(node2clique.out$CpG.Count[r5.rown]),
                     CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r5.rown]),
                     Expression = node2clique.out$EBV.medianexpression[r5.rown])
r6.rown <- sample(1:nrow(node2clique.out), size = 6000, replace = TRUE)
rands4 <- data.table(Clique = paste("r6", rep(1:1000, 6), sep = "."), 
                     CliqueSize = rep("r6", 6000), 
                     Gene = node2clique.out$Gene[r6.rown],
                     CpGCount = as.numeric(node2clique.out$CpG.Count[r6.rown]),
                     CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r6.rown]),
                     Expression = node2clique.out$EBV.medianexpression[r6.rown])
r7.rown <- sample(1:nrow(node2clique.out), size = 7000, replace = TRUE)
rands5 <- data.table(Clique = paste("r7", rep(1:1000, 7), sep = "."), 
                     CliqueSize = rep("r7", 7000), 
                     Gene = node2clique.out$Gene[r7.rown],
                     CpGCount = as.numeric(node2clique.out$CpG.Count[r7.rown]),
                     CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r7.rown]),
                     Expression = node2clique.out$EBV.medianexpression[r7.rown])
r8.rown <- sample(1:nrow(node2clique.out), size = 8000, replace = TRUE)
rands6 <- data.table(Clique = paste("r8", rep(1:1000, 8), sep = "."), 
                     CliqueSize = rep("r8", 8000),
                     Gene = node2clique.out$Gene[r8.rown],
                     CpGCount = as.numeric(node2clique.out$CpG.Count[r8.rown]),
                     CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r8.rown]),
                     Expression = node2clique.out$EBV.medianexpression[r8.rown])
rands <- rbindlist(list(rands1, rands2, rands3, rands4, rands5, rands6))
rands[, MedianCpG:=median(CpGCount, na.rm = T), by = Clique]
rands[, MedianCpGR:=median(CpGRatio, na.rm = T), by = Clique]
rands[, MedianExpr:=median(Expression, na.rm = T), by = Clique]

datax <- node2clique.out[!is.na(node2clique.out$Gene),]
datax$Type <- factor(datax$Type, levels = c("HKG.1", "HKG.2", "other", "EBV.TSG", "other.TSG"))
datax[is.na(datax$Type),8] <- "other"

TypeCount.3 <- table(datax[datax$CliqueSize==3, c(1,8)])
pdf(file = paste(args[1], "pheatmap.CliqueSize.3.pdf", sep = "."), width = 6, height = 8)
pheatmap(TypeCount.3, show_rownames = F, color = colorRampPalette(c("white", "red"))(3))
dev.off()
TypeCount.4 <- table(datax[datax$CliqueSize==4, c(1,8)])
pdf(file = paste(args[1], "pheatmap.CliqueSize.4.pdf", sep = "."), width = 6, height = 8)
pheatmap(TypeCount.4, show_rownames = F, color = colorRampPalette(c("white", "red"))(4))
dev.off()
TypeCount.5 <- table(datax[datax$CliqueSize==5, c(1,8)])
pdf(file = paste(args[1], "pheatmap.CliqueSize.5.pdf", sep = "."), width = 6, height = 8)
pheatmap(TypeCount.5, show_rownames = F, color = colorRampPalette(c("white", "red"))(5))
dev.off()

datay <- data.table(Clique = datax$Clique, CliqueSize = datax$CliqueSize, CpGCount = datax$CpG.Count)
datay[,CpGquantile := cut(CpGCount, breaks = quantile(CpGCount, probs = seq(0, 1, 0.25)), 
                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
CpGCount.count <- table(datay[CliqueSize==4, c(1,4)])
pdf(file = paste(args[1], "pheatmap.CliqueSize.CpGCount.pdf", sep = "."), width = 6, height = 8)
pheatmap(CpGCount.count, show_rownames = F, scale = "none", color = colorRampPalette(c("white", "red"))(6))
dev.off()


datay <- data.table(CliqueSize = as.numeric(datax$CliqueSize), Type = datax$Type)
myplot <- ggplot(datay, aes(x = Type, y = CliqueSize, fill = Type)) + 
    geom_boxplot(width = 0.8)
pdf(file = paste(args[1], "CliqueSize.by.Type.pdf", sep = "."), width = 6, height = 8)
print(myplot)  
dev.off()  

datay <- data.table(Clique = datax$Clique, 
                    CliqueSize = as.numeric(datax$CliqueSize), 
                    Gene = datax$Gene,
                    CpGCount = as.numeric(datax$CpG.Count),
                    CpGRatio = as.numeric(datax$CpG.Ratio),
                    Expression = datax$EBV.medianexpression)
datay[, MedianCpG:=median(CpGCount, na.rm = T), by = Clique]
datay[, MedianCpGR:=median(CpGRatio, na.rm = T), by = Clique]
datay[, MedianExpr:=median(Expression, na.rm = T), by = Clique]
total <- data.table(Clique = paste("All", 1:nrow(nodes), sep = "."), 
                    CliqueSize = rep("Total57820", nrow(nodes)), 
                    Gene = nodes$Gene, 
                    CpGCount = as.numeric(nodes$CpG.Count),
                    CpGRatio = as.numeric(nodes$CpG.Ratio),
                    Expression = nodes$EBV.medianexpression,
                    MedianCpG = as.numeric(nodes$CpG.Count),
                    MedianCpGR = as.numeric(nodes$CpG.Ratio),
                    MedianExpr = nodes$EBV.medianexpression)

datayy <- rbindlist(list(datay, rands))
dataz <- unique(datayy[,c(3,2,4)])
myplot <- ggplot(dataz, aes(x = factor(CliqueSize), y = CpGCount, fill = factor(CliqueSize))) + geom_boxplot()
pdf(file = paste(args[1], "CliqueSize.CpGCount.pdf", sep = "."), width = 8, height = 8)
print(myplot)  
dev.off()
dataz <- unique(datayy[,c(3,2,5)])
myplot <- ggplot(dataz, aes(x = factor(CliqueSize), y = CpGRatio, fill = factor(CliqueSize))) + geom_boxplot()
pdf(file = paste(args[1], "CliqueSize.CpGRatio.pdf", sep = "."), width = 8, height = 8)
print(myplot)  
dev.off()
dataz <- unique(datayy[,c(3,2,6)])
myplot <- ggplot(dataz, aes(x = factor(CliqueSize), y = Expression, fill = factor(CliqueSize))) + geom_boxplot()
pdf(file = paste(args[1], "CliqueSize.Expression.pdf", sep = "."), width = 8, height = 8)
print(myplot)  
dev.off()

cliquestats <- as.data.table(node2clique.out)
cliquestats[is.na(Gene), Type:="Enh"]
cliquestats[is.na(Type), Type:="other"]
CountMat <- table(cliquestats[CliqueSize >= 3, c(1,8)])
pdf(file = paste(args[1], "pheatmap.Clique.pdf", sep = "."), width = 6, height = 8)
pheatmap(CountMat, show_rownames = F, scale = "row") #, color = colorRampPalette(c("blue", "white", "red"))(11))
dev.off()

cliquestats[is.na(Gene),CpGquantile := "Enh"]
# datastats[!is.na(Gene),CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25)), 
#                           labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
cliquestats[!is.na(Gene),CpGquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25)), 
                                          labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
CountMat <- table(cliquestats[CliqueSize >= 3, c(1,10)])
ppcount <- apply(CountMat[,-1], 1, sum)
counts_filtered_df <- CountMat[apply(CountMat, MARGIN = 1, FUN = function(x) sd(x) != 0),]
pdf(file = paste(args[1], "pheatmap.CliqueCpG.pdf", sep = "."), width = 6, height = 8)
pheatmap(counts_filtered_df, show_rownames = F, scale = "row") #, color = colorRampPalette(c("blue", "white", "red"))(11))
# pheatmap(CountMat[,1:4], show_rownames = F, scale = "row")
# pheatmap(CountMat[ppcount > 1,1:4], show_rownames = F, scale = "row", cluster_cols = T)
dev.off()
combination <- paste(round(CountMat[,2]/ppcount,1), round(CountMat[,3]/ppcount,1), 
                     round(CountMat[,4]/ppcount,1), round(CountMat[,5]/ppcount,1), sep = "|")
sort(table(combination), decreasing = T)[1:10]

cliquestats[is.na(Gene),Exprquantile := "Enh"]
cliquestats[!is.na(EBV.medianexpression),Exprquantile := cut(EBV.medianexpression, 
                                                           breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T), 
                                                           labels = paste("Q", 1:4, sep = ""), include.lowest = TRUE)]
CountMat <- table(cliquestats[CliqueSize >= 3, c(1,11)])
counts_filtered_df <- CountMat[apply(CountMat, MARGIN = 1, FUN = function(x) sd(x) != 0),]
pdf(file = paste(args[1], "pheatmap.CliqueExpr.pdf", sep = "."), width = 6, height = 8)
pheatmap(counts_filtered_df, show_rownames = F, scale = "row") 
#, clustering_distance_rows = "euclidean", clustering_method = "mcquitty") #, color = colorRampPalette(c("blue", "white", "red"))(11))
# heatmap(CountMat, scale = "row", labRow = NA, ke)
dev.off()

# r3.rown <- sample(1:nrow(node2clique.out), size = 3000, replace = TRUE)
# rands1 <- data.table(Clique = paste("r3", rep(1:1000, 3), sep = "."), 
#                      CliqueSize = rep("r3", 3000), 
#                      CpGCount = as.numeric(node2clique.out$CpG.Count[r3.rown]),
#                      CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r3.rown]),
#                      Expression = node2clique.out$EBV.medianexpression[r3.rown])
# r4.rown <- sample(1:nrow(node2clique.out), size = 4000, replace = TRUE)
# rands2 <- data.table(Clique = paste("r4", rep(1:1000, 4), sep = "."), 
#                      CliqueSize = rep("r4", 4000), 
#                      CpGCount = as.numeric(node2clique.out$CpG.Count[r4.rown]),
#                      CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r4.rown]),
#                      Expression = node2clique.out$EBV.medianexpression[r4.rown])
# r5.rown <- sample(1:nrow(node2clique.out), size = 5000, replace = TRUE)
# rands3 <- data.table(Clique = paste("r5", rep(1:1000, 5), sep = "."), 
#                      CliqueSize = rep("r5", 5000), 
#                      CpGCount = as.numeric(node2clique.out$CpG.Count[r5.rown]),
#                      CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r5.rown]),
#                      Expression = node2clique.out$EBV.medianexpression[r5.rown])
# r6.rown <- sample(1:nrow(node2clique.out), size = 6000, replace = TRUE)
# rands4 <- data.table(Clique = paste("r6", rep(1:1000, 6), sep = "."), 
#                      CliqueSize = rep("r6", 6000), 
#                      CpGCount = as.numeric(node2clique.out$CpG.Count[r6.rown]),
#                      CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r6.rown]),
#                      Expression = node2clique.out$EBV.medianexpression[r6.rown])
# r7.rown <- sample(1:nrow(node2clique.out), size = 7000, replace = TRUE)
# rands5 <- data.table(Clique = paste("r7", rep(1:1000, 7), sep = "."), 
#                      CliqueSize = rep("r7", 7000), 
#                      CpGCount = as.numeric(node2clique.out$CpG.Count[r7.rown]),
#                      CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r7.rown]),
#                      Expression = node2clique.out$EBV.medianexpression[r7.rown])
# r8.rown <- sample(1:nrow(node2clique.out), size = 8000, replace = TRUE)
# rands6 <- data.table(Clique = paste("r8", rep(1:1000, 8), sep = "."), 
#                      CliqueSize = rep("r8", 8000), 
#                      CpGCount = as.numeric(node2clique.out$CpG.Count[r8.rown]),
#                      CpGRatio = as.numeric(node2clique.out$CpG.Ratio[r8.rown]),
#                      Expression = node2clique.out$EBV.medianexpression[r8.rown])
# rands <- rbindlist(list(rands1, rands2, rands3, rands4, rands5, rands6))
# rands[, MedianCpG:=median(CpGCount, na.rm = T), by = Clique]
# rands[, MedianCpGR:=median(CpGRatio, na.rm = T), by = Clique]
# rands[, MedianExpr:=median(Expression, na.rm = T), by = Clique]
# 
# datax <- node2clique.out[!is.na(node2clique.out$Gene),]
# datax$Type <- factor(datax$Type, levels = c("HKG.1", "HKG.2", "other", "EBV.TSG", "other.TSG"))
# datax[is.na(datax$Type),8] <- "other"
# 
# datay <- data.table(CliqueSize = as.numeric(datax$CliqueSize), Type = datax$Type)
# myplot <- ggplot(datay, aes(x = Type, y = CliqueSize, fill = Type)) + 
#     geom_boxplot(width = 0.8)
# pdf(file = paste(args[1], "CliqueSize.by.Type.pdf", sep = "."), width = 6, height = 8)
# print(myplot)  
# dev.off()  
# 
# datay <- data.table(Clique = datax$Clique, 
#                     CliqueSize = as.numeric(datax$CliqueSize), 
#                     CpGCount = as.numeric(datax$CpG.Count),
#                     CpGRatio = as.numeric(datax$CpG.Ratio),
#                     Expression = datax$EBV.medianexpression)
# datay[, MedianCpG:=median(CpGCount, na.rm = T), by = Clique]
# datay[, MedianCpGR:=median(CpGRatio, na.rm = T), by = Clique]
# datay[, MedianExpr:=median(Expression, na.rm = T), by = Clique]
# total <- data.table(Clique = paste("All", 1:nrow(nodes), sep = "."), 
#                     CliqueSize = rep("Total57820", nrow(nodes)), 
#                     CpGCount = as.numeric(nodes$CpG.Count),
#                     CpGRatio = as.numeric(nodes$CpG.Ratio),
#                     Expression = nodes$EBV.medianexpression,
#                     MedianCpG = as.numeric(nodes$CpG.Count),
#                     MedianCpGR = as.numeric(nodes$CpG.Ratio),
#                     MedianExpr = nodes$EBV.medianexpression)
# 
# datayy <- rbindlist(list(datay, rands))
# dataz <- unique(datayy[,c(1,2,6)])
# myplot <- ggplot(dataz, aes(x = factor(CliqueSize), y = MedianCpG, fill = factor(CliqueSize))) + geom_boxplot()
# pdf(file = paste(args[1], "CliqueSize.medianCpGCount.pdf", sep = "."), width = 8, height = 8)
# print(myplot)  
# dev.off()
# dataz <- unique(datayy[,c(1,2,7)])
# myplot <- ggplot(dataz, aes(x = factor(CliqueSize), y = MedianCpGR, fill = factor(CliqueSize))) + geom_boxplot()
# pdf(file = paste(args[1], "CliqueSize.medianCpGRatio.pdf", sep = "."), width = 8, height = 8)
# print(myplot)  
# dev.off()
# dataz <- unique(datayy[,c(1,2,8)])
# myplot <- ggplot(dataz, aes(x = factor(CliqueSize), y = MedianExpr, fill = factor(CliqueSize))) + geom_boxplot()
# pdf(file = paste(args[1], "CliqueSize.MedianExpr.pdf", sep = "."), width = 8, height = 8)
# print(myplot)  
# dev.off()
