library(igraph)
library(ggplot2)
library(reshape2)
library(data.table)
library(pheatmap)
library(corrplot)
args <- commandArgs(TRUE)
# add CpG island lenth for node feature
# setwd("~/Data/myworkData/igraph")
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

# basic gene features ------------------------------------------------------------
gene.feature <- data.table(nodes, node.feature[match(nodes$Gene, node.feature$node.name),])
gene.feature[is.na(Type), Type:="Other"]

# community features ------------------------------------------------------------
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

# node edge attritubes ------------------------------------------------------------
# get all the edges connected to the specific nodes
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
node2edges.feature <- as.data.table(node2edges)
node2edges.feature[,Thick.percentage:=as.numeric(ThickCount)/(as.numeric(ThickCount)+as.numeric(ThinCount))]
node2edges.feature[,Long.percentage:=as.numeric(LongCount)/(as.numeric(LongCount)+as.numeric(ShortCount))]
node2edges.feature[,PP.percentage:=as.numeric(PPCount)/(as.numeric(PPCount)+as.numeric(PECount))]

# Fature summary ------------------------------------------------------------
node.feature.out <- data.table(loopnodes.feature, 
                               node2community[match(loopnodes.feature$node.name, node2community[,1]),],
                               node2edges.feature[match(loopnodes.feature$node.name, node2edges.feature$Node),])
node.feature.out[,Lengthquantile := cut(Length, breaks = quantile(Length, probs = seq(0, 1, 0.25), na.rm = T), 
                                        labels = paste("LenQ", 1:4, sep = ""), include.lowest = TRUE)]
node.feature.out[,CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25), na.rm = T), 
                                     labels = paste("CpGQ", 1:4, sep = ""), include.lowest = TRUE)]
node.feature.out[,CpGRquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25), na.rm = T), 
                                      labels = paste("CpGRQ", 1:4, sep = ""), include.lowest = TRUE)]
node.feature.out[EBV.medianexpression > 0,
                 Exprquantile := cut(EBV.medianexpression, 
                                     breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T),
                                     labels = paste("ExprQ", 1:4, sep = ""), include.lowest = TRUE)]
node.feature.out[,CGIquantile := cut(cpgIslandLenth, breaks = quantile(cpgIslandLenth, probs = seq(0, 1, 0.25), na.rm = T), 
                                     labels = paste("cgiQ", 1:4, sep = ""), include.lowest = TRUE)]
node.feature.out[EBV.medianexpression == 0, Exprquantile := "NULL"]
node.feature.out[is.na(EBV.medianexpression), Exprquantile := "NULL"]
node.feature.out[is.na(Lengthquantile), Lengthquantile := "NULL"]
node.feature.out[is.na(CpGquantile), CpGquantile := "NULL"]
node.feature.out[is.na(CpGRquantile), CpGRquantile := "NULL"]
node.feature.out[is.na(CGIquantile), CGIquantile := "NULL"]
node.feature.out[, Communityquantile:= cut(as.numeric(CommunitySize), 
                                           breaks = quantile(as.numeric(CommunitySize), probs = seq(0, 1, 0.25), na.rm = T), 
                                           labels = paste("commQ", 1:4, sep = ""), include.lowest = TRUE)]
node.feature.out[,PPFlag := cut(as.numeric(PPCount), 
                                  breaks = quantile(as.numeric(node.feature.out[Type!="Enh"]$PPCount), probs = seq(0, 1, 0.25), na.rm = T), 
                                  labels = paste("PPQ", 1:4, sep = ""), include.lowest = TRUE)]
node.feature.out[is.na(Communityquantile), Communityquantile := "Lonely"]
node.feature.out[is.na(PPFlag), PPFlag := "NULL"]
node.feature.out[is.na(Gene), Lengthquantile:="Enh"]
node.feature.out[is.na(Gene), CpGquantile:="Enh"]
node.feature.out[is.na(Gene), CpGRquantile:="Enh"]
node.feature.out[is.na(Gene), Exprquantile:="Enh"]
node.feature.out[is.na(Gene), CGIquantile:="Enh"]

gene.feature.out <- data.table(gene.feature, 
                               node2community[match(gene.feature$Gene, node2community[,1]),],
                               node2edges.feature[match(gene.feature$Gene, node2edges.feature$Node),])
gene.feature.out[,Lengthquantile := cut(Length, breaks = quantile(Length, probs = seq(0, 1, 0.25), na.rm = T), 
                                       labels = paste("LenQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature.out[,CpGquantile := cut(CpG.Count, breaks = quantile(CpG.Count, probs = seq(0, 1, 0.25), na.rm = T), 
                                    labels = paste("CpGQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature.out[,CpGRquantile := cut(CpG.Ratio, breaks = quantile(CpG.Ratio, probs = seq(0, 1, 0.25), na.rm = T), 
                                     labels = paste("CpGRQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature.out[EBV.medianexpression > 0,
                      Exprquantile := cut(EBV.medianexpression, 
                                          breaks = quantile(EBV.medianexpression, probs = seq(0, 1, 0.25), na.rm = T),
                                          labels = paste("ExprQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature.out[,CGIquantile := cut(cpgIslandLenth, breaks = quantile(cpgIslandLenth, probs = seq(0, 1, 0.25), na.rm = T), 
                                    labels = paste("cgiQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature.out[EBV.medianexpression == 0, Exprquantile := "NULL"]
gene.feature.out[is.na(EBV.medianexpression), Exprquantile := "NULL"]
gene.feature.out[is.na(Lengthquantile), Lengthquantile := "NULL"]
gene.feature.out[is.na(CpGquantile), CpGquantile := "NULL"]
gene.feature.out[is.na(CpGRquantile), CpGRquantile := "NULL"]
gene.feature.out[is.na(CGIquantile), CGIquantile := "NULL"]
gene.feature.out[, Communityquantile:= cut(as.numeric(CommunitySize), 
                                            breaks = quantile(as.numeric(CommunitySize), probs = seq(0, 1, 0.25), na.rm = T), 
                                            labels = paste("commQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature.out[,PPFlag := cut(as.numeric(gene.feature.out$PPCount), 
                                       breaks = quantile(as.numeric(gene.feature.out$PPCount), 
                                                         probs = seq(0, 1, 0.25), na.rm = T), 
                                       labels = paste("PPQ", 1:4, sep = ""), include.lowest = TRUE)]
gene.feature.out[is.na(Communityquantile), Communityquantile := "Lonely"]
gene.feature.out[is.na(PPFlag), PPFlag := "NULL"]

# node clique features ------------------------------------------------------------
# because one node could belong to many different cliquens, 
# so we cannot assign node to unique clique
node2clique <- data.frame(Clique = NA, CliqueSize = NA, Node = NA)
p <- 0
for(i in 1:length(dg)){
    cliq <- largest_cliques(dg[[i]])
    if(length(cliq) > 0 && length(cliq) < 1000){
        print(paste(i, "/", length(dg), ": ",  length(cliq), sep = ""))
        for(j in 1:length(cliq)){
            vx <- names(cliq[[j]])
            if(length(vx[grep("ENSG", vx)]) >= 2){
                for(k in 1:length(vx)){
                    p <- p + 1
                    cname <- paste("component", i, "clique", j, sep = ".")
                    node2clique[p,] <- c(cname, length(vx), vx[k])
                }
            }
        }
    }
}

node2clique.out <- data.frame(node2clique, node.feature.out[match(node2clique$Node, node.feature.out$node.name),])
write.table(node2clique.out, file = paste(args[1], "igraph.node2clique.out.tsv", sep = "."), 
            quote = F, sep = "\t", row.names = F, col.names = T)

# save data ------------------------------------------------------------
write.table(node.feature.out, file = paste(args[1], "igraph.node.feature.out.tsv", sep = "."), 
            quote = F, sep = "\t", row.names = F, col.names = T)
write.table(gene.feature.out, file = paste(args[1], "igraph.gene.feature.out.tsv", sep = "."), 
            quote = F, sep = "\t", row.names = F, col.names = T)

# boxplot ------------------------------------------------------------
myboxplot <- function(datay = datay, ylabs = ylab, pre = args[1]){
    require(ggplot2)
    colnames(datay) <- c("Class", "value")
    ymin <- quantile(as.numeric(datay$value), probs = 0.01, na.rm = T)
    ymax <- quantile(as.numeric(datay$value), probs = 0.99, na.rm = T)
    myplot <- ggplot(datay, aes(x = Class, y = as.numeric(value), fill = Class)) + 
        geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
        labs(y = ylabs)
    pdf(file = paste("batchboxplot", pre, "pdf", sep = "."), width = 6, height = 8)
    print(myplot)
    dev.off()
}
node.feature.out$Type <- factor(node.feature.out$Type, 
                                levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
for(i in c(2:3,6:8,10:11,14,16:26)){
    for(j in c(9,27:33)){
        mydata <- data.table(Class = node.feature.out[,.SD,.SDcols=j], 
                             value = node.feature.out[,.SD,.SDcols=i])
        myboxplot(datay = mydata, ylabs = colnames(node.feature.out)[i],
                  pre = paste(args[1], "node", 
                              colnames(node.feature.out)[i], "by",
                              colnames(node.feature.out)[j], sep = "."))
    }
}
for(i in c(4)){
    for(j in c(9,27:33)){
        mydata <- data.table(Class = node.feature.out[,.SD,.SDcols=j], 
                             value = log2(node.feature.out[,.SD,.SDcols=i] + 1))
        myboxplot(datay = mydata, ylabs = colnames(node.feature.out)[i],
                  pre = paste(args[1], "node", 
                              colnames(node.feature.out)[i], "by",
                              colnames(node.feature.out)[j], sep = "."))
    }
}

gene.feature.out$Type <- factor(gene.feature.out$Type, 
                                levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
for(i in c(2:4,6:7,9:10,14,16:26)){
    for(j in c(5,27:33)){
        mydata <- data.table(Class = gene.feature.out[,.SD,.SDcols=j], 
                             value = gene.feature.out[,.SD,.SDcols=i])
        myboxplot(datay = mydata, ylabs = colnames(gene.feature.out)[i],
                  pre = paste(args[1], "gene", 
                              colnames(gene.feature.out)[i], "by",
                              colnames(gene.feature.out)[j], sep = "."))
    }
}
for(i in c(11)){
    for(j in c(5,27:33)){
        mydata <- data.table(Class = gene.feature.out[,.SD,.SDcols=j], 
                             value = log2(gene.feature.out[,.SD,.SDcols=i] + 1))
        myboxplot(datay = mydata, ylabs = colnames(gene.feature.out)[i],
                  pre = paste(args[1], "gene", 
                              colnames(gene.feature.out)[i], "by",
                              colnames(gene.feature.out)[j], sep = "."))
    }
}

for(i in c(2:4,6:7,9:10,14,16:26)){
    for(j in c(5,27:33)){
        mydata <- data.table(Class = gene.feature.out[,.SD,.SDcols=j], 
                             value = gene.feature.out[,.SD,.SDcols=i])
        colnames(mydata) <- c("Class", "value")
        mydata[is.na(value), value:=0]
        myboxplot(datay = mydata, ylabs = colnames(gene.feature.out)[i],
                  pre = paste(args[1], "NAas0.gene", 
                              colnames(gene.feature.out)[i], "by",
                              colnames(gene.feature.out)[j], sep = "."))
    }
}
for(i in c(11)){
    for(j in c(5,27:33)){
        mydata <- data.table(Class = gene.feature.out[,.SD,.SDcols=j], 
                             value = log2(gene.feature.out[,.SD,.SDcols=i] + 1))
        colnames(mydata) <- c("Class", "value")
        mydata[is.na(value), value:=0]
        myboxplot(datay = mydata, ylabs = colnames(gene.feature.out)[i],
                  pre = paste(args[1], "NAas0.gene", 
                              colnames(gene.feature.out)[i], "by",
                              colnames(gene.feature.out)[j], sep = "."))
    }
}

# scatter plot ------------------------------------------------------------
myscatterplot <- function(datay = datay, xlabs = xlabs, ylabs = ylabs, pre = args[1]){
    require(ggplot2)
    colnames(datay) <- c("SeqFeature", "GraphFeature", "Group")
    xmin <- quantile(as.numeric(datay$SeqFeature), probs = 0.01, na.rm = T)
    xmax <- quantile(as.numeric(datay$SeqFeature), probs = 0.99, na.rm = T)
    ymin <- quantile(as.numeric(datay$GraphFeature), probs = 0.01, na.rm = T)
    ymax <- quantile(as.numeric(datay$GraphFeature), probs = 0.99, na.rm = T)
    myplot <- ggplot(datay, aes(x = as.numeric(datay$SeqFeature), 
                                y = as.numeric(GraphFeature), colour = Group)) + 
        geom_point(size = 0.8) + scale_x_continuous(limits = c(xmin, xmax)) + 
        scale_y_continuous(limits = c(ymin, ymax)) + 
        labs(x = xlabs, y = ylabs)
    pdf(file = paste("batchscatterplot", pre, "pdf", sep = "."), width = 9, height = 8)
    print(myplot)
    dev.off()
}
node.feature.out$Type <- factor(node.feature.out$Type, 
                                levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
for(i in c(6:8,10:11,14)){
    for(j in c(2:3,16:26)){
        for(k in c(9,27:33)){
            mydata <- data.table(SeqFeature = node.feature.out[,.SD,.SDcols=i], 
                                 GraphFeature = node.feature.out[,.SD,.SDcols=j],
                                 Group = node.feature.out[,.SD,.SDcols=k])
            myscatterplot(datay = mydata, xlabs = colnames(node.feature.out)[i], 
                          ylabs = colnames(node.feature.out)[j], pre = paste(args[1], "node", 
                                                  colnames(node.feature.out)[i], 
                                                  colnames(node.feature.out)[j], "by", 
                                                  colnames(node.feature.out)[k], sep = "."))
        }
    }
}
for(i in c(6:8,10:11,14)){
    for(j in c(4)){
        for(k in c(9,27:33)){
            mydata <- data.table(SeqFeature = node.feature.out[,.SD,.SDcols=i], 
                                 GraphFeature = log2(node.feature.out[,.SD,.SDcols=j] + 1),
                                 Group = node.feature.out[,.SD,.SDcols=k])
            myscatterplot(datay = mydata, xlabs = colnames(node.feature.out)[i], 
                          ylabs = colnames(node.feature.out)[j], pre = paste(args[1], "node", 
                                                                             colnames(node.feature.out)[i], 
                                                                             colnames(node.feature.out)[j], "by", 
                                                                             colnames(node.feature.out)[k], sep = "."))
        }
    }
}

gene.feature.out$Type <- factor(gene.feature.out$Type, 
                                levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
for(i in c(2:4,6:7)){
    for(j in c(9:10,14,16:26)){
        for(k in c(5,27:33)){
            mydata <- data.table(SeqFeature = gene.feature.out[,.SD,.SDcols=i], 
                                 GraphFeature = gene.feature.out[,.SD,.SDcols=j],
                                 Group = gene.feature.out[,.SD,.SDcols=k])
            myscatterplot(datay = mydata, xlabs = colnames(gene.feature.out)[i], 
                          ylabs = colnames(gene.feature.out)[j], pre = paste(args[1], "node", 
                                                                             colnames(gene.feature.out)[i], 
                                                                             colnames(gene.feature.out)[j], "by", 
                                                                             colnames(gene.feature.out)[k], sep = "."))
        }
    }
}
for(i in c(2:4,6:7)){
    for(j in c(11)){
        for(k in c(5,27:33)){
            mydata <- data.table(SeqFeature = gene.feature.out[,.SD,.SDcols=i], 
                                 GraphFeature = log2(gene.feature.out[,.SD,.SDcols=j] + 1),
                                 Group = gene.feature.out[,.SD,.SDcols=k])
            myscatterplot(datay = mydata, xlabs = colnames(gene.feature.out)[i], 
                          ylabs = colnames(gene.feature.out)[j], pre = paste(args[1], "node", 
                                                                             colnames(gene.feature.out)[i], 
                                                                             colnames(gene.feature.out)[j], "by", 
                                                                             colnames(gene.feature.out)[k], sep = "."))
        }
    }
}

# correlation test ------------------------------------------------------------
dataM <- data.matrix(node.feature.out[Type!="Enh",c(2:4,6:8,10:11,14,16:26)])
dataM[is.na(dataM)] <- 0
node.cor <- cor(dataM, method = "spearman")
# corrplot(node.cor, method = "circle", type = "upper", 
#          col = colorRampPalette(c("blue", "white", "red"))(100))
pdf(file = paste("corrplot",args[1], "node.NAas0.pdf", sep = "."), width = 8, height = 8)
corrplot(node.cor, method = "circle", type = "upper", order = "hclust")
dev.off()
# corrplot(node.cor, method = "circle", type = "upper")
pdf(file = paste("corrplot",args[1], "node.omitNA.pdf", sep = "."), width = 8, height = 8)
corrplot(cor(na.omit(data.matrix(node.feature.out[Type!="Enh",c(2:4,6:8,10:11,14,16:26)])), 
             method = "spearman"), method = "circle", type = "upper", order = "hclust")
dev.off()
# corrplot(cor(na.omit(data.matrix(node.feature.out[Type!="Enh",c(9,27:33)])), 
             # method = "kendall"), method = "circle", type = "upper")


dataM2 <- data.matrix(gene.feature.out[,c(2:4,6:7,9:11,14,16:26)])
dataM2[is.na(dataM2)] <- 0
gene.cor <- cor(dataM2, method = "spearman")
pdf(file = paste("corrplot",args[1], "gene.NAas0.pdf", sep = "."), width = 8, height = 8)
corrplot(gene.cor, method = "circle", type = "upper", order = "hclust")
dev.off()
pdf(file = paste("corrplot",args[1], "gene.omitNA.pdf", sep = "."), width = 8, height = 8)
corrplot(cor(na.omit(data.matrix(gene.feature.out[,c(2:4,6:7,9:11,14,16:26)])), 
             method = "spearman"), method = "circle", type = "upper", order = "hclust")
dev.off()

# dataM3 <- gene.feature.out[,c(5,2:4,6:7,9:11,14,16:26)]
# dataM3[Type=="HKG.1", TypeRank:=1]
# dataM3[Type=="HKG.2", TypeRank:=2]
# dataM3[Type=="Other", TypeRank:=3]
# dataM3[Type=="EBV.TSG", TypeRank:=4]
# dataM3[Type=="other.TSG", TypeRank:=4]
# corrplot(cor(na.omit(data.matrix(dataM3[,-1])), 
#              method = "spearman"), method = "circle", type = "upper", order = "hclust")

# datax <- na.omit(data.matrix(dataM3[,-1]))
# v <- cor(x = datax[,-ncol(datax)], y = datax[,ncol(datax)], method = "spearman")
