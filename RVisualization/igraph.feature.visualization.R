library(ggplot2)
library(reshape2)
library(data.table)
library(pheatmap)
library(corrplot)
# args <- commandArgs(TRUE)
# setwd("~/Google Drive File Stream/My Drive/housekeeping_genes/interaction.community/igraph")
args <- c("interaction.SRR5831489", "interaction.SRR5831490", 
          "interaction.SRR5831511", "interaction.SRR5831512")

# boxplot ------------------------------------------------------------
# gene.filename <- paste(args[i], "igraph.gene.feature.out.tsv", sep = ".")
# node.features <- fread(gene.filename, header = T, sep = "\t")
myboxplot <- function(datay = datay, ylabs = ylab, pre = args[1]){
    require(ggplot2)
    colnames(datay) <- c("Class", "Quant", "value")
    ymin <- quantile(as.numeric(datay$value), probs = 0.01, na.rm = T)
    ymax <- quantile(as.numeric(datay$value), probs = 0.99, na.rm = T)
    myplot <- ggplot(datay, aes(x = Class, y = as.numeric(value), fill = Quant)) + 
        geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
        labs(y = ylabs) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    pdf(file = paste("sidebyside.boxplot", pre, "pdf", sep = "."), width = 10, height = 8)
    print(myplot)
    dev.off()
}
for(i in 1:length(args)){
    node.filename <- paste(args[i], "igraph.node.feature.out.tsv", sep = ".")
    node.features <- fread(node.filename, header = T, sep = "\t")
    node.features[,Flag:=paste(Type, CpGquantile, sep = " & ")]
    node.features$Type <- factor(node.features$Type, 
                                 levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
    for(j in c(2,10,14,21,26)){
        mydata <- data.table(Class = node.features$Type, 
                             Quant = node.features$CpGquantile, 
                             value = node.features[,.SD,.SDcols=j])
        myboxplot(datay = mydata, ylabs = colnames(node.features)[j],
                  pre = paste(args[i], "node", 
                              colnames(node.features)[j], "by",
                              "Type_CpGCount", sep = "."))
    }
    for(j in c(2,10,14,21,26)){
        mydata <- data.table(Class = node.features$Exprquantile, 
                             Quant = node.features$CpGquantile, 
                             value = node.features[,.SD,.SDcols=j])
        myboxplot(datay = mydata, ylabs = colnames(node.features)[j],
                  pre = paste(args[i], "node", 
                              colnames(node.features)[j], "by",
                              "Expr_CpGCount", sep = "."))
    }
}

# node.features$Flag <- factor(node.features$Flag, levels = sort(unique(node.features$Flag)))

for(i in 1:length(args)){
    node.filename <- paste(args[i], "igraph.gene.feature.out.tsv", sep = ".")
    node.features <- fread(node.filename, header = T, sep = "\t")
    node.features <- node.features[grep("\\|MT-", Gene, invert = T)]
    node.features$Type <- factor(node.features$Type, 
                                 levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
    for(j in c(6,9,14,21,26)){
        mydata <- data.table(Class = node.features$Type, 
                             Quant = node.features$CpGquantile, 
                             value = node.features[,.SD,.SDcols=j])
        colnames(mydata) <- c("Class", "Quant", "value")
        mydata[is.na(value), value:=0]
        myboxplot(datay = mydata, ylabs = colnames(node.features)[j],
                  pre = paste(args[i], "gene", 
                              colnames(node.features)[j], "by",
                              "Type_CpGCount", sep = "."))
    }
    for(j in c(6,9,14,21,26)){
        mydata <- data.table(Class = node.features$Exprquantile, 
                             Quant = node.features$CpGquantile, 
                             value = node.features[,.SD,.SDcols=j])
        colnames(mydata) <- c("Class", "Quant", "value")
        mydata[is.na(value), value:=0]
        myboxplot(datay = mydata, ylabs = colnames(node.features)[j],
                  pre = paste(args[i], "gene", 
                              colnames(node.features)[j], "by",
                              "Expr_CpGCount", sep = "."))
    }
}

# pheatmap ------------------------------------------------------------
for(i in 1:length(args)){
    node.filename <- paste(args[i], "igraph.node.feature.out.tsv", sep = ".")
    node.features <- fread(node.filename, header = T, sep = "\t")
    for(j in c(9,27:28,30)){
        countdata <- table(node.features[CommunitySize > 5,.SD,.SDcols=c(13,j)])
        pdf(file = paste("community.component", args[i], colnames(node.features)[j], "pdf", sep = "."),
            width = 6, height = 10)
        pheatmap(countdata, show_rownames = F, scale = "row")
        dev.off()
    }
}
for(i in 1:length(args)){
    node.filename <- paste(args[i], "igraph.node.feature.out.tsv", sep = ".")
    node.features <- fread(node.filename, header = T, sep = "\t")
    for(j in c(9)){
        countdata <- table(node.features[CommunitySize > 5,.SD,.SDcols=c(13,j)])
        rmax <- apply(countdata[,-c(2,5)], 1, max)
        heatdata <- data.table(HKG.1 = countdata[rmax > 1,3]/typeCount[2],
                               HKG.2 = countdata[rmax > 1,4]/typeCount[3],
                               Other = countdata[rmax > 1,5]/typeCount[4],
                               EBV.TSG = countdata[rmax > 1,1]/typeCount[1],
                               other.TSG = countdata[rmax > 1,6]/typeCount[5])
        pdf(file = paste("community.component.percentage", args[i], colnames(node.features)[j], "pdf", sep = "."),
            width = 6, height = 10)
        pheatmap(heatdata, show_rownames = F, scale = "none")
        dev.off()
    }
}


# combine ------------------------------------------------------------
community.all <- data.frame()
for(i in 1:length(args)){
    node.filename <- paste(args[i], "igraph.gene.feature.out.tsv", sep = ".")
    node.features <- fread(node.filename, header = T, sep = "\t")
    node.features[, expriment := args[i]]
    community.all <- rbind(community.all, node.features)
}
typeCount <- table(node.features$Type)
community.all[,comFlag:=ifelse(is.na(Community), "Lonely", "InCommunity")]
datay <- community.all[,c(5,35,34)]
dataz <- melt(table(datay)/rowSums(table(datay)/4))
dataz$Type <- factor(dataz$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
dataz$comFlag <- factor(dataz$comFlag, levels = c("Lonely", "InCommunity"))
myplot <- ggplot(dataz, aes(x = Type, y = value, fill = comFlag)) + 
    geom_bar(stat = "identity") + facet_grid(.~expriment)
    # theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste("stack.barplot", "InCommunity", "pdf", sep = "."), width = 16, height = 8)
print(myplot)
dev.off()

datay <- community.all[CommunitySize > 5,c(5,35,34)]
dataz <- as.data.table(melt(table(datay)))
for(i in 1:length(typeCount)){
    dataz[Type==names(typeCount)[i], ratio:=value/typeCount[i]]
}
dataz$Type <- factor(dataz$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
myplot <- ggplot(dataz, aes(x = Type, y = ratio, fill = Type)) + 
    geom_bar(stat = "identity") + facet_grid(.~expriment)
pdf(file = paste("stack.barplot", "InCommunity", "comsizegt5.pdf", sep = "."), width = 16, height = 8)
print(myplot)
dev.off()

datay <- community.all[CommunitySize > 10,c(5,35,34)]
dataz <- as.data.table(melt(table(datay)))
for(i in 1:length(typeCount)){
    dataz[Type==names(typeCount)[i], ratio:=value/typeCount[i]]
}
dataz$Type <- factor(dataz$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
myplot <- ggplot(dataz, aes(x = Type, y = ratio, fill = Type)) + 
    geom_bar(stat = "identity") + facet_grid(.~expriment)
pdf(file = paste("stack.barplot", "InCommunity", "comsizegt10.pdf", sep = "."), width = 16, height = 8)
print(myplot)
dev.off()

datay <- community.all[,c(5,28,35,34)]
dataz <- as.data.table(melt(table(datay)))
for(i in 1:length(typeCount)){
    dataz[Type==names(typeCount)[i], ratio:=value/typeCount[i]]
}
dataz[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
dataz$TypeCpG <- factor(dataz$TypeCpG, levels = sort(unique(dataz$TypeCpG)))
dataz$comFlag <- factor(dataz$comFlag, levels = c("Lonely", "InCommunity"))
myplot <- ggplot(dataz, aes(x = TypeCpG, y = ratio, fill = comFlag)) + 
    geom_bar(stat = "identity") + facet_wrap(expriment ~ ., nrow = 2) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste("stack.barplot", "InCommunity", "byTypeCpG.pdf", sep = "."), width = 12, height = 10)
print(myplot)
dev.off()

datay <- community.all[,c(5,14,34)]
datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
ymax <- quantile(datay$CommunitySize, probs = 0.95, na.rm = T)
myplot <- ggplot(datay, aes(x = expriment, y = CommunitySize, fill = Type)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste("boxplot", "CommunitySize", "byType.pdf", sep = "."), width = 10, height = 8)
print(myplot)
dev.off()
datay <- community.all[,c(28,14,34)]
ymax <- quantile(datay$CommunitySize, probs = 0.95, na.rm = T)
myplot <- ggplot(datay, aes(x = expriment, y = CommunitySize, fill = CpGquantile)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
pdf(file = paste("boxplot", "CommunitySize", "byCpG.pdf", sep = "."), width = 10, height = 8)
print(myplot)
dev.off()
datay <- community.all[,c(5,28,14,34)]
datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
ymax <- quantile(datay$CommunitySize, probs = 0.95, na.rm = T)
myplot <- ggplot(datay, aes(x = expriment, y = CommunitySize, fill = TypeCpG)) + 
    geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
pdf(file = paste("boxplot", "CommunitySize", "byTypeCpG.pdf", sep = "."), width = 15, height = 10)
print(myplot)
dev.off()

# combine plot ------------------------------------------------------------
community.all <- data.frame()
for(i in 1:length(args)){
    node.filename <- paste(args[i], "igraph.gene.feature.out.tsv", sep = ".")
    node.features <- fread(node.filename, header = T, sep = "\t")
    node.features[, expriment := args[i]]
    community.all <- rbind(community.all, node.features)
}
typeCount <- table(node.features$Type)
community.all[,comFlag:=ifelse(is.na(Community), "Lonely", "InCommunity")]
community.all[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
community.all$Type <- factor(community.all$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
community.all$TypeCpG <- factor(community.all$TypeCpG, levels = sort(unique(community.all$TypeCpG)))
#### boxplot function
combinboxplot <- function(datap = datay, ylabs = ylab, pre = args[i]){
    require(ggplot2)
    colnames(datap) <- c("Class", "Group", "value")
    ymin <- quantile(as.numeric(datap$value), probs = 0.01, na.rm = T)
    ymax <- quantile(as.numeric(datap$value), probs = 0.99, na.rm = T)
    myplot <- ggplot(datap, aes(x = Group, y = as.numeric(value), fill = Class)) + 
        geom_boxplot(width = 0.8, outlier.alpha = 0.1, outlier.size = 0.1) + 
        scale_y_continuous(limits = c(ymin, ymax)) +
        labs(y = ylabs) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    pdf(file = paste("boxplot", pre, "pdf", sep = "."), width = 15, height = 10)
    print(myplot)
    dev.off()
}
### node degree 
datay <- community.all[,c(5,9,28,34,36)]
combinboxplot(datap = datay[,c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], sep = "."))
combinboxplot(datap = datay[,c(3,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[3], sep = "."))
combinboxplot(datap = datay[,c(5,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[5], sep = "."))
### node closeness
datay <- community.all[,c(5,10,28,34,36)]
combinboxplot(datap = datay[,c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], sep = "."))
combinboxplot(datap = datay[,c(3,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[3], sep = "."))
combinboxplot(datap = datay[,c(5,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[5], sep = "."))
combinboxplot(datap = datay[expriment=="interaction.SRR5831489",c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], "interaction.SRR5831489", sep = "."))
### node betweenness
datay <- community.all[,c(5,11,28,34,36)]
datay[,node.betweenness:=log2(node.betweenness + 1)]
combinboxplot(datap = datay[,c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], sep = "."))
combinboxplot(datap = datay[,c(3,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[3], sep = "."))
combinboxplot(datap = datay[,c(5,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[5], sep = "."))
### PPcount
datay <- community.all[,c(5,22,28,34,36)]
combinboxplot(datap = datay[,c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], sep = "."))
combinboxplot(datap = datay[,c(3,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[3], sep = "."))
combinboxplot(datap = datay[,c(5,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[5], sep = "."))
### PP.percetage
datay <- community.all[,c(5,26,28,34,36)]
combinboxplot(datap = datay[,c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], sep = "."))
combinboxplot(datap = datay[,c(3,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[3], sep = "."))
combinboxplot(datap = datay[,c(5,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[5], sep = "."))
### medPETs
datay <- community.all[,c(5,20,28,34,36)]
combinboxplot(datap = datay[,c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], sep = "."))
combinboxplot(datap = datay[,c(3,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[3], sep = "."))
combinboxplot(datap = datay[,c(5,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[5], sep = "."))
### medDisBin
datay <- community.all[,c(5,21,28,34,36)]
combinboxplot(datap = datay[,c(1,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[1], sep = "."))
combinboxplot(datap = datay[,c(3,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[3], sep = "."))
combinboxplot(datap = datay[,c(5,4,2)], ylabs = colnames(datay)[2], 
              pre = paste(colnames(datay)[2], "by", colnames(datay)[5], sep = "."))

# community component ------------------------------------------------------------
# with enhancer
network.all <- data.frame()
for(i in 1:length(args)){
    node.filename <- paste(args[i], "igraph.node.feature.out.tsv", sep = ".")
    node.features <- fread(node.filename, header = T, sep = "\t")
    node.features[, expriment := args[i]]
    network.all <- rbind(network.all, node.features)
}
network.all$Type <- factor(network.all$Type, 
                           levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
network.all$CpGquantile <- factor(network.all$CpGquantile, 
                           levels = c("Enh", "CpGQ1", "CpGQ2", "CpGQ3", "CpGQ4"))
network.all$Exprquantile <- factor(network.all$Exprquantile, 
                                   levels = c("Enh", "NULL", "ExprQ1", "ExprQ2", "ExprQ3", "ExprQ4"))
#### boxplot function
combinbarplot <- function(datap = datay, ylabs = ylab, pre = args[i]){
    require(ggplot2)
    colnames(datap) <- c("Community", "Type", "value", "experiment")
    myplot <- ggplot(datap, aes(x = Community, y = as.numeric(value), fill = Type)) + 
        geom_bar(stat = "identity") + facet_grid( . ~ experiment) +
        labs(y = ylabs) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    pdf(file = paste("stack.barplot", pre, "pdf", sep = "."), width = 20, height = 6)
    print(myplot)
    dev.off()
}

sidebysidebarplot <- function(datap = datay, ylabs = ylab, pre = args[i]){
    require(ggplot2)
    colnames(datap) <- c("Community", "Type", "value", "experiment")
    myplot <- ggplot(datap, aes(x = Community, y = as.numeric(value), fill = Type)) + 
        geom_bar(stat = "identity", position = "dodge") + facet_grid(experiment ~ .) +
        labs(y = ylabs) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    pdf(file = paste("sidebyside.barplot", pre, "pdf", sep = "."), width = 15, height = 10)
    print(myplot)
    dev.off()
}

for(j in c(9,28,30)){
    topcoms <- data.table()
    # botcoms <- data.table()
    randcoms <- data.table()
    for(i in 1:length(args)){
        tempdata <- network.all[expriment==args[i]]
        countdata <- table(tempdata[,.SD,.SDcols=c(13,j)])
        comsize <- rowSums(countdata)
        countdata <- countdata[,grep("Enh", colnames(countdata), invert = T)]
        topcom <- as.data.table(melt(countdata[order(comsize, decreasing = T)[1:30],]
                                     /comsize[order(comsize, decreasing = T)[1:30]]))
        topcom[,Community:=rep(paste("top",1:30), ncol(countdata))]
        topcoms <- rbind(topcoms, data.table(topcom, experiment = args[i]))
        
        # botcom <- as.data.table(melt(countdata[order(comsize, decreasing = F)[1:30],]
        #                              /comsize[order(comsize, decreasing = F)[1:30]]))
        # botcom[,Community:=rep(paste("bottom",1:30), ncol(countdata))]
        # botcoms <- rbind(botcoms, data.table(botcom, experiment = args[i]))

        # rands <- sample(1:nrow(countdata), 30)
        # randcom <- as.data.table(melt(countdata[rands,]/comsize[rands]))
        # randcom[,Community:=rep(paste("rand",1:30), ncol(countdata))]
        # randcoms <- rbind(randcoms, data.table(randcom, experiment = args[i]))
        
        total <- colSums(countdata)
        randcom <- as.data.table(data.frame(Community = "total", 
                                            Type = names(total), value = total/nrow(tempdata)))
        randcoms <- rbind(randcoms, data.table(randcom, experiment = args[i]))
    }
    topcoms$Community <- factor(topcoms$Community, levels = unique(topcoms$Community))
    combinbarplot(datap = topcoms, ylabs = colnames(network.all)[j], 
                  pre = paste(colnames(network.all)[j], "by", "community", sep = "."))
    tmp <- topcoms[grep("Other", as.matrix(topcoms[,2]), invert = T)]
    tmp <- tmp[grep("NULL", as.matrix(tmp[,2]), invert = T)]
    sidebysidebarplot(datap = tmp, 
                      ylabs = colnames(network.all)[j], 
                  pre = paste(colnames(network.all)[j], "by", "community", sep = "."))
    colnames(randcoms) <- colnames(topcoms)
    mixcoms <- rbind(topcoms, randcoms)
    # mixcoms <- rbind(topcoms, botcoms)
    mixcoms$Community <- factor(mixcoms$Community, levels = unique(mixcoms$Community))
    tmp <- mixcoms[grep("Other", as.matrix(mixcoms[,2]), invert = T)]
    tmp <- tmp[grep("NULL", as.matrix(tmp[,2]), invert = T)]
    combinbarplot(datap = mixcoms, ylabs = colnames(network.all)[j], 
                  pre = paste(colnames(network.all)[j], "by", "community.mix", sep = "."))
    sidebysidebarplot(datap = tmp, ylabs = colnames(network.all)[j], 
                  pre = paste(colnames(network.all)[j], "by", "community.mix", sep = "."))
}


######## no meaning
for(j in c(9)){
    topcoms <- data.table()
    for(i in 1:length(args)){
        tempdata <- network.all[expriment==args[i]]
        countdata <- table(tempdata[,.SD,.SDcols=c(13,j)])
        comsize <- rowSums(countdata)
        countdata <- countdata[,grep("Enh", colnames(countdata), invert = T)]
        topcom <- as.data.table(melt(countdata[order(comsize, decreasing = T)[1:30],]))
        topcom[,value:=value/typeCount[match(Type,names(typeCount))]]
        topcom[,Community:=rep(paste("top",1:30), ncol(countdata))]
        topcoms <- rbind(topcoms, data.table(topcom, experiment = args[i]))
    }
    topcoms$Community <- factor(topcoms$Community, levels = unique(topcoms$Community))
    combinbarplot(datap = topcoms, ylabs = colnames(network.all)[j], 
                  pre = paste(colnames(network.all)[j], "percentage", "by", "community", sep = "."))
    tmp <- topcoms[grep("Other", as.matrix(topcoms[,2]), invert = T)]
    sidebysidebarplot(datap = tmp, 
                      ylabs = colnames(network.all)[j], 
                      pre = paste(colnames(network.all)[j], "percentage", "by", "community", sep = "."))
}


# datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
# ymax <- quantile(datay$node.degree, probs = 0.99, na.rm = T)
# myplot <- ggplot(datay, aes(x = expriment, y = node.degree, fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.degree", "byType.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = node.degree, fill = CpGquantile)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.degree", "byCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = node.degree, fill = TypeCpG)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.degree", "byTypeCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# ### node closeness 
# datay <- community.all[,c(5,10,28,34)]
# datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
# ymin <- quantile(as.numeric(datay[expriment=="interaction.SRR5831489"]$node.closeness), probs = 0.01, na.rm = T)
# ymax <- quantile(as.numeric(datay[expriment=="interaction.SRR5831489"]$node.closeness), probs = 0.99, na.rm = T)
# myplot <- ggplot(datay[expriment=="interaction.SRR5831489"], aes(x = expriment, y = as.numeric(node.closeness), fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.closeness", "byType.pdf", sep = "."), width = 6, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay[expriment=="interaction.SRR5831489"], aes(x = expriment, y = as.numeric(node.closeness), fill = CpGquantile)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.closeness", "byCpG.pdf", sep = "."), width = 6, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay[expriment=="interaction.SRR5831489"], aes(x = expriment, y = as.numeric(node.closeness), fill = TypeCpG)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.closeness", "byTypeCpG.pdf", sep = "."), width = 6, height = 10)
# print(myplot)
# dev.off()
# ### node betweenness 
# datay <- community.all[,c(5,11,28,34)]
# datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
# ymin <- quantile(log2(as.numeric(datay[expriment=="interaction.SRR5831489"]$node.betweenness)+1), probs = 0.01, na.rm = T)
# ymax <- quantile(log2(as.numeric(datay[expriment=="interaction.SRR5831489"]$node.betweenness)+1), probs = 0.99, na.rm = T)
# myplot <- ggplot(datay[expriment=="interaction.SRR5831489"], aes(x = expriment, y = log2(as.numeric(node.betweenness)), fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.betweenness", "byType.pdf", sep = "."), width = 6, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay[expriment=="interaction.SRR5831489"], aes(x = expriment, y = log2(as.numeric(node.betweenness)), fill = CpGquantile)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.betweenness", "byCpG.pdf", sep = "."), width = 6, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay[expriment=="interaction.SRR5831489"], aes(x = expriment, y = log2(as.numeric(node.betweenness)), fill = TypeCpG)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(ymin, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "node.betweenness", "byTypeCpG.pdf", sep = "."), width = 6, height = 10)
# print(myplot)
# dev.off()
# ### PPcount
# datay <- community.all[,c(5,22,28,34)]
# datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
# ymax <- quantile(datay$PPCount, probs = 0.99, na.rm = T)
# myplot <- ggplot(datay, aes(x = expriment, y = PPCount, fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "PPCount", "byType.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = PPCount, fill = CpGquantile)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "PPCount", "byCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = PPCount, fill = TypeCpG)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "PPCount", "byTypeCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# ### PP.percetage
# datay <- community.all[,c(5,26,28,34)]
# datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
# ymax <- quantile(datay$PP.percetage, probs = 0.99, na.rm = T)
# myplot <- ggplot(datay, aes(x = expriment, y = PP.percetage, fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "PP.percetage", "byType.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = PP.percetage, fill = CpGquantile)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "PP.percetage", "byCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = PP.percetage, fill = TypeCpG)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "PP.percetage", "byTypeCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# 
# ### medPETs
# datay <- community.all[,c(5,20,28,34)]
# datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
# ymax <- quantile(datay$medPETs, probs = 0.95, na.rm = T)
# myplot <- ggplot(datay, aes(x = expriment, y = medPETs, fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "medPETs", "byType.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = medPETs, fill = CpGquantile)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "medPETs", "byCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = medPETs, fill = TypeCpG)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "medPETs", "byTypeCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# 
# ### medDisBin
# datay <- community.all[,c(5,21,28,34)]
# datay[,TypeCpG:=paste(Type, CpGquantile, sep = "&")]
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# datay$TypeCpG <- factor(datay$TypeCpG, levels = sort(unique(datay$TypeCpG)))
# ymax <- quantile(datay$medDisBin, probs = 0.95, na.rm = T)
# myplot <- ggplot(datay, aes(x = expriment, y = medDisBin, fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "medDisBin", "byType.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = medDisBin, fill = CpGquantile)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "medDisBin", "byCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = medDisBin, fill = TypeCpG)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax)) +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
# pdf(file = paste("boxplot", "medDisBin", "byTypeCpG.pdf", sep = "."), width = 15, height = 10)
# print(myplot)
# dev.off()
# # clique combine ------------------------------------------------------------
# clique.all <- data.frame()
# for(i in 1:length(args)){
#     node.filename <- paste(args[i], "igraph.node2clique.out.tsv", sep = ".")
#     node.features <- fread(node.filename, header = T, sep = "\t")
#     node.features[, expriment := args[i]]
#     clique.all <- rbind(clique.all, node.features)
# }
# ### nothing interesting
# datay <- unique(clique.all[,c(2,4,12,37)])
# datay$Type <- factor(datay$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG", "Enh"))
# ymax <- quantile(datay$CliqueSize, probs = 0.95, na.rm = T)
# myplot <- ggplot(datay, aes(x = expriment, y = CliqueSize, fill = Type)) + 
#     geom_boxplot(width = 0.8) + scale_y_continuous(limits = c(0, ymax))
# pdf(file = paste("boxplot", "CliqueSize", "byType.pdf", sep = "."), width = 10, height = 8)
# print(myplot)
# dev.off()
# 
# # if(length(vx[grep("ENSG", vx)]) >= 3)
# # so the clieuq is alredy PP-clique
# datay <- unique(clique.all[,c(4,12,37)])
# dataz <- as.data.table(melt(table(datay[Type!="Enh",2:3])))
# for(i in 1:length(typeCount)){
#     dataz[Type==names(typeCount)[i], ratio:=value/typeCount[i]]
# }
# dataz$Type <- factor(dataz$Type, levels = c("HKG.1", "HKG.2", "Other", "EBV.TSG", "other.TSG"))
# myplot <- ggplot(dataz, aes(x = Type, y = ratio, fill = expriment)) + 
#     geom_bar(stat = "identity", position = "dodge") 
# pdf(file = paste("stack.barplot", "InClique", "pdf", sep = "."), width = 10, height = 8)
# print(myplot)
# dev.off()
# 
# datay <- data.table(Clique = clique.all$Clique, 
#                     CliqueSize = as.numeric(clique.all$CliqueSize), 
#                     Gene = clique.all$node.name, 
#                     CpG.Count = as.numeric(clique.all$CpG.Count),
#                     Expression = as.numeric(clique.all$EBV.medianexpression),
#                     expriment = clique.all$expriment)
# datay[, MedianCpG:=median(CpG.Count, na.rm = T), by = list(Clique, expriment)]
# datay[, MedianExpr:=median(Expression, na.rm = T), by = list(Clique, expriment)]
# datay[, sdCpG:=sd(CpG.Count, na.rm = T), by = list(Clique, expriment)]
# datay[, sdExpr:=sd(Expression, na.rm = T), by = list(Clique, expriment)]
# dataz <- unique(datay[,c(1:2,6:10)])
# myplot <- ggplot(datay, aes(x = expriment, y = MedianCpG, fill = factor(CliqueSize))) + 
#     geom_boxplot(width = 0.8)
# pdf(file = paste("boxplot", "MedianCpG", "byCliuqeSize.pdf", sep = "."), width = 10, height = 8)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = MedianExpr, fill = factor(CliqueSize))) + 
#     geom_boxplot(width = 0.8)
# pdf(file = paste("boxplot", "MedianExpr", "byCliuqeSize.pdf", sep = "."), width = 10, height = 8)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = sdCpG, fill = factor(CliqueSize))) + 
#     geom_boxplot(width = 0.8)
# pdf(file = paste("boxplot", "sdCpG", "byCliuqeSize.pdf", sep = "."), width = 10, height = 8)
# print(myplot)
# dev.off()
# myplot <- ggplot(datay, aes(x = expriment, y = sdExpr, fill = factor(CliqueSize))) + 
#     geom_boxplot(width = 0.8)
# pdf(file = paste("boxplot", "sdExpr", "byCliuqeSize.pdf", sep = "."), width = 10, height = 8)
# print(myplot)
# dev.off()

