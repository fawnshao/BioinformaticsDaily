# source("https://bioconductor.org/biocLite.R")
# biocLite("RUVSeq")
library(RUVSeq)
library(RColorBrewer)
library(tximport)
library(readr)
tx2gene <- read_tsv("/home1/04935/shaojf/scratch/Salmon.index/tx2gene.gencode.v19.tsv")
colors <- brewer.pal(3, "Set2")
args <- commandArgs(TRUE)
# args <- c("tximport.samples.txt", "YY1del")
input <- read.table(args[1], header = T, sep = " ", row.names = 1)
myfiles <- as.vector(input$files)
names(myfiles) <- rownames(input)
txi <- tximport(myfiles, type = "salmon", tx2gene = tx2gene)
#### Salmon TPM
write.table(x = data.frame(rownames(txi$abundance), txi$abundance), 
	file = paste(args[2],"salmon.tximport.abundance.tsv", sep = "."), 
	sep = "\t", row.names = F, quote = F)
estimatecount <- txi$counts
# filter <- apply(input, 1, function(x) length(x[x>5])>=2)
# filtered <- zfGenes[filter,]
# genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
# spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
genes <- rownames(estimatecount)[grep("^ERCC-", rownames(estimatecount), invert = TRUE)]
spikes <- rownames(estimatecount)[grep("^ERCC-", rownames(estimatecount))]
x <- as.factor(c(rep("ctrl", 2), rep("YY1del", 2)))
set <- newSeqExpressionSet(as.matrix(round(estimatecount,0)), 
	phenoData = data.frame(x, row.names = colnames(estimatecount)))
# set

set <- betweenLaneNormalization(set, which = "upper")
set1 <- RUVg(set, spikes, k = 1)
pData(set1)
pdf(file = paste(args[2],"RUV.plotRLE.pdf", sep = "."))
plotRLE(set1, outline = FALSE, ylim = c(-4, 4), col = colors[x])
dev.off()
pdf(file = paste(args[2],"RUV.plotPCA.pdf", sep = "."))
plotPCA(set1, col = colors[x], cex = 1.2)
dev.off()
nc <- normCounts(set1)
write.table(x = data.frame(rownames(nc), nc), 
	file = paste(args[2],"RUV.normCounts.tsv", sep = "."), 
	sep = "\t", row.names = F, quote = F)


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                               colData = pData(set1),
                               design = ~ W_1 + x)
dds <- DESeq(dds)
res <- results(dds)
# res

out <- data.frame(counts(dds, normalized = F), counts(dds, normalized = T), res)
write.table(x = data.frame(rownames(out), out), 
	file = paste(args[2],"RUV.tsv", sep = "."), 
	sep = "\t", row.names = F, quote = F)
