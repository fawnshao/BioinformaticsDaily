args <- commandArgs(TRUE)
# args <- c("rep12.Hela.gene.eRNA.counts", "siCTL", "siCTL", "siTIP60", "siTIP60")
# cds = newCountDataSet(CountTable, Design$condition )
# estimateSizeFactors(cds[which(cds$feature=="snoRNA"])
# sizeFactors( cds )
# head(counts( cds, normalized=TRUE))
library(DESeq2)
inputfile <- args[1]
inputdata <- read.table(inputfile, header = T, row.names = 1, sep = "\t")

coldata <- data.frame(condition = args[2:length(args)])
rownames(coldata) <- colnames(inputdata)

dds <- DESeqDataSetFromMatrix(countData = inputdata,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
cmp <- paste("condition", coldata[nrow(coldata),1], "vs", coldata[1,1], sep="_")
res <- results(dds, name = cmp)
# or to shrink log fold changes association with condition:
# res.1 <- lfcShrink(dds, coef = cmp, type = "apeglm")
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
normres <- counts(dds, normalized = TRUE)
output <- data.frame(rownames(inputdata), inputdata, assay(vsd), assay(rld), normres, res)
colnames(output) <- c("ID", paste("raw", colnames(inputdata), sep = "."), 
	paste("vsd", colnames(inputdata), sep = "."), paste("rld", colnames(inputdata), sep = "."),
	paste("normalized", colnames(inputdata), sep = "."), colnames(res))
write.table(output, file = paste(args[1], "DESeq2.out.tsv", sep = "."), 
	quote = F, sep = "\t", row.names = F)