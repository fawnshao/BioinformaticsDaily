# wilcox.test(x = rnorm(100, 50, 5), y = rnorm(60, 30, 5))$p.value
library(data.table)
args <- commandArgs(TRUE)
# args <- c("oncogene.200kb.nononcogene.neighbors.pMUT.forPermutation.txt",
#           "rand273.200kb.nonpseudoncogene.neighbors.pMUT.forPermutation.txt",
#           "onco.pvals.txt")
datax <- fread(args[1], header = F, na.strings = "/")
datay <- fread(args[2], header = F, na.strings = "/")
datax[,freq:=V3/V4]
datay[,freq:=V3/V4]
types <- unique(datax$V1)
pvals <- c()
for(i in 1:length(types)){
    pvals[i] <- wilcox.test(x = datax[V1==types[i]]$freq, y = datay[V1==types[i]]$freq)$p.value
}
res <- data.frame(Types = types, pvals = pvals)
write.table(res, file = args[3], append = T,  sep = "\t", row.names = F, col.names = F)