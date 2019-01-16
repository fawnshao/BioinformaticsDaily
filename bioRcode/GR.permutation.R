library(GenomicRanges)
library(data.table)
bedfile <- fread("hg19.refGene.tss.uniq.srt.bed", header = F)
bedfile <- bedfile[grep("MIR",V4, invert = T)]
oncogenes <- fread("hg19.refGene.tss.Cosmic.CancerGeneCensus.oncogene.gene.bed", header = F)
pmuts <- fread("13cancers.pMUT", header = F)
allmuts <- fread("13cancers.allmut.stats")

my.reftss <- GRanges(seqnames = bedfile$V1, ranges = IRanges(start = bedfile$V2, 
                                                             end = bedfile$V3, names = bedfile$V4), 
                     strand = bedfile$V6)
oncogenetss <- GRanges(seqnames = oncogenes$V1, ranges = IRanges(start = oncogenes$V2, 
                                                                 end = oncogenes$V3, names = oncogenes$V4), 
                     strand = oncogenes$V6)

