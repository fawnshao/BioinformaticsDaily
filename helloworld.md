st-2.compute.amazonaws.com
ssh -i "amazonKey.pem" ubuntu@18.217.196.213
# sample Markdown
library(data.table)
library(ggplot2)
args <- list.files(pattern = ".tsv.gz.txt.specimen$")
for(z in 1:length(args)){
    datax <- fread(args[z], sep = "\t", header = F)
    colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues", "types")

    tumors <- datax[!grep("Normal", types),]
    if(length(unique(datax$samples)) > 100){
        groups <- unique(datax$genes)
        for(i in 1:length(groups)){
            a <- tumors[ genes == groups[i] ]
            tumors.upper <- a[ NormalizedValues > quantile(NormalizedValues, probs = 0.8),]
            tumors.lower <- a[ NormalizedValues < quantile(NormalizedValues, probs = 0.2),]

            write.table(tumors.upper, row.names = F, file = paste(args[z], groups[i], "tumors.upper", "tsv", sep = "."), sep = "\t")
            write.table(tumors.lower, row.names = F, file = paste(args[z], groups[i], "tumors.lower", "tsv", sep = "."), sep = "\t")
        }
    }
}
