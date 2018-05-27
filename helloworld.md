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

rquantile=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/data.quantile.R
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

for f in exp_*.tsv.gz
do
    gunzip -c $f | cut -f 2-3,8-9 | grep -wf NET1.ids > NET1.$f.txt
done

for f in `ll NET1.exp_*.tsv.gz.txt | awk '$5 > 0 {print $9}'`
do
    perl $myperl specimen.info $f 0 1 | cut -f 1-4,7 > $f.specimen
done

for f in NET1.exp_*.tsv.gz.txt.specimen
do
    Rscript $rquantile $f
done

library(data.table)
library(ggplot2)
args <- c("NET1.exp_array.BRCA-US.tsv.gz.txt.specimen", "NET1.exp_seq.BRCA-US.tsv.gz.txt.specimen")
for(i in 1:2){
    datax <- fread(args[i], sep = "\t", header = F)
    colnames(datax) <- c("cancers", "samples", "genes", "NormalizedValues", "types")

    ymax <- quantile(datax$NormalizedValues, probs = 0.95)
    datax$types <- factor(datax$types, levels = c("Normal - tissue adjacent to primary", "Primary tumour - solid tissue", "Metastatic tumour - metastasis to distant location"))
    myplot <- ggplot(data = datax, aes(x = types, y = NormalizedValues, fill = types)) + 
        geom_boxplot(position = position_dodge(1), outlier.alpha = 0.1, outlier.size = 0.1) + 
        scale_y_continuous(limits = c(0, ymax)) +
        # facet_wrap(. ~ genes) +
        labs(title = args[1], caption = date()) + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
    png(filename = paste(args[i], "boxplot.png", sep = "."), width = 600, height = 800)
    print(myplot)
    dev.off()
}


