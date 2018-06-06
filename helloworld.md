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

bedtools closest -D a -a <(awk '$2~/hkg/' hkg.tsg.srtbyPCA.transcript.bed | cut -f 3- | awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' | bedtools sort -i -) -b <(bedtools sort -i hg19-tRNAs.bed) > hkg.tRNA.txt
awk '$7!="." && $19 < 1000 && $19 > -1000' hkg.tRNA.txt

bedtools closest -D a -a <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed | bedtools sort -i -) -b <(bedtools sort -i hg19-tRNAs.bed) > gencode.v19.transcript.tRNA.txt
awk '$7!="." && $19 < 1000 && $19 > -1000' gencode.v19.transcript.tRNA.txt | wc -l

bedtools closest -D a -a <(bedtools sort -i hg19-tRNAs.bed) -b <(awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' gencode.v19.transcript.bed | bedtools sort -i -) > tRNA.gencode.v19.transcript.txt
bedtools closest -D a -a <(bedtools sort -i hg19-tRNAs.bed) -b <(cut -f 3- hkg.tsg.srtbyPCA.transcript.bed | awk -vOFS="\t" '{a=$2-1000;b=$2+1000;if($6=="-"){a=$3-1000;b=$3+1000}if(a<0){a=0;}print $1,a,b,$4,$5,$6}' | bedtools sort -i -) > tRNA.hkg.tsg.txt

perl /home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/make_binarymatrix_from_single_file.pl <(awk '{print $4"%"$1":"$2"-"$3"\t"$5}' gencode.v19.enhanceratlas.txt) > gencode.v19.enhanceratlas.ep.mat

head -1 gencode.v19.enhanceratlas.ep.mat | awk '{print "Gene\tType\t"$0}' > hkg.tsg.srtbyPCA.enhanceratlas.ep.mat
perl $myperl <(sed 's/%/\t/' gencode.v19.enhanceratlas.ep.mat | tail -n +2) <(tail -n +2 hkg.tsg.srtbyPCA.class) 0 0 | cut -f 1-2,4- >> hkg.tsg.srtbyPCA.enhanceratlas.ep.mat
awk -F"\t" -vOFS="\t" '{sum=0;for(i=4;i<=NF;i++){sum+=$i}print $1,$2,$3,sum}' hkg.tsg.srtbyPCA.enhanceratlas.ep.mat | sort -k4,4nr > hkg.tsg.srtbyPCA.enhanceratlas.ep.mat.stats
awk -F"\t" '{sum=0;for(i=4;i<=NF;i++){sum+=$i}if(sum > 10 || NR==1 || ($2!/hkg/ && sum >3)){print $0}}' hkg.tsg.srtbyPCA.enhanceratlas.ep.mat > hkg.tsg.srtbyPCA.enhanceratlas.ep.high.mat

# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/enhanceratlas.ep.interaction/ mygoogle:hkg_tsg/both.pc.and.nc.genes/enhanceratlas.ep.interaction/

library(data.table)
library(pheatmap)
args <- c("hkg.tsg.srtbyPCA.enhanceratlas.ep.high.mat")
input <- fread(args[1], sep = "\t", header = T, na.strings = "/")
scores <- data.matrix(input[,-c(1:3)])
rownames(scores) <- paste(as.matrix(input[,1]), as.matrix(input[,3]), sep = "%")
annosR <- input[,2]
rownames(annosR) <- rownames(scores)
colors <- colorRampPalette(c("white", "blue"))(100)
png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
myplot <- pheatmap(scores, scale = "none", annotation_row = annosR,
    show_rownames = F, show_colnames = T, color = colors, 
    cluster_cols = F, cluster_rows = F)
dev.off()

awk '$NF==0' hkg.tsg.srtbyPCA.enhanceratlas.ep.mat.stats | grep -v "MT-" > hkg.tsg.srtbyPCA.enhanceratlas.ep.na.stat
awk '$NF>20' hkg.tsg.srtbyPCA.enhanceratlas.ep.mat.stats | grep "hkg" > hkg.tsg.srtbyPCA.enhanceratlas.ep.hkggt20.stat
awk '$NF<2' hkg.tsg.srtbyPCA.enhanceratlas.ep.mat.stats | grep "hkg" | grep -v "MT-"  | grep -vwf <(cut -f 1 hkg.tsg.srtbyPCA.enhanceratlas.ep.gt20.stat | sort | uniq) > hkg.tsg.srtbyPCA.enhanceratlas.ep.hkglt2.stat

# cut -f 1 hkg.tsg.srtbyPCA.enhanceratlas.ep.hkglt2.stat | sort | uniq | wc -l
perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(cut -f 1 hkg.tsg.srtbyPCA.enhanceratlas.ep.hkggt20.stat | sort | uniq) 0 0 | cut -f 4- > hkg.tsg.srtbyPCA.enhanceratlas.ep.hkggt20.gene.bed

perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(cut -f 1 hkg.tsg.srtbyPCA.enhanceratlas.ep.hkglt2.stat | sort | uniq) 0 0 | cut -f 4- > hkg.tsg.srtbyPCA.enhanceratlas.ep.hkglt2.gene.bed

perl $myperl hkg.tsg.srtbyPCA.transcript.bed <(cut -f 1 hkg.tsg.srtbyPCA.enhanceratlas.ep.na.stat | sort | uniq) 0 0 | cut -f 4- > hkg.tsg.srtbyPCA.enhanceratlas.ep.na.gene.bed

# cut -f 3 hkg.tsg.srtbyPCA.enhanceratlas.ep.hkggt20.stat | sort | uniq -c | awk '$1>1' | more
cut -f 3 hkg.tsg.srtbyPCA.enhanceratlas.ep.hkggt20.stat | sort | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print $0"\t"NR}' > hkg.tsg.srtbyPCA.enhanceratlas.ep.hkggt20.enhancer.bed

cut -f 3 hkg.tsg.srtbyPCA.enhanceratlas.ep.hkglt2.stat | sort | uniq | sed 's/:/\t/;s/-/\t/' | awk '{print $0"\t"NR}' > hkg.tsg.srtbyPCA.enhanceratlas.ep.hkglt2.enhancer.bed

# rclone sync /home1/04935/shaojf/stampede2/housekeeping_genes/both.pc.and.nc.genes/PSYCHIC/ mygoogle:hkg_tsg/both.pc.and.nc.genes/PSYCHIC/
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
for f in *.bed
do
    pre=`echo $f | sed 's/.bed//'`
    perl $myperl <(sed 's/:/\t/' $f) <(tail -n +2 hkg.tsg.srtbyPCA.class | sed 's/|/\t/') 3 1 | cut -f 1-6 | sed 's/\t/|/' > hkg.tsg.srtbyPCA.$pre
done

for f in *.bed
do
    pre=`echo $f | sed 's/.bed//'`
    # awk -F"\t" '$3!="/"{print $2"\t"$3":"$4"-"$5}' hkg.tsg.srtbyPCA.$pre
    awk -F"\t" '$2~/hkg/ && $3!="/"{print $1}' hkg.tsg.srtbyPCA.$pre | sort | uniq | wc -l
done > hkg.stats

for f in *.bed
do
    pre=`echo $f | sed 's/.bed//'`
    cut -f 4 $f | cut -f 1 -d":" | sort | uniq | wc -l
done > all.stats

243
366
1089
1044
766
906
957
409
494

2193
2906
8373
8982
6889
8105
7893
4359
6373