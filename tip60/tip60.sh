#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/TIP60.KAT5/ mygoogle:TIP60.KAT5/
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

# cd /home1/04935/shaojf/stampede2/TIP60.KAT5/PANCAN
gunzip -c cohort.TCGA.Pan-Cancer.PANCAN/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz | grep -w -e KAT5 -e EP400 -e TCGA > Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.KAT5.EP400
sed -n '2p' Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.KAT5.EP400 | head_line | awk '$2<-1{print $1}' | tr "\n" "," | sed 's/,$//' | xargs -n 1 -I mycol cut -f mycol <(head -1 Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.KAT5.EP400) | tr "\t" "\n" > KAT5.del.sample
sed -n '2p' Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.KAT5.EP400 | head_line | awk '$2>1 && $2!="KAT5"{print $1}' | tr "\n" "," | sed 's/,$//' | xargs -n 1 -I mycol cut -f mycol <(head -1 Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.KAT5.EP400) | tr "\t" "\n" > KAT5.amp.sample
sed -n '2p' Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.KAT5.EP400 | head_line | awk '$2==0{print $1}' | tr "\n" "," | sed 's/,$//' | xargs -n 1 -I mycol cut -f mycol <(head -1 Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.KAT5.EP400) | tr "\t" "\n" > KAT5.nocnv.sample

grep -wf KAT5.del.sample <(head_line tcga_RSEM_gene_tpm) | awk '{print $1}' | tr "\n" "," | sed 's/,$//' | xargs -n 1 -I mycol cut -f 1,mycol tcga_RSEM_gene_tpm > KAT5.del.tpm
grep -wf KAT5.amp.sample <(head_line tcga_RSEM_gene_tpm) | awk '{print $1}' | tr "\n" "," | sed 's/,$//' | xargs -n 1 -I mycol cut -f 1,mycol tcga_RSEM_gene_tpm > KAT5.amp.tpm
# grep -wf KAT5.nocnv.sample <(head_line tcga_RSEM_gene_tpm) | awk '{print $1}' | tr "\n" "," | sed 's/,$//' | xargs -n 1 -I mycol cut -f 1,mycol tcga_RSEM_gene_tpm > KAT5.nocnv.tpm
paste KAT5.del.tpm <(cut -f 2- KAT5.amp.tpm) > KAT5.del.amp.tpm
echo "sample Type" | tr " " "\t" > KAT5.del.amp.samples
head_line KAT5.del.tpm | tail -n +2 | awk '{print $2"\tDel"}' >> KAT5.del.amp.samples
head_line KAT5.amp.tpm | tail -n +2 | awk '{print $2"\tAmp"}' >> KAT5.del.amp.samples

gunzip -c gencode.v23.annotation.gff3.gz | awk '$3=="gene"{print $9}' | cut -f2-3 -d";" | sed 's/gene_id=//;s/gene_type=//;s/;/\t/' > gencode.v23.gene_type
sort -k2 gencode.v23.gene_type > gencode.v23.gene_type.srt
head -1 KAT5.del.amp.tpm > KAT5.del.amp.tpm.srt
perl $myperl KAT5.del.amp.tpm gencode.v23.gene_type.srt 0 0 | cut -f 3- >> KAT5.del.amp.tpm.srt

perl $myperl <(sort tcga_RSEM_gene_tpm.samples.sim | uniq) KAT5.del.amp.samples 0 0 | cut -f 1-2,4 > KAT5.del.amp.samples.2
perl $myperl <(sort tcga_RSEM_gene_tpm.samples.sim | uniq) KAT5.nocnv.sample 0 0 | grep "Head and Neck region|Primary Tumor" | cut -f 1,3 > Skin.KAT5.nocnv.sample.2

grep -e "Skin" -e "sample" KAT5.del.amp.samples.2 > Skin.KAT5.del.amp.samples.2
grep -n "Skin" KAT5.del.amp.samples.2 | cut -f 1 -d":" | tr "\n" "," | sed 's/,$//' | xargs -n 1 -I mycol cut -f 1,mycol KAT5.del.amp.tpm.srt > Skin.KAT5.del.amp.tpm.srt

# Head and Neck region|Primary Tumor

######
# library(data.table)
# library(pheatmap)
# # args <- c("KAT5.del.amp.tpm.srt", "KAT5.del.amp.samples.2", "gencode.v23.gene_type.srt")
# args <- c("Skin.KAT5.del.amp.tpm.srt", "Skin.KAT5.del.amp.samples.2", "gencode.v23.gene_type.srt")
# input <- fread(args[1], sep = "\t", header = T)
# scores <- data.matrix(input[,-1])
# colnames(scores) <- 1:ncol(scores)
# rownames(scores) <- 1:nrow(scores) # 1:nrow(scores) # as.matrix(input[,1])
# class <- fread(args[2], sep = "\t", header = T)
# annosC <- class[,-1]
# rownames(annosC) <- colnames(scores)
# colnames(annosC) <- c("Type", "Tissue")
# types <- fread(args[3], sep = "\t", header = F)
# annosR <- types[,2]
# rownames(annosR) <- rownames(scores)

# colors <- colorRampPalette(c("blue", "white", "red"))(100)
# scores[!is.na(scores) & scores < 0 ] <- 0

# del.mean <- apply(scores[,class$Type=="Del"], 1, mean)
# amp.mean <- apply(scores[,class$Type=="Amp"], 1, mean)
# logfc <- amp.mean - del.mean
# scores.diff <- scores[abs(logfc) > 1,]
# annosR.diff <- annosR[abs(logfc) > 1,]
# rownames(annosR.diff) <- rownames(scores.diff)
# scores.diff[!is.na(scores.diff) & scores.diff > 15] <- 15
# scores.diff.o <- scores.diff[,order(annosC$Type, annosC$Tissue)]
# annosC.o <- annosC[order(annosC$Type, annosC$Tissue),]
# rownames(annosC.o) <- colnames(scores.diff.o)
# png(filename = paste(args[1], "diff.pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(scores.diff.o, scale = "none", annotation_col = annosC.o, annotation_row = annosR.diff, 
# 	show_rownames = F, show_colnames = F, color = colors, 
# 	cluster_cols = F, cluster_rows = T)
# dev.off()

# annosR.diff.nc <- annosR.diff[annosR.diff$V2!="protein_coding",]
# rownames(annosR.diff.nc) <- rownames(scores.diff.o[annosR.diff$V2!="protein_coding",])
# png(filename = paste(args[1], "diff.ncRNA.pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(scores.diff.o[annosR.diff$V2!="protein_coding",], scale = "none", annotation_col = annosC.o, annotation_row = annosR.diff.nc, 
# 	show_rownames = F, show_colnames = F, color = colors, 
# 	cluster_cols = F, cluster_rows = T)
# dev.off()

# scores[!is.na(scores) & scores > 15] <- 15
# png(filename = paste(args[1], "pheatmap.png", sep = "."), width = 1500, height = 1200)
# myplot <- pheatmap(scores[,order(annosC$Type, annosC$Tissue)], scale = "none", annotation_col = annosC.o, annotation_row = annosR, 
# 	show_rownames = F, show_colnames = F, color = colors, 
# 	cluster_cols = F, cluster_rows = F)
# dev.off()
######

############# for DSB
# cd /home1/04935/shaojf/stampede2/TIP60.KAT5/sv
for f in structural_somatic_mutation.*.tsv.gz
do
	pre=`echo $f | cut -f 2 -d"."`
	gunzip -c $f | cut -f 2-3,7,12-13,17-18 > DSB.$pre.txt
done
for f in `ls DSB.*.txt`
do
	pre=`echo $f | cut -f 2 -d"."`
	ln -s ../exprs/KAT5.exp_*.$pre.tsv.gz.txt.specimen .
done


# head_line human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt | awk -vOFS="\t" '{print $1,$2}' | perl $myperl Human.sample_name2library_id.txt /dev/stdin 1 1 | cut -f 1-3 > fantom.enhancer.tissues.txt
# ll DSB.*.txt | cut -f2 -d"." | grep -wf /dev/stdin ../ICGC.id 
# BRCA-FR Breast Cancer - FR
# CLLE-ES Chronic Lymphocytic Leukemia - ES
# LIRI-JP Liver Cancer - RIKEN, JP
# MALY-DE Malignant Lymphoma - DE
# OV-AU Ovarian Cancer - AU
# PACA-AU Pancreatic Cancer - AU
# PACA-CA Pancreatic Cancer - CA
# PAEN-AU Pancreatic Cancer Endocrine neoplasms - AU
# PRAD-CA Prostate Adenocarcinoma - CA
# PRAD-FR Prostate Cancer - Adenocarcinoma - FR
for f in `ls DSB.*.txt`
do
	pre=`echo $f | cut -f 2 -d"."`
	tissues=`grep -w $pre ICGC.tissues | awk '{print $2}'`
	grep -i -e $tissues fantom.enhancer.tissues.txt | cut -f 1 | tr "\n" "," | sed 's/,$//' | xargs -I mycol cut -f 1,mycol human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt | awk '{sum=0;for(i=2;i<=NF;i++){sum+=$i}if(sum>0){print $1"\t"$1}}' | sed 's/:/\t/;s/-/\t/' > $pre.active.enhancer.bed
done

for f in `ls DSB.*.txt`
do
	pre=`echo $f | cut -f 2 -d"."`
	bedtools intersect -wo -a <(tail -n +2 DSB.$pre.txt | awk -F"\t" -vOFS="\t" '{a=$5-5000;b=$7-5000;if(a<0){a=0}if(b<0){b=0}print "chr"$4,a,$5+5000,$2"\nchr"$6,b,$7+5000,$2}' | uniq) -b <(awk -vOFS="\t" '{print $1,$2-5000,$3+5000,$4}' $pre.active.enhancer.bed) > DSB.enhancer.$pre
done

for f in `ls DSB.*.txt`
do
	pre=`echo $f | cut -f 2 -d"."`
	cut -f 1-4 DSB.enhancer.$pre | sort | uniq | cut -f 4 | sort | uniq	-c | awk '{print $2"\t"$1}' > DSB.enhancer.$pre.stats 
done

for f in `ls DSB.*.txt`
do
	pre=`echo $f | cut -f 2 -d"."`
	if [ -f KAT5.exp_array.$pre.tsv.gz.txt.specimen ]
	then
		perl $myperl DSB.enhancer.$pre.stats KAT5.exp_array.$pre.tsv.gz.txt.specimen 0 1 | grep -v "/" | cut -f 1-5,7 > exp_array.DBS.$pre.tsv
	fi
	if [ -f KAT5.exp_seq.$pre.tsv.gz.txt.specimen ]
	then
		perl $myperl DSB.enhancer.$pre.stats KAT5.exp_seq.$pre.tsv.gz.txt.specimen 0 1 | grep -v "/" | cut -f 1-5,7 > exp_seq.DBS.$pre.tsv
	fi
done

rscatterplot=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/exp.mut.scatter.plot.R

Rscript $rscatterplot $f



for f in copy_number_somatic_mutation.*.tsv.gz
do
	pre=`echo $f | sed 's/copy_number_somatic_mutation.//;s/.tsv.gz//'`
	if [ ! -f KAT5.exp_array.$pre.tsv.gz.txt.specimen ] && [ ! -f KAT5.exp_seq.$pre.tsv.gz.txt.specimen ]
	then
		rm $f
	else
		gunzip -c $f | cut -f 2-3,8,20 | tail -n +2 | sed 's/ /./g' | sort | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$5,$1}' >> cnv.byspeciman.stats
	fi
done


