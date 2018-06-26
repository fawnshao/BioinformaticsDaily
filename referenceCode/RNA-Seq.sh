#!/bin/bash

### for STAR
hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
hg19indexdir=/home1/04935/shaojf/scratch/star_index/hg19.star
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
erccfile=/home1/04935/shaojf/scratch/star_index/ERCC92.fa
erccindexdir=/home1/04935/shaojf/scratch/star_index/ERCC92.star
# mkdir $erccindexdir
# STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $erccindexdir --genomeFastaFiles $erccfile

workDir=/home1/04935/shaojf/scratch/YY1.mut.RNAseq
fastqDir=$workDir/fastqDir
erccDIR=$workDir/ERCC.map
genomeDir=$workDir/genome.map
mkdir -p $erccDIR
mkdir -p $genomeDir
cd $workDir

# cut -f 1,9 SRP064870.runinfo.tsv | grep OTHER | cut -f 2 | xargs -n 1 prefetch
# for f in *.sra
# do
# 	fastq-dump --split-files $f &
# done
# while read line
# do
# 	srr=`echo $line | awk '{print $1}'`
# 	name=`echo $line | awk '{print $4"."$3}'`
# 	cat ${srr}_1.fastq >> ${name}_1.fastq &
# 	cat ${srr}_2.fastq >> ${name}_2.fastq &
# 	wait
# done < mysamples
# mv LC_* ../fastqDir/
date
for f in $fastqDir/*_1.fastq
do
	pre=`echo $f | awk -F"/" '{print $NF}' | sed 's/_1.fastq//'`
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $erccindexdir --genomeLoad LoadAndRemove \
	--readFilesIn $fastqDir/${pre}_1.fastq $fastqDir/${pre}_2.fastq \
	--outSAMunmapped None \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix $erccDIR/${pre}.toERCC. \
	--outSAMattributes All \
	--outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local

	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn $erccDIR/${pre}.toERCC.Unmapped.out.mate1 \
	$erccDIR/${pre}.toERCC.Unmapped.out.mate2 \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFilterScoreMinOverLread  0.5 \
	--outFilterMatchNminOverLread 0.5 \
	--outFileNamePrefix $genomeDir/${pre}.toGenome. \
	--outSAMattributes All \
	--outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local
done

