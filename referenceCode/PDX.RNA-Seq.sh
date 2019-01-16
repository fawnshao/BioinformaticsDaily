#!/bin/bash

### for STAR
hg19file=/data/shaojf/myReference/star_index/hg19.fa
hg19indexdir=/data/shaojf/myReference/star_index/hg19.star
genomesize=/data/shaojf/myReference/star_index/hg19.chrom.sizes

workDir=/data/shaojf/BCM.PDX.RNA-seq/STAR.map
fastqDir=$workDir/RawFastq
genomeDir=$workDir/genome.map
mkdir -p $genomeDir
cd $workDir

for f in $fastqDir/*_R1_001.fastq.gz
do
	date
	pre=`echo $f | awk -F"/" '{print $NF}' | sed 's/_R1_001.fastq.gz//'`
	echo $pre

	STAR --runMode alignReads --runThreadN 48 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn $fastqDir/${pre}_R1_001.fastq.gz $fastqDir/${pre}_R2_001.fastq.gz \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFilterScoreMinOverLread  0.5 \
	--outFilterMatchNminOverLread 0.5 \
	--outFileNamePrefix $genomeDir/${pre}.toGenome. \
	--readFilesCommand zcat \
	--outSAMattributes All \
	--outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:$pre \
	--alignEndsType Local

	samtools sort -@ 12 -o $genomeDir/${pre}.toGenome.srt.bam $genomeDir/${pre}.toGenome.Aligned.out.bam &
done
wait
for f in $genomeDir/*.toGenome.srt.bam
do
	samtools index -@ 36 $f
done



