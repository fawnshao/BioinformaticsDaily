#!/bin/bash
# hg19indexdir=/data/shaojf/myReference/star_index/hg19.star
hg19indexdir=/data/shaojf/myReference/star_index/STAR-hg19
WorkDir=/data/shaojf/LiLab.Seq/NS35
FastqDir=$WorkDir/FastqDir
CleanDir=$WorkDir/CleanDir
BamDir=$WorkDir/BamDir
# BwDir=$WorkDir/BwDir
# hubname="MCF7-shNUP93-PRO"
# runname="MCF7-shNUP93-PRO"
# Threads=272

cd $WorkDir
mkdir -p $BamDir
# mkdir -p $BwDir
# mkdir -p $CleanDir

# --readFilesCommand zcat \
for f in $FastqDir/*_R1_001.fastq
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's?_R1_001.fastq??'`
	STAR --runMode alignReads --runThreadN 60 \
		--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
		--readFilesIn $FastqDir/${i}_R1_001.fastq \
			$FastqDir/${i}_R2_001.fastq \
		--outSAMunmapped Within \
		--outFilterMultimapNmax 1 \
		--outFilterMultimapScoreRange 1 \
		--outFilterScoreMinOverLread  0.3 \
		--outFilterMatchNminOverLread 0.3 \
		--outFileNamePrefix $BamDir/$i. \
		--outSAMattributes All --outSAMtype BAM Unsorted \
		--outFilterType BySJout --outReadsUnmapped Fastx \
		--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
		--alignEndsType Local

	# bowtie2 -5 3 -3 1 -p $Threads -x $HG19 \
	# 	-1 $CleanDir/$i.R1.umiprocessed.fastq.gz \
	# 	-2 $CleanDir/$i.R2.umiprocessed.fastq.gz -S $BamDir/$i.bt2.sam
	# samtools view -q 10 -@ 272 -bo $BamDir/$i.q10.bam $BamDir/$i.Aligned.out.bam
	# samtools sort -@ 68 -o $BamDir/$i.q10.srt.bam $BamDir/$i.q10.bam
	# samtools index $BamDir/$i.q10.srt.bam
	# umi_tools dedup -I $BamDir/$i.q10.srt.bam -S $BamDir/$i.q10.dedup.bam -L $BamDir/$i.umi_tools.log \
	# 	--paired --output-stats $i.umi_tools.stats
	# makeTagDirectory $i.mTD -sspe -flip $BamDir/$i.q10.dedup.bam &
done

