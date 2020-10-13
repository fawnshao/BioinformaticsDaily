#!/bin/bash

WorkDir=/data/shaojf/LiLab.Seq/NS36
FastqDir=$WorkDir/FastqDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
hubname="NS36.ChIP"
runname="NS36.ChIP"
fastqpre="NS36-Feng-TT"
fastqtail="_R1_001.fastq.gz"

mkdir -p $BamDir
mkdir -p $BwDir
Threads=50

# new star
HG19=/data/shaojf/myReference/star_index/STAR-hg19
genomesize=/data/shaojf/myReference/star_index/hg19.chrom.sizes
blackListFileName=/data/shaojf/myReference/hg19.blacklist.ENCFF001TDO.merged.bed

# FastqcDir=$WorkDir/FastqcDir
# mkdir $FastqcDir
# fastqc --threads $Threads --outdir $FastqcDir $FastqDir/*.fastq.gz &
###### my pipeline
for f in $FastqDir/$fastqpre*$fastqtail
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed "s?$fastqtail??"`
	echo $i

	STAR --runMode alignReads --runThreadN $Threads \
		--genomeDir $HG19 \
		--readFilesIn $FastqDir/${i}_R1_001.fastq.gz \
			$FastqDir/${i}_R2_001.fastq.gz \
		--outSAMunmapped Within \
		--outFilterMultimapNmax 1 \
		--outFilterMultimapScoreRange 1 \
		--outFilterScoreMinOverLread  0.5 \
		--outFilterMatchNminOverLread 0.5 \
		--outFileNamePrefix $BamDir/$i. \
		--outSAMattributes All --outSAMtype BAM Unsorted \
		--outFilterType BySJout --outReadsUnmapped Fastx \
		--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
		--alignEndsType Local --readFilesCommand zcat
	samtools view -1 -q 10 -bo $BamDir/$i.q10.bam --threads $Threads $BamDir/$i.Aligned.out.bam
	samtools sort --threads 10 $BamDir/$i.q10.bam -o $BamDir/$i.sorted.bam
	samtools index $BamDir/$i.sorted.bam

	bamCoverage -b $BamDir/$i.sorted.bam -o CPM.$i.plus.bw -of bigwig -bs 1 -p 30 \
		--blackListFileName $blackListFileName --filterRNAstrand forward --extendReads 150 \
		--effectiveGenomeSize 2864785220 --normalizeUsing CPM --maxFragmentLength 500 &
	bamCoverage -b $BamDir/$i.sorted.bam -o CPM.$i.minus.bw -of bigwig -bs 1 -p 30 \
		--blackListFileName $blackListFileName --filterRNAstrand reverse --extendReads 150 \
		--effectiveGenomeSize 2864785220 --normalizeUsing CPM --maxFragmentLength 500 &
done
wait
