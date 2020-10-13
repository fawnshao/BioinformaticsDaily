#!/bin/bash
module load samtools-1.9
starindexdir=/public/workspace/shaojf/myref/STAR.hg38
WorkDir=/public/workspace/shaojf/Wangcc
FastqDir=$WorkDir/FastqDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
Threads=32

cd $WorkDir
mkdir -p $BamDir
mkdir -p $BwDir
# mkdir -p $CleanDir

blacklist=/public/workspace/shaojf/myref/UCSC/hg38/enschr.hg38.blacklist.ENCFF356LFX.bed
for f in $FastqDir/*.1.fq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's?.1.fq.gz??'`
	STAR --runMode alignReads --runThreadN $Threads \
		--genomeDir $starindexdir --genomeLoad LoadAndRemove \
		--readFilesIn $FastqDir/${i}.1.fq.gz $FastqDir/${i}.2.fq.gz \
		--outSAMunmapped Within \
		--readFilesCommand zcat \
		--outFilterMultimapNmax 1 \
		--outFilterMultimapScoreRange 1 \
		--outFilterScoreMinOverLread  0.6 \
		--outFilterMatchNminOverLread 0.6 \
		--outFileNamePrefix $BamDir/$i. \
		--outSAMattributes All --outSAMtype BAM Unsorted \
		--outFilterType BySJout --outReadsUnmapped Fastx \
		--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
		--alignEndsType Local

	samtools view -q 30 -@ $Threads -bo $BamDir/$i.q30.bam $BamDir/$i.Aligned.out.bam
	samtools sort -@ $Threads -o $BamDir/$i.q30.srt.bam $BamDir/$i.q30.bam
	samtools index $BamDir/$i.q30.srt.bam
	bamCoverage -b $BamDir/$i.q30.srt.bam -o $BwDir/$i.bw \
		-bs 1 --blackListFileName $blacklist -p $Threads \
		--effectiveGenomeSize 2913022398 --normalizeUsing CPM &
done
wait

