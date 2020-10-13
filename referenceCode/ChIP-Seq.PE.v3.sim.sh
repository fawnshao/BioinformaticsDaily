#!/bin/bash
#SBATCH -J myjob           # Job name
#SBATCH -o myjob.out       # Name of stdout output file
#SBATCH -e myjob.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 8:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

WorkDir=/data/shaojf/LiLab.Seq/NS36
FastqDir=$WorkDir/FastqDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
hubname="NS36.ChIP"
runname="NS36.ChIP"
fastqpre="NS36-ChIP-Xiaoyu-ChIP"
fastqtail="_R1_001.fastq.gz"

mkdir -p $BamDir
mkdir -p $BwDir
Threads=50

HG19=/data/shaojf/myReference/bowtie2-index/hg19
genomesize=/data/shaojf/myReference/star_index/hg19.chrom.sizes
blackListFileName=/data/shaojf/myReference/hg19.blacklist.ENCFF001TDO.merged.bed

# FastqcDir=$WorkDir/FastqcDir
# mkdir $FastqcDir
# fastqc --threads $Threads --outdir $FastqcDir $FastqDir/*.fastq.gz &
###### my pipeline
for f in $FastqDir/$fastqpre*$fastqtail
do
	echo $f
	i=`basename $f | sed "s?$fastqtail??"`
	echo $i

	bowtie2 --local --very-sensitive -p $Threads -x $HG19 \
		-1 $FastqDir/${i}_R1_001.fastq.gz \
		-2 $FastqDir/${i}_R2_001.fastq.gz -S $BamDir/$i.bt2.sam
	samtools view -1 -q 10 -bo $BamDir/$i.q10.bam --threads $Threads $BamDir/$i.bt2.sam
	samtools sort --threads 10 $BamDir/$i.q10.bam -o $BamDir/$i.sorted.bam

	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 -fragLength pe -format sam $BamDir/$i.sorted.bam
	findPeaks $i.mTD -style factor -o auto &

	samtools fixmate --threads $Threads -m $BamDir/$i.q10.bam $BamDir/$i.fixmate.bam
	samtools sort --threads 20 $BamDir/$i.fixmate.bam -o $BamDir/$i.fixmate.srt.bam
	samtools markdup -r -s -O BAM -@ $Threads $BamDir/$i.fixmate.srt.bam $BamDir/$i.rmdup.bam
	samtools index $BamDir/$i.rmdup.bam
	macs2 callpeak --verbose 3 -t $BamDir/$i.rmdup.bam -f BAMPE -g hs -n $i --keep-dup all 1> $i.log 2>&1 &
	bamCoverage -b $BamDir/$i.rmdup.bam -o $i.bamCoverage.pe.CPM.bw -of bigwig \
		-bs 1 -p 20 --normalizeUsing CPM -e 150 --blackListFileName $blackListFileName &
done
makeMultiWigHub.pl $runname hg19 -url $BwDir -webdir $BwDir -d $fastqpre*.mTD
wait

# for f in NS*.q10.bam
# do
# 	echo $f >> samtools.stats.txt
# 	samtools flagstat -@ 60 $f >> samtools.stats.txt
# done
# for f in  NS*.rmdup.bam
# do
# 	echo $f >> samtools.stats.txt
# 	samtools flagstat -@ 60 $f >> samtools.stats.txt
# done


