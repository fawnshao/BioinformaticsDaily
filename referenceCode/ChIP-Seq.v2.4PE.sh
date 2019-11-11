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

runname="NS31.SSBP1.ChIP"
workdir=/data/shaojf/LiLab.Seq/NS31.SSBP1.ChIP
FastqDir=$workdir/FastqDir
FastqcDir=$workdir/FastqcDir
VisDir=$workdir/bigwigs
Threads=70

HG19=/data/shaojf/myReference/bowtie2-index/hg19
genomesize=/data/shaojf/myReference/star_index/hg19.chrom.sizes

mkdir $FastqcDir
mkdir $VisDir

fastqc --threads $Threads --outdir $FastqcDir $FastqDir/*.fastq.gz &
###### my pipeline
for f in $FastqDir/*_R1_001.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/_R1_001.fastq.gz//'`
	echo $i

	bowtie2 --local --very-sensitive -p $Threads -x $HG19 \
		-1 $FastqDir/${i}_R1_001.fastq.gz -2 $FastqDir/${i}_R2_001.fastq.gz -S $i.sam
	samtools view -1 -q 10 -bo $i.q10.bam --threads $Threads $i.sam
	samtools sort --threads 20 $i.q10.bam -o $i.sorted.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 -fragLength pe -format sam $i.sorted.bam
	findPeaks $i.mTD -style factor -o auto &
	samtools fixmate --threads $Threads -m $i.q10.bam $i.fixmate.bam
	samtools sort --threads 20 $i.fixmate.bam -o $i.fixmate.srt.bam
	samtools markdup -r -s -O BAM -@ $Threads $i.fixmate.srt.bam $i.rmdup.bam
	samtools index $i.rmdup.bam
	macs2 callpeak --verbose 3 -t $i.rmdup.bam -f BAMPE -g hs -n $i --keep-dup all 1> $i.log 2>&1 &
	bedtools bamtobed -bedpe -i $i.rmdup.bam > $i.rmdup.bedpe
	macs2 callpeak --verbose 3 -t $i.rmdup.bedpe -f BEDPE -g hs -n $i.bedpe --keep-dup all 1> $i.bedpe.log 2>&1 &
	bamCoverage -b $i.rmdup.bam -o $i.bamCoverage.pe.CPM.bw -of bigwig -bs 1 -p 20 --normalizeUsing CPM -e 150 &
done
makeMultiWigHub.pl $runname hg19 -url $VisDir -webdir $VisDir -d *.mTD
wait
