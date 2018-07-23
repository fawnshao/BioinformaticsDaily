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
#### additional analysis for ATAC-Seq
runname="CFP1.ChIP"
workdir=/home1/04935/shaojf/scratch/CFP1
FastqDir=$workdir/FastqsDir
# FastqcDir=$workdir/FastqcDir
VisDir=$workdir/bigwigs
Threads=200
pairs=SRR2932629:SRR2932630

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

# mkdir $FastqcDir
mkdir $VisDir

# fastqc --threads $Threads --outdir $FastqcDir $FastqDir/*.fastq.gz &

###### my pipeline
for f in $FastqDir/*_1.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/_1.fastq.gz//'`
	echo $i

	bowtie2 --local --very-sensitive -p $Threads -x $HG19 \
	-1 $FastqDir/${i}_1.fastq.gz -2 $FastqDir/${i}_2.fastq.gz -S $i.sam
	samtools view -1 -q 10 -bo $i.q10.bam --threads $Threads $i.sam
	samtools sort --threads 20 $i.q10.bam -o $i.sorted.bam
	samtools fixmate --threads $Threads -m $i.q10.bam $i.fixmate.bam
	samtools sort --threads 20 $i.fixmate.bam -o $i.fixmate.srt.bam
	samtools markdup -r -s -O BAM -@ $Threads $i.fixmate.srt.bam $i.rmdup.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.rmdup.bam
	findPeaks $i.mTD -style factor -o auto &
done
t=`echo $pairs | awk -F":" '{print $1}'`
c=`echo $pairs | awk -F":" '{print $2}'`
macs2 callpeak --verbose 3 -t $t.rmdup.bam -c $c.rmdup.bam -f BAMPE -g hs -n CFP1.$t 1> CFP1.$t.bampe.log 2>&1 &
makeMultiWigHub.pl $runname hg19 -url $VisDir -webdir $VisDir -d *.mTD
wait
