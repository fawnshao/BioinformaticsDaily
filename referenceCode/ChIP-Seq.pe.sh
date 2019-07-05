#!/bin/bash
#----------------------------------------------------
# Sample SLURM job script
#   for TACC Stampede2 KNL nodes
#
#   *** Serial Job on Normal Queue ***
#
# Last revised: 27 Jun 2017
#
# Notes:
#
#   -- Copy/edit this script as desired.  Launch by executing
#      "sbatch knl.serial.slurm" on a Stampede2 login node.
#
#   -- Serial codes run on a single node (upper case N = 1).
#        A serial code ignores the value of lower case n,
#        but slurm needs a plausible value to schedule the job.
#
#   -- For a good way to run multiple serial executables at the
#        same time, execute "module load launcher" followed
#        by "module help launcher".

#----------------------------------------------------

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.out       # Name of stdout output file
#SBATCH -e myjob.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1             # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 10:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19

WorkDir=/home1/04935/shaojf/scratch/NUP_related_others/NS25-PROSeq_ChIPSeq/
FastqDir=$WorkDir/FastqDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
hubname="MCF7-NUP93-ChIP"
runname="MCF7-NUP93-ChIP"
Threads=272

cd $WorkDir
# mkdir -p $BamDir
# mkdir -p $BwDir

for f in $FastqDir/NS25-Xiaoyu-ChIP*_R1_001.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's?_R1_001.fastq.gz??'`
	bowtie2 --local --very-sensitive -p $Threads -x $HG19 \
		-1 $FastqDir/${i}_R1_001.fastq.gz \
		-2 $FastqDir/${i}_R1_001.fastq.gz -S $BamDir/$i.bt2.sam
	samtools view -1 -q 10 -bo $BamDir/$i.q10.bam --threads $Threads $BamDir/$i.bt2.sam
	samtools sort --threads 20 $BamDir/$i.q10.bam -o $BamDir/$i.sorted.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 -format sam $BamDir/$i.sorted.bam
	findPeaks $i.mTD -style factor -o auto &
	samtools fixmate --threads $Threads -m $BamDir/$i.q10.bam $BamDir/$i.fixmate.bam
	samtools sort --threads 20 $BamDir/$i.fixmate.bam -o $BamDir/$i.fixmate.srt.bam
	samtools markdup -r -s -O BAM -@ $Threads $BamDir/$i.fixmate.srt.bam $BamDir/$i.rmdup.bam
	macs2 callpeak --verbose 3 -t $BamDir/$i.rmdup.bam -f BAM -g hs -n $i \
		--keep-dup all 1> $i.log 2>&1 &
done
wait
makeMultiWigHub.pl $hubname hg19 -url $BwDir -webdir $BwDir -d *ChIP*.mTD -force &
wait
