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
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 8:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job
runname="NS14-Joo-gRNA-Era"
workdir=/home1/04935/shaojf/scratch/NS14-Joo-gRNA-Era
fastqDir=$workdir/fastqs
fastqcDir=$workdir/fastqcRes
cleanDir=$workdir/cleanfastq
visdir=$workdir/bigwigs
Threads=270
mkdir $fastqcDir
mkdir $cleanDir
HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

fastqc --threads $Threads --outdir $fastqcDir $fastqDir/*.fastq.gz
# AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT NNNNNNNN GCTTTAT ATATCTTGTGGAAAGGA CGAAACACC
# CGAAACACCG – 20 nt gRNA - gttttagagc
# AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNGCTTTATATATCTTGTGGAAAGGACGAAACACCG
# -g NNNNNNNNGCTTTATATATCTTGTGGAAAGGACGAAACACCG \
	# -b GTTTTAGAGCTAGGCCA \
module load python3/3.6.3
for f in $fastqDir/*.fastq.gz
do
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	~/.local/bin/cutadapt --nextseq-trim=20 \
	--match-read-wildcards \
	-a NNNNNNNNGCTTTATATATCTTGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGGCCA \
	-m 20 -j $Threads \
	-o $cleanDir/$i.clean.fastq $f 1>$i.cutadapt.log 2>&1 &
done
wait

for f in $fastqDir/*.fastq.gz
do
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	echo $i
	# bowtie2 --local --very-sensitive -p $Threads -x $HG19 -U $cleanDir/$i.clean.fastq -S $i.sam
	bowtie2 -L 18 -p $Threads -x $HG19 -U $cleanDir/$i.clean.fastq -S $i.sam
	samtools view -1 -q 10 --threads 30 $i.sam | samtools sort --threads 30 > $i.sorted.bam
	samtools markdup -r --threads $Threads $i.sorted.bam $i.rmdup.bam
	samtools index $i.rmdup.bam
	bamCoverage -b $i.rmdup.bam -o $i.bigWig
	# makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
done
# makeMultiWigHub.pl $runname hg19 -url $visdir -webdir $visdir -d *.mTD

# TTTTCNTGGCTTTATATATCTTGTGGAAAGGA
# TTTTCNTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGAGACGGCATCTTCACCGTAC GTTTTAGAGCTAGGCCA
# TTTTCNTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGAGACGGCATCTTCACCGTAC GTTTTAGAGC TAGGCCA
# CGCGCCATGCTTTATATATCTTGTGGAAAGGACGAAACACCGGAGACGGCATCTTCCCCGTAC GTTTAAGAGCTAGGCCA


# @NS500669:263:HVG5CBGX5:1:11101:7361:1044 1:N:0:TCGAAG
# 					                      GAGACGGCATCTTCACCGTACGTTTTAGAGCTAGGCC
#                                           GAGACGGCATCTTCACCGTACGTTTTAGAGCTAGGCC
#                                           GAGACGGCATCTTCACCGTAC
# TTTTCNTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGAGACGGCATCTTCACCGTACGTTTTAGAGCTAGGCCA
# GAGACGGCATCTTCACCGTAC


# TTTTCNTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGAGACGGCATCTTCACCGTACGTTTTAGAGCTAGGCCA
#                                           GAGACGGCATCTTCACCGTAC
