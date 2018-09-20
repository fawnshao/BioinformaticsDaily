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

WorkDir=/home1/04935/shaojf/scratch/NS17.GRO-seq
FastqDir=$WorkDir/FastqDir
CleanDir=$WorkDir/CleanDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
hubname="MCF7-KD-GRO"
runname="MCF7-KD-GRO"
Threads=272

cd $WorkDir
mkdir -p $BamDir
mkdir -p $BwDir
mkdir -p $CleanDir

fastqc --threads $Threads --outdir $FastqDir $FastqDir/*.fastq.gz &
module load python3/3.6.3
### --local 
for f in $FastqDir/*.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | cut -d"-" -f 1-5`
	
	~/.local/bin/cutadapt --nextseq-trim=20 -m 18 -j $Threads \
	-o $CleanDir/$i.clean.fastq $f
	bowtie2 -5 3 -3 1 -p $Threads -x $HG19 -U $CleanDir/$i.clean.fastq -S $BamDir/$i.bt2.sam
	samtools view -q 10 -@ 272 -bo $BamDir/$i.q10.bam $BamDir/$i.bt2.sam
	makeTagDirectory $i.mTD -tbp 3 -fragLength 200 $BamDir/$i.q10.bam &
done
wait
makeMultiWigHub.pl $hubname hg19 -url $BwDir -webdir $BwDir -d *.mTD -force -strand

analyzeRepeats.pl rna hg19 -rpkm -strand + -condenseGenes -d *.mTD > $runname.gene.rpkm.txt &
# awk -vOFS="\t" '{print $4,$1,$2,$3,$6,$5}' MCF-7.p300.55905_countedputative.eRNAs.siCTL.bed > MCF-7.p300.55905_countedputative.eRNAs.txt
analyzeRepeats.pl MCF-7.p300.55905_countedputative.eRNAs.txt hg19 -rpkm -strand + -d *.mTD > $runname.eRNA.rpkm.txt &
