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

# Launch serial code...
runname="NS15-Xiaoyu-MCF7-TT"
workdir=/home1/04935/shaojf/scratch/NS15_xiaoyu_TTseq
fastqDir=$workdir/raw.data
cleanDir=$workdir/clean.data
fastqcDir=$workdir/fastqc.res
visdir=$workdir/bigwigs
Threads=272

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

mkdir -p $fastqcDir
mkdir -p $cleanDir
mkdir -p $visdir
###### my pipeline
fastqc --threads $Threads --outdir $fastqcDir $fastqDir/*.fastq.gz &
module load python3/3.6.3
for f in $fastqDir/*.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	
	~/.local/bin/cutadapt --nextseq-trim=20 \
	-a tag1=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
	-a tag1rev=CTGTCTCTTATACACATCTGACGCTGCCGACGA \
	-a index1=GTCTCGTGGGCTCGG \
	-a index1rev=CCGAGCCCACGAGAC \
	-a tag2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
	-a tag2rev=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
	-a index2=TCGTCGGCAGCGTC \
	-a index2rev=GACGCTGCCGACGA\
	-m 18 -j $Threads \
	-o $cleanDir/$i.clean.fastq $f &
done
wait

for f in $fastqDir/*.fastq.gz
do
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	echo $i

	bowtie2 --local --very-sensitive -p $Threads -x $HG19 -U $cleanDir/$i.clean.fastq -S $i.sam
	samtools view -1 -q 10 --threads $Threads $i.sam | samtools sort --threads 48 > $i.sorted.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 3 -flip $i.sorted.bam &
done
wait
makeMultiWigHub.pl $runname hg19 -url $visdir -webdir $visdir -d *.mTD -force -strand &
analyzeRepeats.pl rna hg19 -rpkm -strand + -condenseGenes -d *.mTD > $runname.rpkm.txt &
wait
