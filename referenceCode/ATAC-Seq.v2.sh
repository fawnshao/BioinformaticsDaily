#!/bin/sh
#### additional analysis for ATAC-Seq
runname="NS13-Xiaoyu-NSC-T21"
workdir=/home1/04935/shaojf/scratch/NS13_xiaoyu_ATAC
FastqDir=$workdir/FastqsDir
FastqcDir=$workdir/FastqcDir
TrimadapterDir=$workdir/TrimmedFastqsDir
VisDir=$workdir/bigwigs
Threads=200

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

mkdir $FastqcDir
mkdir $TrimadapterDir

fastqc --threads $Threads --outdir $FastqcDir $FastqDir/*.fastq.gz &

for f in $FastqDir/*.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`

	cutadapt --nextseq-trim=20 \
	-a tag1=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
	-a tag1rev=CTGTCTCTTATACACATCTGACGCTGCCGACGA \
	-a index1=GTCTCGTGGGCTCGG \
	-a index1rev=CCGAGCCCACGAGAC \
	-a tag2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
	-a tag2rev=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
	-a index2=TCGTCGGCAGCGTC \
	-a index2rev=GACGCTGCCGACGA\
	-m 18 \
	-o $TrimadapterDir/$i.clean.fastq.gz $f &
done
wait
###### my pipeline
for f in $FastqDir/*.fastq.gz
do
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	echo $i

	bowtie2 --local --very-sensitive -p $Threads -x $HG19 -U $TrimadapterDir/$i.clean.fastq.gz -S $i.sam
	samtools view -1 -q 10 --threads $Threads $i.sam | samtools sort --threads 20 > $i.sorted.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
	findPeaks $i.mTD -style factor -o auto &
	samtools markdup -r -s -O BAM -@ $Threads $i.sorted.bam $i.rmdup.bam
	macs2 callpeak --verbose 3 -t $i.rmdup.bam -f BAM -g hs -n $i --keep-dup 1 --nomodel --shift -100 --extsize 200 --bdg --SPMR 1> $i.log 2>&1 &
done
makeMultiWigHub.pl $runname hg19 -url $VisDir -webdir $VisDir -d *.mTD
wait
