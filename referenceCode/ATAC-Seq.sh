#!/bin/bash
# this code will generate .sorted.bam files first, and then generate Homer tag directory and UCSC file as track hubs
workdir=/home1/04935/shaojf/scratch/NUP_related_others/ATAC-seq
fastq_dir=$workdir/fastqs
fastqc_dir=$workdir/fastqc_res
trimadapterdir=$workdir/fastqs_trimmed
mappingtogenome=$workdir/reads2genome
bw_dir=$workdir/bigwigs
dhspeaks=$workdir/dhspeaks

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
visdir=/home1/04935/shaojf/stampede2/UCSCvis

fastqc --threads 68 --outdir $fastqc_dir $fastq_dir/*.fastq &
echo    ***____running fastqc in background____***

# Illumina Nextera Adapters
# Nextera Transposase Adapters
# (Used for Nextera tagmentation)
# Read 1
# 5’ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
# Read 2
# 5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
# Nextera Index Kit – PCR Primers
# Index 1 Read
# 5’ CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG
# Index 2 Read
# 5’ AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC

for f in $fastq_dir/*.fastq
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq//'`

	cutadapt --nextseq-trim=20 \
	-a tag1=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
	-a tag1rev=CTGTCTCTTATACACATCTGACGCTGCCGACGA \
	-a index1=GTCTCGTGGGCTCGG \
	-a index1rev=CCGAGCCCACGAGAC \
	-a tag2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
	-a tag2rev=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
	-a index2=TCGTCGGCAGCGTC \
	-a index2rev=GACGCTGCCGACGA\
	-o $trimadapterdir/$i.clean.fastq $f &
	echo    ***____step0: running cutadapt successful____***
done
wait

for f in $fastq_dir/*.fastq
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq//'`

	bowtie2 --local --very-sensitive -p 68 -x $HG19 -U $trimadapterdir/$i.clean.fastq | samtools view -1 -q 10 --threads 68 | samtools sort --threads 68 > $i.sorted.bam
	echo    ***____step1: bowtie2 successful____***

	nice makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
	echo    ***____step2: home makeTagDirectory successful____***
	echo    ***--- these samples were normed to 10million total reads____***

	nice makeBigWig.pl $i.mTD/ hg19 -url $visdir -webdir $visdir
	echo    ***____step3: minus trand UCSC bedgraph.gz file of $i.mTD folder completed____***

	nice findPeaks $i.mTD -style factor -o auto
	echo    ***____step4: findPeaks for $i.mTD folder completed____***
	echo    ++++++++++++++++++++++++ Yahoo - samples all finished  +++++++++++++++++++++++++++++++++++

	samtools rmdup $i.sorted.bam $i.sorted.rmdup.bam
	macs2 callpeak --verbose 3 -t $i.sorted.rmdup.bam -f BAM -g hs -n $i --call-summits
done

makeMultiWigHub.pl NUP_ChIP hg19 -url $visdir -webdir $visdir -d *.mTD
