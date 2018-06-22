#!/bin/sh
#### additional analysis for ATAC-Seq
runname="NS13-Xiaoyu-NSC-T21"
workdir=/home1/04935/shaojf/scratch/NS13_xiaoyu_ATAC
fastq_dir=$workdir/fastqs
fastqc_dir=$workdir/fastqc_res
trimadapterdir=$workdir/fastqs_trimmed
visdir=$workdir/bigwigs

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

mkdir $fastqc_dir
mkdir $trimadapterdir

fastqc --threads 68 --outdir $fastqc_dir $fastq_dir/*.fastq.gz &

for f in $fastq_dir/*.fastq.gz
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
	--info-file=$i.cutadapt.metrics \
	-o $trimadapterdir/$i.clean.fastq.gz $f &
done
wait
###### my pipeline
for f in $fastq_dir/*.fastq.gz
do
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	echo $i

	bowtie2 --local --very-sensitive -p 68 -x $HG19 -U $trimadapterdir/$i.clean.fastq.gz -S $i.sam
	samtools view -1 -q 10 --threads 68 $i.sam | samtools sort --threads 68 > $i.sorted.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
	# makeBigWig.pl $i.mTD/ hg19 -url $visdir -webdir $visdir &
	findPeaks $i.mTD -style factor -o auto &
	# picard MarkDuplicates I=$i.sorted.bam O=marked_duplicates.$i.sorted.bam M=$i.marked_dup_metrics.txt
	samtools markdup -r -s -O BAM -@ 68 $i.sorted.bam $i.rmdup.bam
	macs2 callpeak --verbose 3 -t $i.rmdup.bam -f BAM -g hs -n $i --keep-dup 1 --nomodel --shift -100 --extsize 200 --bdg --SPMR 1> $i.log 2>&1 &
done
makeMultiWigHub.pl $runname hg19 -url $visdir -webdir $visdir -d *.mTD
wait
