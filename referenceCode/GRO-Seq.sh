#!/bin/sh
HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
visdir=/home1/04935/shaojf/stampede2/UCSCvis
url=http://sjf.dingding.biz/trackhubs
for i in MCF7-sh*.fastq.gz
do
	cutadapt --nextseq-trim=28 -o ${i%.f*}_trimed.fastq $i
	cutadapt -a 'A{100}' ${i%.f*}_trimed.fastq > ${i%.f*}_trimed2.fastq
	rm ${i%.f*}_trimed.fastq
	bowtie2 -5 3 -3 1 -p 60  -x $HG19 ${i%.f*}_trimed2.fastq | samtools view -1 | samtools sort > ${i%.f*}.sorted.bam
	makeTagDirectory ${i%.f*}/ -tbp 3 -fragLength 200 ${i%.f*}.sorted.bam
done
makeMultiWigHub.pl MCF7-shNUP53_GroSeq_E2 hg19 -url $url -webdir $visdir -d MCF7-sh*GRO-E2 -force -strand