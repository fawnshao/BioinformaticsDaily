#!/bin/sh
runname="NS13-Xiaoyu-MCF7-TT.flip"
workdir=/home1/04935/shaojf/scratch/NS13_xiaoyu_TTseq
fastqDir=$workdir/raw.data
fastqcDir=$workdir/fastqc.res
visdir=$workdir/bigwigs

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

###### my pipeline
for f in $fastqDir/*.fastq.gz
do
        i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
        echo $i

        bowtie2 --local --very-sensitive -p 68 -x $HG19 -U $f -S $i.sam
        samtools view -1 -q 10 --threads 68 $i.sam | samtools sort --threads 68 > $i.sorted.bam
        makeTagDirectory flip.$i.mTD -genome hg19 -checkGC -tbp 3 -flip $i.sorted.bam &
done
wait
makeMultiWigHub.pl $runname hg19 -url $visdir -webdir $visdir -d flip.*.mTD -force -strand &
wait
