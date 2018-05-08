#!/bin/bash
# this code will generate .sorted.bam files first, and then generate Homer tag directory and UCSC file as track hubs

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
visdir=/home1/04935/shaojf/stampede2/UCSCvis

for f in MCF7-*.fastq.gz
do
	echo $f

	i=`echo $f | sed 's/.fastq.gz//'`

	nice bowtie2 -p 48 -x $HG19 -U $f | samtools view -1 -q 10 --threads 48 | samtools sort --threads 48 > $i.sorted.bam
	echo    ***____step1: bowtie2 successful ***____

	nice makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
	echo    ***____step2: home makeTagDirectory successful ***____
	echo    ***--- these samples were normed to 10million total reads ****------

	# nice makeUCSCfile $i.mTD -o > $i.mTD/$i.mTD.bedgraph.gz
	nice makeBigWig.pl $i.mTD/ hg19 -url $visdir -webdir $visdir
	echo    ***____step3: minus trand UCSC bedgraph.gz file of $i.mTD folder completed ***____

	nice findPeaks $i.mTD -style factor -o auto
	echo    ***____step4: findPeaks for $i.mTD folder completed ***____
	echo    ++++++++++++++++++++++++ Yahoo - samples all finished  +++++++++++++++++++++++++++++++++++

	samtools rmdup $i.sorted.bam $i.sorted.rmdup.bam
	macs2 callpeak -t $i.sorted.rmdup.bam -f BAM -g hs -n $i --fix-bimodal --call-summits
done

for f in MDAMB231-*.fastq.gz
do
	echo $f

	i=`echo $f | sed 's/.fastq.gz//'`

	nice bowtie2 -p 48 -x $HG19 -U $f | samtools view -1 -q 10 --threads 48 | samtools sort --threads 48 > $i.sorted.bam
	echo    ***____step1: bowtie2 successful ***____

	nice makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
	echo    ***____step2: home makeTagDirectory successful ***____
	echo    ***--- these samples were normed to 10million total reads ****------

	# nice makeUCSCfile $i.mTD -o > $i.mTD/$i.mTD.bedgraph.gz
	nice makeBigWig.pl $i.mTD/ hg19 -url $visdir -webdir $visdir
	echo    ***____step3: minus trand UCSC bedgraph.gz file of $i.mTD folder completed ***____

	samtools rmdup $i.sorted.bam $i.sorted.rmdup.bam
done
nice findPeaks MDAMB231-NUP93-bethyl-DC_S13.mTD -style factor -o auto -i MDAMB231-Input-DC_S15.mTD
macs2 callpeak -t MDAMB231-NUP93-bethyl-DC_S13.sorted.rmdup.bam \
-c MDAMB231-Input-DC_S15.sorted.rmdup.bam -f BAM -g hs \
-n MDAMB231-NUP93-bethyl-DC_S13 --fix-bimodal --call-summits

makeMultiWigHub.pl NUP_ChIP hg19 -url $visdir -webdir $visdir -d *.mTD
