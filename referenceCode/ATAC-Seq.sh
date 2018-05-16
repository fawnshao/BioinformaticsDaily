#!/bin/bash
# this code will generate .sorted.bam files first, and then generate Homer tag directory and UCSC file as track hubs
####
# rclone sync /scratch/04935/shaojf/NUP_related_others/ATAC-seq/ mygoogle:NUP_project/NUP_related_others/ATAC-seq/
####
workdir=/home1/04935/shaojf/scratch/NUP_related_others/ATAC-seq
fastq_dir=$workdir/fastqs
fastqc_dir=$workdir/fastqc_res
trimadapterdir=$workdir/fastqs_trimmed
# mappingtogenome=$workdir/reads2genome
# bw_dir=$workdir/bigwigs
# dhspeaks=$workdir/dhspeaks

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
visdir=/home1/04935/shaojf/stampede2/UCSCvis

mkdir $fastqc_dir
mkdir $trimadapterdir

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
	--info-file=$i.cutadapt.metrics \
	-o $trimadapterdir/$i.clean.fastq $f &
	echo    ***____step0: running cutadapt successful____***
done
wait
# grep -e "=== Adapter" -e "Trimmed" ATAC.out > cutadapt.stats

###### esATAC
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("esATAC")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")


###### my pipeline
for f in $fastq_dir/*.fastq
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq//'`

	bowtie2 --local --very-sensitive -p 68 -x $HG19 -U $trimadapterdir/$i.clean.fastq | \
	samtools view -1 -q 10 --threads 68 | samtools sort --threads 68 > $i.sorted.bam
	echo    ***____step1: bowtie2 successful____***

	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
	echo    ***____step2: home makeTagDirectory successful____***
	echo    ***--- these samples were normed to 10million total reads____***

	makeBigWig.pl $i.mTD/ hg19 -url $visdir -webdir $visdir &
	echo    ***____step3: minus trand UCSC bedgraph.gz file of $i.mTD folder completed____***

	findPeaks $i.mTD -style factor -o auto &
	echo    ***____step4: findPeaks for $i.mTD folder completed____***

	# samtools rmdup $i.sorted.bam $i.sorted.rmdup.bam
	picard MarkDuplicates I=$i.sorted.bam O=marked_duplicates.$i.sorted.bam M=$i.marked_dup_metrics.txt
	macs2 callpeak --verbose 3 -t marked_duplicates.$i.sorted.bam -f BAM -g hs -n $i --call-summits --keep-dup all &
	
	echo    ++++++++++++++++++++++++ Yahoo - samples all finished  +++++++++++++++++++++++++++++++++++
done

makeMultiWigHub.pl NUP_ATAC hg19 -url $visdir -webdir $visdir -d *.mTD
wait

echo "Experiment Length" | tr " " "\t" > NUP_ATAC.peak.length.tsv
for peaks in *peaks.narrowPeak
do
	expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	cut -f 1-3 $peaks | uniq | awk -vvar=$expname -vOFS="\t" '{print var,$3-$2+1}'  >> NUP_ATAC.peak.length.tsv
done
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/MultipleDistribution.R NUP_ATAC.peak.length.tsv 1 2

echo "Experiment Length" | tr " " "\t" > NUP_ATAC.nup53.peak.length.tsv
echo "Experiment Length" | tr " " "\t" > NUP_ATAC.nup93.peak.length.tsv
for peaks in *shNUP53*peaks.narrowPeak
do
	expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	cut -f 1-3 $peaks | uniq | awk -vvar=$expname -vOFS="\t" '{print var,$3-$2+1}'  >> NUP_ATAC.nup53.peak.length.tsv
done
for peaks in *shNUP93*peaks.narrowPeak
do
	expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	cut -f 1-3 $peaks | uniq | awk -vvar=$expname -vOFS="\t" '{print var,$3-$2+1}'  >> NUP_ATAC.nup93.peak.length.tsv
done
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/MultipleDistribution.R NUP_ATAC.nup53.peak.length.tsv 1 2
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/MultipleDistribution.R NUP_ATAC.nup93.peak.length.tsv 1 2


for peaks in *peaks.narrowPeak
do
	expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{print $1,$2,$3,$4,"1000","+"}' $peaks | uniq) hg19 $expname.homer.motifs -size given -cpu 68 -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs > $expname.findMotifsGenome.txt &
	annotatePeaks.pl <(awk -vOFS="\t" '{print $1,$2,$3,$4,"1000","+"}' $peaks | uniq) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $expname.motifs.txt 2> $expname.motifs.log &
done

for peaks in *peaks.narrowPeak
do
	expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	annotatePeaks.pl <(awk -vOFS="\t" '{print $1,$2,$3,$4,"1000","+"}' $peaks | uniq) hg19 -size given -d *.mTD/ 1> $expname.simanno.txt 2> $expname.simanno.log &
done

for pairs in NS9-Xiaoyu-MCF7-shNUP53-2_Dox_E2_ATAC2_S2:NS9-Xiaoyu-MCF7-shNUP53-2_E2_ATAC1_S1 NS9-Xiaoyu-MCF7-shNUP93-1_Dox_E2_ATAC4_S4:NS9-Xiaoyu-MCF7-shNUP93-1_E2_ATAC3_S3
do
	treats=`echo $pairs | cut -f 1 -d":"`
	controls=`echo $pairs | cut -f 2 -d":"`
	mergePeaks ${treats}_peaks.narrowPeak ${controls}_peaks.narrowPeak 1> $treats.$controls.mergePeaks.txt 2> $treats.$controls.mergePeaks.log
done

