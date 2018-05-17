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
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
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


######### call peaks seperately, with MACS2 #########
for f in $fastq_dir/*.fastq
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq//'`
	simpre=`echo $i | sed 's/NS9-Xiaoyu-MCF7-//;s/_ATAC[0-9]_S[0-9]//;'`
	macs2 callpeak --verbose 3 -t marked_duplicates.$i.sorted.bam -f BAM -g hs -n cutsites.$simpre --keep-dup all --nomodel --shift -100 --extsize 200 --bdg --SPMR 1> cutsites.$simpre.log 2>&1 &
	# macs2 callpeak --verbose 3 -t marked_duplicates.$i.sorted.bam -f BAM -g hs -n mononucleosome.$simpre --keep-dup all --nomodel --shift -37 --extsize 73 --bdg --SPMR 1> mononucleosome.$simpre.log 2>&1 &
done
#####################################################
wait

for f in cutsites.*.bdg
do
	bedSort $f $f.srt
	bedGraphToBigWig $f.srt $genomesize $f.bw &
done
wait

######### find the differential chromatin accessability #########
# https://doi.org/10.1016/j.cell.2016.05.052
# Nfib Promotes Metastasis through a Widespread Increase in Chromatin Accessibility
# bedtools multicov
# DESeq2
# log2 fold change was 0.5 at an FDR <0.1
for pairs in shNUP53-2_Dox_E2:shNUP53-2_E2 shNUP93-1_Dox_E2:shNUP93-1_E2
do
	treats=`echo $pairs | cut -f 1 -d":"`
	controls=`echo $pairs | cut -f 2 -d":"`
	mergePeaks cutsites.${treats}_peaks.narrowPeak cutsites.${controls}_peaks.narrowPeak 1> $treats.$controls.mergePeaks.txt 2> $treats.$controls.mergePeaks.log &
done
wait
for peaks in *.mergePeaks.txt
do
	expname=`echo $peaks | sed 's/.mergePeaks.txt//'`
	annotatePeaks.pl $peaks hg19 -size given -d *.mTD/ 1> $expname.depth.txt 2> $expname.depth.log &
done
expname="shNUP53-2_Dox_E2.shNUP53-2_E2"
head -1 $expname.depth.txt | awk '{print "Flag\t"$0}' > $expname.diff.txt
awk -F"\t" '$21/($20+1) > 2' <(tail -n +2 $expname.depth.txt) | sort -k 21n,21nr | awk '{print "Down\t"$0}' >> $expname.diff.txt
awk -F"\t" '$20/($21+1) > 2' <(tail -n +2 $expname.depth.txt) | sort -k 20n,20n | awk '{print "Up\t"$0}' >> $expname.diff.txt
expname="shNUP93-1_Dox_E2.shNUP93-1_E2"
head -1 $expname.depth.txt | awk '{print "Flag\t"$0}' > $expname.diff.txt
awk -F"\t" '$21/($20+1) > 2' <(tail -n +2 $expname.depth.txt) | sort -k 21n,21nr | awk '{print "Down\t"$0}' >> $expname.diff.txt
awk -F"\t" '$20/($21+1) > 2' <(tail -n +2 $expname.depth.txt) | sort -k 20n,20n | awk '{print "Up\t"$0}' >> $expname.diff.txt
#################################################################

######### plot ATAC-Seq profile for NUP peaks #########
annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{print $1,$2,$3,$4,$5,"+"}' Nup53blrp_E2_peaks.narrowPeak) hg19 -size 6000 -hist 25 -d *.mTD 1> Nup53blrp_E2_peaks.ATAC.profile 2>Nup53blrp_E2_peaks.ATAC.log
annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{print $1,$2,$3,$4,$5,"+"}' Nup53blrp_E2_peaks.narrowPeak) hg19 -size 6000 -hist 25 -ghist -d *.mTD 1> Nup53blrp_E2_peaks.ATAC.profile.mat 2>Nup53blrp_E2_peaks.ATAC.mat.log
annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{print $1,$2,$3,$4,$5,"+"}' Nup53blrp_E2_peaks.narrowPeak | sort -k7,7nr) hg19 -size 6000 -hist 25 -ghist -d *.mTD 1> Nup53blrp_E2_peaks.srtbylogFC.ATAC.profile.mat 2>Nup53blrp_E2_peaks.srtbylogFC.ATAC.mat.log
cut -f 1,2,5,8,11 Nup53blrp_E2_peaks.ATAC.profile | sed '1s/NS9-Xiaoyu-MCF7-//g;1s/_ATAC[0-9]_S[0-9]//g;1s/.mTD Coverage//g' > Nup53blrp_E2_peaks.ATAC.profile.coverage
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/plotHomerProfile.R Nup53blrp_E2_peaks.ATAC.profile.coverage Nup53blrp_E2_peaks.ATAC.profile.mat
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/plotHomerProfile.R Nup53blrp_E2_peaks.ATAC.profile.coverage Nup53blrp_E2_peaks.srtbylogFC.ATAC.profile.mat

annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{print $1,$2,$3,$4,$5,"+"}' Nup53blrp_E2_peaks.narrowPeak) hg19 -size 2000 -hist 5 -norm 1e6 -d *.mTD 1> 2k.Nup53blrp_E2_peaks.ATAC.profile 2> 2k.Nup53blrp_E2_peaks.ATAC.log
annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{print $1,$2,$3,$4,$5,"+"}' Nup53blrp_E2_peaks.narrowPeak | sort -k7,7nr) hg19 -size 2000 -hist 5 -ghist -norm 1e6 -d *.mTD 1> 2k.Nup53blrp_E2_peaks.srtbylogFC.ATAC.profile.mat 2> 2k.Nup53blrp_E2_peaks.srtbylogFC.ATAC.mat.log
cut -f 1,2,5,8,11 2k.Nup53blrp_E2_peaks.ATAC.profile | sed '1s/NS9-Xiaoyu-MCF7-//g;1s/_ATAC[0-9]_S[0-9]//g;1s/.mTD Coverage//g' > 2k.Nup53blrp_E2_peaks.ATAC.profile.coverage
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/plotHomerProfile.R 2k.Nup53blrp_E2_peaks.ATAC.profile.coverage 2k.Nup53blrp_E2_peaks.srtbylogFC.ATAC.profile.mat

chippeaks=Nup53blrp_E2_peaks.narrowPeak
pre=Nup53blrp_E2
for peaks in cutsites.sh*_peaks.narrowPeak
do
	tail=`echo $peaks | sed 's/cutsites.//;s/_peaks.narrowPeak//;'`
	mergePeaks $chippeaks $peaks 1> $pre.$tail.mergePeaks.txt 2> $pre.$tail.mergePeaks.log
	tmp=`tail -3 $pre.$tail.mergePeaks.log | cut -f 3 | tr "\n" " "`
	Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/vennplot_from_a_stdin.R $tail $pre $tmp
done
mergePeaks $chippeaks cutsites.shNUP53-2_*_peaks.narrowPeak 1> $pre.shNUP53.mergePeaks.txt 2> $pre.shNUP53.mergePeaks.log
# tmp=`tail -7 $pre.shNUP53.mergePeaks.log | cut -f 4 | tr "\n" " "`
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/vennplot_from_a_stdin.R shNUP53-2_Dox_E2 shNUP53-2_E2 $pre 3161 13743 2165 36778 884 192 10097
#######################################################


######### QC and visualization #########
echo "Experiment Length" | tr " " "\t" > NUP_ATAC.peak.length.tsv
# for peaks in *peaks.narrowPeak
for peaks in cutsites.*_peaks.broadPeak
do
	# expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	expname=`echo $peaks | sed 's/_peaks.broadPeak//;s/cutsites.//'`
	cut -f 1-3 $peaks | uniq | awk -vvar=$expname -vOFS="\t" '{print var,$3-$2+1}'  >> NUP_ATAC.peak.length.tsv
done
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/MultipleDistribution.R NUP_ATAC.peak.length.tsv 1 2

echo "Experiment Length" | tr " " "\t" > NUP_ATAC.nup53.peak.length.tsv
echo "Experiment Length" | tr " " "\t" > NUP_ATAC.nup93.peak.length.tsv
# for peaks in *shNUP53*peaks.narrowPeak
for peaks in cutsites.shNUP53*_peaks.broadPeak
do
	# expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	expname=`echo $peaks | sed 's/_peaks.broadPeak//;s/cutsites.//'`
	cut -f 1-3 $peaks | uniq | awk -vvar=$expname -vOFS="\t" '{print var,$3-$2+1}'  >> NUP_ATAC.nup53.peak.length.tsv
done
for peaks in cutsites.shNUP93*_peaks.broadPeak
do
	# expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	expname=`echo $peaks | sed 's/_peaks.broadPeak//;s/cutsites.//'`
	cut -f 1-3 $peaks | uniq | awk -vvar=$expname -vOFS="\t" '{print var,$3-$2+1}'  >> NUP_ATAC.nup93.peak.length.tsv
done
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/MultipleDistribution.R NUP_ATAC.nup53.peak.length.tsv 1 2
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/MultipleDistribution.R NUP_ATAC.nup93.peak.length.tsv 1 2


# for peaks in *peaks.narrowPeak
for peaks in cutsites.shNUP93*_peaks.broadPeak
do
	# expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	expname=`echo $peaks | sed 's/_peaks.broadPeak//;s/cutsites.//'`
	findMotifsGenome.pl <(awk -vOFS="\t" '{print $1,$2,$3,$4,"1000","+"}' $peaks | uniq) hg19 $expname.homer.motifs -size given -cpu 68 -mknown /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs > $expname.findMotifsGenome.txt &
	annotatePeaks.pl <(awk -vOFS="\t" '{print $1,$2,$3,$4,"1000","+"}' $peaks | uniq) hg19 -m /home1/04935/shaojf/myTools/HOMER/data/knownTFs/vertebrates/known.motifs -size given -cpu 68 1> $expname.motifs.txt 2> $expname.motifs.log &
done

# for peaks in *peaks.narrowPeak
for peaks in cutsites.shNUP93*_peaks.broadPeak
do
	# expname=`echo $peaks | sed 's/_peaks.narrowPeak//'`
	expname=`echo $peaks | sed 's/_peaks.broadPeak//;s/cutsites.//'`
	annotatePeaks.pl <(awk -vOFS="\t" '{print $1,$2,$3,$4,"1000","+"}' $peaks | uniq) hg19 -size given -d *.mTD/ 1> $expname.simanno.txt 2> $expname.simanno.log &
done

for pairs in NS9-Xiaoyu-MCF7-shNUP53-2_Dox_E2_ATAC2_S2:NS9-Xiaoyu-MCF7-shNUP53-2_E2_ATAC1_S1 NS9-Xiaoyu-MCF7-shNUP93-1_Dox_E2_ATAC4_S4:NS9-Xiaoyu-MCF7-shNUP93-1_E2_ATAC3_S3
do
	treats=`echo $pairs | cut -f 1 -d":"`
	controls=`echo $pairs | cut -f 2 -d":"`
	mergePeaks ${treats}_peaks.narrowPeak ${controls}_peaks.narrowPeak 1> $treats.$controls.mergePeaks.txt 2> $treats.$controls.mergePeaks.log
done
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/vennplot_from_a_stdin.R shNUP53-2_Dox_E2 shNUP53-2_E2 9534 17120 61058
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/vennplot_from_a_stdin.R shNUP93-2_Dox_E2 shNUP93-2_E2 3562 30895 51030
