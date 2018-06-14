#!/bin/sh
# rclone sync /home1/04935/shaojf/scratch/TIP60.project/ mygoogle:TIP60.KAT5/TIP60.project/
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

# cd /home1/04935/shaojf/scratch/TIP60.project/

#### additional analysis for ATAC-Seq
workdir=/home1/04935/shaojf/scratch/TIP60.project/ATAC-seq
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
	-o $trimadapterdir/$i.clean.fastq $f &
done
wait
###### my pipeline
for f in $fastq_dir/*.fastq.gz
do
	# echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	echo $i

	bowtie2 --local --very-sensitive -p 68 -x $HG19 -U $trimadapterdir/$i.clean.fastq -S $i.sam
	samtools view -1 -q 10 --threads 68 $i.sam | samtools sort --threads 68 > $i.sorted.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 $i.sorted.bam
	makeBigWig.pl $i.mTD/ hg19 -url $visdir -webdir $visdir &
	findPeaks $i.mTD -style factor -o auto &
	picard MarkDuplicates I=$i.sorted.bam O=marked_duplicates.$i.sorted.bam M=$i.marked_dup_metrics.txt
	# macs2 callpeak --verbose 3 -t marked_duplicates.$i.sorted.bam -f BAM -g hs -n $i --call-summits --keep-dup all &
	macs2 callpeak --verbose 3 -t marked_duplicates.$i.sorted.bam -f BAM -g hs -n $i --call-summits --keep-dup 1 &
done
makeMultiWigHub.pl Feng_ATAC hg19 -url $visdir -webdir $visdir -d *.mTD
wait

for f in $fastq_dir/*.fastq.gz
do
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	echo $i
	samtools markdup -r -s -O BAM -@ 68 $i.sorted.bam $i.rmdup.bam &
	# samtools index -@ 68 $i.rmdup.bam &
done

for f in $fastq_dir/*.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq.gz//'`
	simpre=`echo $i | sed 's/NS9-Feng-//;s/-ATAC[0-9]//;'`
	macs2 callpeak --verbose 3 -t marked_duplicates.$i.sorted.bam -f BAM -g hs -n cutsites.$simpre --keep-dup 1 --nomodel --shift -100 --extsize 200 --bdg --SPMR 1> cutsites.$simpre.log 2>&1 &
done

mergePeaks cutsites.Hela-*_peaks.narrowPeak 1> Hela.ATAC.mergePeaks.txt 2> Hela.ATAC.mergePeaks.log
Rscript /home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/vennplot_from_a_stdin.R Hela-siCTL Hela-siP400 Hela-siTIP60 5007 5997 16092 1796 5955 6431 31599
# login2(1245)$ tail -8 Hela.ATAC.mergePeaks.log 
# cutsites.Hela-siCTL_peaks.narrowPeak	cutsites.Hela-siP400_peaks.narrowPeak	cutsites.Hela-siTIP60_peaks.narrowPeak	Total	Name
# 		X	16092	cutsites.Hela-siTIP60_peaks.narrowPeak
# 	X		5997	cutsites.Hela-siP400_peaks.narrowPeak
# 	X	X	5955	cutsites.Hela-siP400_peaks.narrowPeak|cutsites.Hela-siTIP60_peaks.narrowPeak
# X			5007	cutsites.Hela-siCTL_peaks.narrowPeak
# X		X	6431	cutsites.Hela-siCTL_peaks.narrowPeak|cutsites.Hela-siTIP60_peaks.narrowPeak
# X	X		1796	cutsites.Hela-siCTL_peaks.narrowPeak|cutsites.Hela-siP400_peaks.narrowPeak
# X	X	X	31599	cutsites.Hela-siCTL_peaks.narrowPeak|cutsites.Hela-siP400_peaks.narrowPeak|cutsites.Hela-siTIP60_peaks.narrowPeak
for peaks in cutsites.Hela-*_peaks.narrowPeak
do
	expname=`echo $peaks | sed 's/cutsites.//;s/_peaks.narrowPeak//'`
	annotatePeaks.pl $peaks hg19 -size given 1> $expname.depth.txt 2> $expname.depth.log &
done
expname=Hela.ATAC
annotatePeaks.pl Hela.ATAC.mergePeaks.txt hg19 -size given -d *.mTD/ 1> $expname.depth.txt 2> $expname.depth.log &


####################################################################################
#### dREG
# https://github.com/Danko-Lab/dREG
# https://github.com/Danko-Lab/dREG-Model
# https://github.com/Danko-Lab/dREG.HD
# export PATH=$HOME/myTools/dREG:$HOME/myTools/dREG/dREG.HD:$PATH

### bam to bigwig
bigwigdirs=/home1/04935/shaojf/scratch/TIP60.project/GRO-Seq/bigwigs
url=http://sjf.dingding.biz/trackhubs
for f in Hela-*-GROseq.sorted.bam
do
	pre=`echo $f | sed 's/.sorted.bam//'`
	makeTagDirectory $pre.mTD -genome hg19 -checkGC -tbp 3 -fragLength 200 $f
	makeBigWig.pl $pre.mTD/ hg19 -url $url -webdir $bigwigdirs -force -strand
done

# Read counts (not normalized) formatted as a bigWig file.
for f in Hela-*-GROseq.sorted.bam
do
	pre=`echo $f | sed 's/.sorted.bam//'`
	samtools markdup -r -s -O BAM -@ 68 $f $pre.rmdup.bam
	samtools index -@ 68 $pre.rmdup.bam
	bamCoverage -b $pre.rmdup.bam -o bamCoverage.notNormalized/$pre.plus.bigwig -of bigwig --filterRNAstrand forward
	bamCoverage -b $pre.rmdup.bam -o bamCoverage.notNormalized/$pre.minus.bigwig -of bigwig --filterRNAstrand reverse
done

genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
for f in Hela-*-GROseq.sorted.bam
do
	pre=`echo $f | sed 's/.sorted.bam//'`
	bedSort bamCoverage.notNormalized/$pre.plus.bigwig.bg bamCoverage.notNormalized/$pre.plus.bigwig.bg.srt
	bedGraphToBigWig bamCoverage.notNormalized/$pre.plus.bigwig.bg.srt $genomesize bamCoverage.notNormalized/bedGraphToBigWig.$pre.plus.bigwig &
	bedSort bamCoverage.notNormalized/$pre.minus.bigwig.bg bamCoverage.notNormalized/$pre.minus.bigwig.bg.srt
	bedGraphToBigWig bamCoverage.notNormalized/$pre.plus.bigwig.bg.srt $genomesize bamCoverage.notNormalized/bedGraphToBigWig.$pre.minus.bigwig &
done

dregmodel=/home1/04935/shaojf/myTools/dREG/dREG-Model/asvm.gdm.6.6M.20170828.rdata
while read line
do
	gro=`echo $line | awk '{print $1}'`
	negbw=bedGraphToBigWig.${gro}.plus.bigwig
	posbw=bedGraphToBigWig.${gro}.minus.bigwig
	bash $HOME/myTools/dREG/run_peakcalling.bsh $posbw $negbw $gro.dREG.peak $dregmodel 68
	bash $HOME/myTools/dREG/run_dREG.bsh $posbw $negbw $gro.dREG $gro.dREG.score $dregmodel 68
	bash $HOME/myTools/dREG/writeBed.bsh 0.25 $gro.dREG.score.bedGraph.gz
done < pairs.txt

hdmodel=/home1/04935/shaojf/myTools/dREG/dREG.HD/dREG.HD/inst/extdata/dREG_HD.model.rdata
while read line
do
	gro=`echo $line | awk '{print $1}'`
	atac=`echo $line | awk '{print $2}'`
	negbw=${gro}.mTDneg.ucsc.bigWig
	posbw=${gro}.mTDpos.ucsc.bigWig
	bash $HOME/myTools/dREG/dREG.HD/run_dREG-HD.bsh ${gro}.dREG.bed $posbw $negbw $hdmodel 30
done < pairs.txt
####################################################################################
# Hela-siCTL-rpt1-GROseq NS9-Feng-Hela-siCTL-ATAC1
# Hela-siCTL-rpt2-GROseq NS9-Feng-Hela-siCTL-ATAC1
# Hela-siCTL-rpt3-GROseq NS9-Feng-Hela-siCTL-ATAC1
# Hela-siP400-rpt1-GROseq NS9-Feng-Hela-siP400-ATAC2
# Hela-siP400-rpt2-GROseq NS9-Feng-Hela-siP400-ATAC2
# Hela-siP400-rpt3-GROseq NS9-Feng-Hela-siP400-ATAC2
# Hela-siTIP60-rpt1-GROseq NS9-Feng-Hela-siTIP60-ATAC3
# Hela-siTIP60-rpt2-GROseq NS9-Feng-Hela-siTIP60-ATAC3
# Hela-siTIP60-rpt3-GROseq NS9-Feng-Hela-siTIP60-ATAC3

####################################################################################
# NRSA
####step 1
pause_PROseq.pl -o siP400vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siP400-rpt1-GROseq.sorted.bam  Hela-siP400-rpt2-GROseq.sorted.bam  Hela-siP400-rpt3-GROseq.sorted.bam -m hg19 &
pause_PROseq.pl -o siTIP60vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siTIP60-rpt1-GROseq.sorted.bam  Hela-siTIP60-rpt2-GROseq.sorted.bam  Hela-siTIP60-rpt3-GROseq.sorted.bam -m hg19 &
wait
####step 2
# -pri 1 -dir 0 -peak ../Hdac3_ZT10_peaks.narrowPeak -lk pp -cf 0.001
eRNA.pl -w siP400vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siP400-rpt1-GROseq.sorted.bam  Hela-siP400-rpt2-GROseq.sorted.bam  Hela-siP400-rpt3-GROseq.sorted.bam -m hg19 &
eRNA.pl -w siP400vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siTIP60-rpt1-GROseq.sorted.bam  Hela-siTIP60-rpt2-GROseq.sorted.bam  Hela-siTIP60-rpt3-GROseq.sorted.bam -m hg19 &
wait
####################################################################################



####################################################################################
######### QC, use deeptools #########
#### add ChIP-Seq, Gro-Seq, as a control
# ls | sed 's/.ucsc.bigWig//;s/NS9-Feng-//;s/[1-3].mTD//' | tr "\n" " "
multiBigwigSummary bins --binSize 1000 --bwfiles *ucsc.bigWig --outFileName Hela.bin1000.out.npz --outRawCounts Hela.bin1000.out.coverage
plotCorrelation -in Hela.bin1000.out.npz -o Hela.bin1000.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --labels Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC 

cut -f 2-4 Hela.ATAC.mergePeaks.txt | tail -n +2 | uniq | awk '{print $0"\t"NR}' > a.bed
multiBigwigSummary BED-file --BED a.bed --bwfiles *ucsc.bigWig --outFileName Hela.ATAC_peaks.out.npz --outRawCounts Hela.ATAC_peaks.out.coverage
plotCorrelation -in Hela.ATAC_peaks.out.npz -o Hela.ATAC_peaks.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --labels Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC 

plotFingerprint --bamfiles *.bam -plot Hela.ATAC.fingerprint.pdf --ignoreDuplicates --plotFileFormat pdf

computeMatrix reference-point -R a.bed -S *ucsc.bigWig -out Hela.ATAC_peaks.computeMatrix.gz --referencePoint center -b 3000 -a 3000 -bs 10 
plotHeatmap -m Hela.ATAC_peaks.computeMatrix.gz -o Hela.ATAC_peaks.computeMatrix.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC --heatmapWidth 20

plotProfile -m Hela.ATAC_peaks.computeMatrix.gz -o Hela.ATAC_peaks.computeMatrix.plotProfile.pdf --kmeans 1 --plotFileFormat pdf --samplesLabel Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC --plotWidth 20 --plotHeight 10 --refPointLabel PeakCenter --plotType fill
####################################################################################

