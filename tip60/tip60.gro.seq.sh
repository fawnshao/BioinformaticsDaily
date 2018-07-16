#!/bin/sh
# rclone sync -L /home1/04935/shaojf/scratch/TIP60.project/ mygoogle:TIP60.KAT5/TIP60.project/
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
mydeseq=/home1/04935/shaojf/myTools/BioinformaticsDaily/referenceCode/runDESeq2.R
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


# ####################################################################################
# #### dREG
# # https://github.com/Danko-Lab/dREG
# # https://github.com/Danko-Lab/dREG-Model
# # https://github.com/Danko-Lab/dREG.HD
# # export PATH=$HOME/myTools/dREG:$HOME/myTools/dREG/dREG.HD:$PATH

# ### bam to bigwig
# bigwigdirs=/home1/04935/shaojf/scratch/TIP60.project/GRO-Seq/bigwigs
# url=http://sjf.dingding.biz/trackhubs
# for f in Hela-*-GROseq.sorted.bam
# do
# 	pre=`echo $f | sed 's/.sorted.bam//'`
# 	makeTagDirectory $pre.mTD -genome hg19 -checkGC -tbp 3 -fragLength 200 $f
# 	makeBigWig.pl $pre.mTD/ hg19 -url $url -webdir $bigwigdirs -force -strand
# done

# # Read counts (not normalized) formatted as a bigWig file.
# for f in Hela-*-GROseq.sorted.bam
# do
# 	pre=`echo $f | sed 's/.sorted.bam//'`
# 	samtools markdup -r -s -O BAM -@ 68 $f $pre.rmdup.bam
# 	samtools index -@ 68 $pre.rmdup.bam
# 	bamCoverage -b $pre.rmdup.bam -o bamCoverage.notNormalized/$pre.plus.bigwig -of bigwig --filterRNAstrand forward
# 	bamCoverage -b $pre.rmdup.bam -o bamCoverage.notNormalized/$pre.minus.bigwig -of bigwig --filterRNAstrand reverse
# done

# genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
# for f in Hela-*-GROseq.sorted.bam
# do
# 	pre=`echo $f | sed 's/.sorted.bam//'`
# 	bedSort bamCoverage.notNormalized/$pre.plus.bigwig.bg bamCoverage.notNormalized/$pre.plus.bigwig.bg.srt
# 	bedGraphToBigWig bamCoverage.notNormalized/$pre.plus.bigwig.bg.srt $genomesize bamCoverage.notNormalized/bedGraphToBigWig.$pre.plus.bigwig &
# 	bedSort bamCoverage.notNormalized/$pre.minus.bigwig.bg bamCoverage.notNormalized/$pre.minus.bigwig.bg.srt
# 	bedGraphToBigWig bamCoverage.notNormalized/$pre.plus.bigwig.bg.srt $genomesize bamCoverage.notNormalized/bedGraphToBigWig.$pre.minus.bigwig &
# done

# dregmodel=/home1/04935/shaojf/myTools/dREG/dREG-Model/asvm.gdm.6.6M.20170828.rdata
# while read line
# do
# 	gro=`echo $line | awk '{print $1}'`
# 	negbw=bedGraphToBigWig.${gro}.plus.bigwig
# 	posbw=bedGraphToBigWig.${gro}.minus.bigwig
# 	bash $HOME/myTools/dREG/run_peakcalling.bsh $posbw $negbw $gro.dREG.peak $dregmodel 68
# 	bash $HOME/myTools/dREG/run_dREG.bsh $posbw $negbw $gro.dREG $gro.dREG.score $dregmodel 68
# 	bash $HOME/myTools/dREG/writeBed.bsh 0.25 $gro.dREG.score.bedGraph.gz
# done < pairs.txt

# hdmodel=/home1/04935/shaojf/myTools/dREG/dREG.HD/dREG.HD/inst/extdata/dREG_HD.model.rdata
# while read line
# do
# 	gro=`echo $line | awk '{print $1}'`
# 	atac=`echo $line | awk '{print $2}'`
# 	negbw=${gro}.mTDneg.ucsc.bigWig
# 	posbw=${gro}.mTDpos.ucsc.bigWig
# 	bash $HOME/myTools/dREG/dREG.HD/run_dREG-HD.bsh ${gro}.dREG.bed $posbw $negbw $hdmodel 30
# done < pairs.txt
# ####################################################################################
# # Hela-siCTL-rpt1-GROseq NS9-Feng-Hela-siCTL-ATAC1
# # Hela-siCTL-rpt2-GROseq NS9-Feng-Hela-siCTL-ATAC1
# # Hela-siCTL-rpt3-GROseq NS9-Feng-Hela-siCTL-ATAC1
# # Hela-siP400-rpt1-GROseq NS9-Feng-Hela-siP400-ATAC2
# # Hela-siP400-rpt2-GROseq NS9-Feng-Hela-siP400-ATAC2
# # Hela-siP400-rpt3-GROseq NS9-Feng-Hela-siP400-ATAC2
# # Hela-siTIP60-rpt1-GROseq NS9-Feng-Hela-siTIP60-ATAC3
# # Hela-siTIP60-rpt2-GROseq NS9-Feng-Hela-siTIP60-ATAC3
# # Hela-siTIP60-rpt3-GROseq NS9-Feng-Hela-siTIP60-ATAC3

# ####################################################################################
# # NRSA
# ####step 1
# pause_PROseq.pl -o siP400vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siP400-rpt1-GROseq.sorted.bam  Hela-siP400-rpt2-GROseq.sorted.bam  Hela-siP400-rpt3-GROseq.sorted.bam -m hg19 &
# pause_PROseq.pl -o siTIP60vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siTIP60-rpt1-GROseq.sorted.bam  Hela-siTIP60-rpt2-GROseq.sorted.bam  Hela-siTIP60-rpt3-GROseq.sorted.bam -m hg19 &
# wait
# ####step 2
# # -pri 1 -dir 0 -peak ../Hdac3_ZT10_peaks.narrowPeak -lk pp -cf 0.001
# eRNA.pl -w siP400vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siP400-rpt1-GROseq.sorted.bam  Hela-siP400-rpt2-GROseq.sorted.bam  Hela-siP400-rpt3-GROseq.sorted.bam -m hg19 &
# eRNA.pl -w siP400vsCTL -in1 Hela-siCTL-rpt1-GROseq.sorted.bam  Hela-siCTL-rpt2-GROseq.sorted.bam  Hela-siCTL-rpt3-GROseq.sorted.bam -in2 Hela-siTIP60-rpt1-GROseq.sorted.bam  Hela-siTIP60-rpt2-GROseq.sorted.bam  Hela-siTIP60-rpt3-GROseq.sorted.bam -m hg19 &
# wait
# ####################################################################################



# ####################################################################################
# ######### QC, use deeptools #########
# #### add ChIP-Seq, Gro-Seq, as a control
# # ls | sed 's/.ucsc.bigWig//;s/NS9-Feng-//;s/[1-3].mTD//' | tr "\n" " "
# multiBigwigSummary bins --binSize 1000 --bwfiles *ucsc.bigWig --outFileName Hela.bin1000.out.npz --outRawCounts Hela.bin1000.out.coverage
# plotCorrelation -in Hela.bin1000.out.npz -o Hela.bin1000.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --labels Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC 

# cut -f 2-4 Hela.ATAC.mergePeaks.txt | tail -n +2 | uniq | awk '{print $0"\t"NR}' > a.bed
# multiBigwigSummary BED-file --BED a.bed --bwfiles *ucsc.bigWig --outFileName Hela.ATAC_peaks.out.npz --outRawCounts Hela.ATAC_peaks.out.coverage
# plotCorrelation -in Hela.ATAC_peaks.out.npz -o Hela.ATAC_peaks.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --labels Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC 

# # plotFingerprint --bamfiles *.bam -plot Hela.ATAC.fingerprint.pdf --ignoreDuplicates --plotFileFormat pdf

# computeMatrix reference-point -R a.bed -S *ucsc.bigWig -out Hela.ATAC_peaks.computeMatrix.gz --referencePoint center -b 3000 -a 3000 -bs 10 
# plotHeatmap -m Hela.ATAC_peaks.computeMatrix.gz -o Hela.ATAC_peaks.computeMatrix.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC --heatmapWidth 20

# plotProfile -m Hela.ATAC_peaks.computeMatrix.gz -o Hela.ATAC_peaks.computeMatrix.plotProfile.pdf --kmeans 1 --plotFileFormat pdf --samplesLabel Hela-siCTL-rpt1-GROseq.mTDneg Hela-siCTL-rpt1-GROseq.mTDpos Hela-siCTL-rpt2-GROseq.mTDneg Hela-siCTL-rpt2-GROseq.mTDpos Hela-siCTL-rpt3-GROseq.mTDneg Hela-siCTL-rpt3-GROseq.mTDpos Hela-siP400-rpt1-GROseq.mTDneg Hela-siP400-rpt1-GROseq.mTDpos Hela-siP400-rpt2-GROseq.mTDneg Hela-siP400-rpt2-GROseq.mTDpos Hela-siP400-rpt3-GROseq.mTDneg Hela-siP400-rpt3-GROseq.mTDpos Hela-siTIP60-rpt1-GROseq.mTDneg Hela-siTIP60-rpt1-GROseq.mTDpos Hela-siTIP60-rpt2-GROseq.mTDneg Hela-siTIP60-rpt2-GROseq.mTDpos Hela-siTIP60-rpt3-GROseq.mTDneg Hela-siTIP60-rpt3-GROseq.mTDpos Hela-siCTL-ATAC Hela-siP400-ATAC Hela-siTIP60-ATAC --plotWidth 20 --plotHeight 10 --refPointLabel PeakCenter --plotType fill
# ####################################################################################

# ####################################################################################
# # enhancer heatmap
# labs=`ls Hela-si*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//'`
# awk '$NF=="-"' hela.enhancer.bed > hela.enhancer.neg.bed
# awk '$NF=="+"' hela.enhancer.bed > hela.enhancer.pos.bed
# # for reg in hela.enhancer.neg.bed hela.enhancer.pos.bed
# for reg in hela.enhancer.bed
# do
# 	pre=hela.enhancer.gro
# 	computeMatrix reference-point -R $reg -S Hela-si*.ucsc.bigWig -out $pre.computeMatrix.gz --referencePoint center -b 3000 -a 3000 -bs 10 -p 48
# 	plotHeatmap -m $pre.computeMatrix.gz -o $pre.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 --refPointLabel EnhancerCenter &
# 	plotProfile -m $pre.computeMatrix.gz -o $pre.plotProfile.pdf --kmeans 3 --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --refPointLabel EnhancerCenter --numPlotsPerRow 2 &
# done
# wait
# ####use feng bigwig files
# for reg in enhancer_siTIP60.bed genes_siTIP60.bed TSS_siTIP60.bed
# do
# 	pre=hela.gro
# 	labs=`ls Hela-si{CTL,TIP60}*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
# 	computeMatrix reference-point -R Increase_$reg Decrease_$reg -S Hela-siCTL*.ucsc.bigWig Hela-siTIP60*.ucsc.bigWig -out $pre.$reg.computeMatrix.gz --regionBodyLength 5000 -b 3000 -a 3000 -bs 10 -p 48
# 	plotHeatmap -m $pre.$reg.computeMatrix.gz -o $pre.$reg.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 --regionsLabel Increase_$reg Decrease_$reg &
# 	plotProfile -m $pre.$reg.computeMatrix.gz -o $pre.$reg.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --numPlotsPerRow 2 --regionsLabel Increase_$reg Decrease_$reg &
# done
# wait


# for reg in *enhancer_siTIP60.bed *TSS_siTIP60.bed
# do
# 	pre=hela.gro
# 	labs=`ls Hela-si{CTL,TIP60}*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
# 	computeMatrix reference-point -R $reg -S Hela-siCTL*.ucsc.bigWig Hela-siTIP60*.ucsc.bigWig -out $pre.$reg.computeMatrix.gz --referencePoint center -b 5000 -a 5000 -bs 100 -p 48
# 	plotHeatmap -m $pre.$reg.computeMatrix.gz -o $pre.$reg.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 &
# 	plotProfile -m $pre.$reg.computeMatrix.gz -o $pre.$reg.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 &
# done
# for reg in *genes_siTIP60.bed
# do
# 	pre=hela.gro
# 	labs=`ls Hela-si{CTL,TIP60}*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
# 	computeMatrix reference-point -R $reg -S Hela-siCTL*.ucsc.bigWig Hela-siTIP60*.ucsc.bigWig -out $pre.$reg.computeMatrix.gz --regionBodyLength 5000 -b 5000 -a 5000 -bs 100 -p 48
# 	plotHeatmap -m $pre.$reg.computeMatrix.gz -o $pre.$reg.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 &
# 	plotProfile -m $pre.$reg.computeMatrix.gz -o $pre.$reg.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup &
# done
# wait

# bedtools merge -d 5000 -i hela.enhancer.bed > merged.hela.enhancer.bed

#################################################################################
# use sjf eRNA list
cat intergenic.eRNA.bed intragenic.eRNA.bed | bedtools sort -i - > Hela.eRNA.basedonENCODE.bed
awk -vOFS="\t" '{print $4,$1,$2,$3,$6,"2000"}' Hela.eRNA.basedonENCODE.bed > Hela.eRNA.forHomer.txt
analyzeRepeats.pl rna hg19 -d Homer.tags/Hela-si*-GROseq -strand + -raw -condenseGenes > Hela.gene.raw.txt
analyzeRepeats.pl rna hg19 -d Homer.tags/Hela-si*-GROseq -strand + -fpkm -condenseGenes > Hela.gene.fpkm.txt
analyzeRepeats.pl Hela.eRNA.forHomer.txt hg19 -d Homer.tags/Hela-si*-GROseq -strand + -raw > Hela.eRNA.raw.txt
analyzeRepeats.pl Hela.eRNA.forHomer.txt hg19 -d Homer.tags/Hela-si*-GROseq -strand + -fpkm > Hela.eRNA.fpkm.txt

head -1 Hela.gene.raw.txt | cut -f 2-5,9-17 | sed 's?Homer.tags/??g' | tr "\t" "\n" | awk '{print $1}' | tr "\n" "\t" | sed 's?\t$??' | awk '{print "ID\t"$0}' > Hela.gene.eRNA.counts
tail -n +2 Hela.gene.raw.txt | cut -f 1-5,9-17 >> Hela.gene.eRNA.counts
tail -n +2 Hela.eRNA.raw.txt | cut -f 1-5,9-17 >> Hela.gene.eRNA.counts
cut -f 1,6-7,12-13 Hela.gene.eRNA.counts > rep12.Hela.gene.eRNA.counts

Rscript $mydeseq rep12.Hela.gene.eRNA.counts siCTL siCTL siTIP60 siTIP60
# Rscript $mydeseq rep12.Hela.gene.counts siCTL siCTL siTIP60 siTIP60 > rep12.Hela.gene.counts.log
# Rscript $mydeseq rep12.Hela.eRNA.counts siCTL siCTL siTIP60 siTIP60 > rep12.Hela.eRNA.counts.log

# awk -vOFS="\t" '$23 < 0.01 {print $1,$14,$15,$16,$17,$19,$22,$23}' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv 
# awk '$23 < 0.01 && $19 > 1 && $18 > 10' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv > rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv
# version 1
# awk '$23 < 0.01 && $19 > 1' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv > rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv
# awk '$23 < 0.01 && $19 < -1' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv > rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv

# version 2
# awk '$23 < 0.05 && $19 > 1.58' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv > rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv
# awk '$23 < 0.05 && $19 < -1.58' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv > rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv

# version 3
awk '$23 < 0.01 && $19 > 1.58' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv > rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv
awk '$23 < 0.01 && $19 < -1.58' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv > rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv
perl $myperl <(awk -vOFS="\t" '{print $2,$3,$4,$1,".",$5}' Hela.gene.eRNA.counts | tail -n +2) <(cut -f 1 rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv | grep genic_) 3 0 | cut -f 2- | sort -k6,6 -k1,4 > up.eRNA.bed
perl $myperl <(awk -vOFS="\t" '{print $2,$3,$4,$1,".",$5}' Hela.gene.eRNA.counts | tail -n +2) <(cut -f 1 rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv | grep -v genic_) 3 0 | cut -f 2- | sort -k6,6 -k1,4 > up.gene.bed
perl $myperl <(awk -vOFS="\t" '{print $2,$3,$4,$1,".",$5}' Hela.gene.eRNA.counts | tail -n +2) <(cut -f 1 rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv | grep genic_) 3 0 | cut -f 2- | sort -k6,6 -k1,4 > down.eRNA.bed
perl $myperl <(awk -vOFS="\t" '{print $2,$3,$4,$1,".",$5}' Hela.gene.eRNA.counts | tail -n +2) <(cut -f 1 rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv | grep -v genic_) 3 0 | cut -f 2- | sort -k6,6 -k1,4 > down.gene.bed

bedtools intersect -wa -a <(cat int??genic.enhancer.p300.peaks | cut -f 1-6) -b up.eRNA.bed | sort | uniq > up.eRNA.p300.peaks
bedtools intersect -wa -a <(cat int??genic.enhancer.p300.peaks | cut -f 1-6) -b down.eRNA.bed | sort | uniq > down.eRNA.p300.peaks
awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3}print $1,a,b,$4,$5,$6}' up.gene.bed > up.tss.bed
awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3}print $1,a,b,$4,$5,$6}' down.gene.bed > down.tss.bed
awk -vOFS="\t" '{a=$3-1;b=$3;if($6=="-"){a=$2-1;b=$2}print $1,a,b,$4,$5,$6}' up.gene.bed > up.tes.bed
awk -vOFS="\t" '{a=$3-1;b=$3;if($6=="-"){a=$2-1;b=$2}print $1,a,b,$4,$5,$6}' down.gene.bed > down.tes.bed

labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.eRNA.p300.peaks \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig \
-out Hela.DESeq2res.up.eRNA.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m Hela.DESeq2res.up.eRNA.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 10 --outFileSortedRegions Hela.DESeq2res.up.eRNA.srt.bed --sortUsingSamples 6 7 --zMax 20
plotProfile -m Hela.DESeq2res.up.eRNA.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 10 --plotHeight 20
plotProfile -m Hela.DESeq2res.up.eRNA.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --perGroup 
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.eRNA.p300.peaks \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig \
-out Hela.DESeq2res.up.eRNA.GRO.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m Hela.DESeq2res.up.eRNA.GRO.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.GRO.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 10 --outFileSortedRegions Hela.DESeq2res.up.eRNA.GRO.srt.bed
plotProfile -m Hela.DESeq2res.up.eRNA.GRO.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.GRO.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 10 --plotHeight 20
plotProfile -m Hela.DESeq2res.up.eRNA.GRO.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.GRO.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup 
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R Hela.DESeq2res.up.eRNA.GRO.srt.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
-out Hela.DESeq2res.up.eRNA.ChIP.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m Hela.DESeq2res.up.eRNA.ChIP.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.ChIP.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 --heatmapWidth 10 --heatmapHeight 10 --sortRegions keep
plotProfile -m Hela.DESeq2res.up.eRNA.ChIP.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.ChIP.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 --plotWidth 10 --plotHeight 20
plotProfile -m Hela.DESeq2res.up.eRNA.ChIP.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.ChIP.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 --plotWidth 20 --plotHeight 10 --perGroup
labs=`ls Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.eRNA.p300.peaks \
-S Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig \
-out Hela.DESeq2res.up.eRNA.GROsep.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 6
plotHeatmap -m Hela.DESeq2res.up.eRNA.GROsep.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.GROsep.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 10 --zMax 8
plotHeatmap -m Hela.DESeq2res.up.eRNA.GROsep.computeMatrix.gz -o Hela.DESeq2res.up.eRNA.GROsep.pG.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 --heatmapHeight 30 --perGroup --zMax 8


labs=`ls Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.tss.bed up.tes.bed \
-S Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig \
-out Hela.DESeq2res.up.gene.GRO.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m Hela.DESeq2res.up.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.up.gene.GRO.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 10 --outFileSortedRegions Hela.DESeq2res.up.gene.GRO.srt.bed
plotHeatmap -m Hela.DESeq2res.up.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.up.gene.GRO.pG.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 30 --perGroup --zMax 8
plotProfile -m Hela.DESeq2res.up.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.up.gene.GRO.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 10 --plotHeight 20
plotProfile -m Hela.DESeq2res.up.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.up.gene.GRO.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup 
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R down.tss.bed down.tes.bed \
-S Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig \
-out Hela.DESeq2res.down.gene.GRO.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m Hela.DESeq2res.down.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.down.gene.GRO.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 10 --outFileSortedRegions Hela.DESeq2res.down.gene.GRO.srt.bed
plotHeatmap -m Hela.DESeq2res.down.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.down.gene.GRO.pG.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 30 --perGroup 
plotProfile -m Hela.DESeq2res.down.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.down.gene.GRO.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 10 --plotHeight 20
plotProfile -m Hela.DESeq2res.down.gene.GRO.computeMatrix.gz -o Hela.DESeq2res.down.gene.GRO.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup 





labs=`ls Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.eRNA.p300.peaks up.tss.bed up.tes.bed down.eRNA.p300.peaks down.tss.bed down.tes.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
-out Hela.DESeq2res.ChIP.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 270
plotHeatmap -m Hela.DESeq2res.ChIP.computeMatrix.gz -o Hela.DESeq2res.ChIP.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 --heatmapWidth 10 --outFileSortedRegions Hela.DESeq2res.ChIP.srt.bed
plotHeatmap -m Hela.DESeq2res.ChIP.computeMatrix.gz -o Hela.DESeq2res.ChIP.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.DESeq2res.ChIP.computeMatrix.gz -o Hela.DESeq2res.ChIP.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 --plotWidth 20 --plotHeight 10
plotProfile -m Hela.DESeq2res.ChIP.computeMatrix.gz -o Hela.DESeq2res.ChIP.pergp.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 --plotWidth 20 --plotHeight 10 --perGroup

labs=`ls Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.eRNA.p300.peaks up.tss.bed up.tes.bed down.eRNA.p300.peaks down.tss.bed down.tes.bed \
-S Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig \
-out Hela.DESeq2res.GRO.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 270
plotHeatmap -m Hela.DESeq2res.GRO.computeMatrix.gz -o Hela.DESeq2res.GRO.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --outFileSortedRegions Hela.DESeq2res.GRO.srt.bed
plotHeatmap -m Hela.DESeq2res.GRO.computeMatrix.gz -o Hela.DESeq2res.GRO.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.DESeq2res.GRO.computeMatrix.gz -o Hela.DESeq2res.GRO.pergp.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10
plotProfile -m Hela.DESeq2res.GRO.computeMatrix.gz -o Hela.DESeq2res.GRO.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup


labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.eRNA.p300.peaks up.tss.bed up.tes.bed down.eRNA.p300.peaks down.tss.bed down.tes.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig \
-out Hela.DESeq2res.GRO.comb.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 270
plotHeatmap -m Hela.DESeq2res.GRO.comb.computeMatrix.gz -o Hela.DESeq2res.GRO.comb.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --outFileSortedRegions Hela.DESeq2res.GRO.comb.srt.bed
plotHeatmap -m Hela.DESeq2res.GRO.comb.computeMatrix.gz -o Hela.DESeq2res.GRO.comb.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.DESeq2res.GRO.comb.computeMatrix.gz -o Hela.DESeq2res.GRO.comb.pergp.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10
plotProfile -m Hela.DESeq2res.GRO.comb.computeMatrix.gz -o Hela.DESeq2res.GRO.comb.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup

##### Feng's DEG
awk '$6=="+"' Increase_genes_siTIP60.bed | bedtools sort -i - > srt.Increase_genes_siTIP60.bed
awk '$6=="-"' Increase_genes_siTIP60.bed | bedtools sort -i - >> srt.Increase_genes_siTIP60.bed
awk '$6=="+"' Decrease_genes_siTIP60.bed | bedtools sort -i - > srt.Decrease_genes_siTIP60.bed
awk '$6=="-"' Decrease_genes_siTIP60.bed | bedtools sort -i - >> srt.Decrease_genes_siTIP60.bed
labs=`ls Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix scale-regions --missingDataAsZero \
-R srt.Increase_genes_siTIP60.bed srt.Decrease_genes_siTIP60.bed \
-S Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig \
-out Hela.srt.FengDEG.GRO.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 271 --regionBodyLength 5000
plotHeatmap -m Hela.srt.FengDEG.GRO.computeMatrix.gz -o Hela.srt.FengDEG.GRO.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --sortRegions keep

awk '$6=="+"' Increase_genes_siTIP60.bed | bedtools sort -i - > srt.pos.Increase_genes_siTIP60.bed
awk '$6=="-"' Increase_genes_siTIP60.bed | bedtools sort -i - > srt.neg.Increase_genes_siTIP60.bed
awk '$6=="+"' Decrease_genes_siTIP60.bed | bedtools sort -i - > srt.pos.Decrease_genes_siTIP60.bed
awk '$6=="-"' Decrease_genes_siTIP60.bed | bedtools sort -i - > srt.neg.Decrease_genes_siTIP60.bed
labs1=`ls Hela-si{CTL,TIP60}-rpt{1,2}-GROseqpos.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
labs2=`ls Hela-si{CTL,TIP60}-rpt{1,2}-GROseqneg.ucsc.bigWig | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix scale-regions --missingDataAsZero \
-R srt.pos.Increase_genes_siTIP60.bed srt.pos.Decrease_genes_siTIP60.bed \
srt.neg.Increase_genes_siTIP60.bed srt.neg.Decrease_genes_siTIP60.bed \
-S Hela-si{CTL,TIP60}-rpt{1,2}-GROseqpos.ucsc.bigWig Hela-si{CTL,TIP60}-rpt{1,2}-GROseqneg.ucsc.bigWig \
-out Hela.srt.sep.FengDEG.GRO.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 271 --regionBodyLength 5000
plotHeatmap -m Hela.srt.sep.FengDEG.GRO.computeMatrix.gz -o Hela.srt.sep.FengDEG.GRO.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs1 $labs2 --heatmapWidth 10
plotHeatmap -m Hela.srt.sep.FengDEG.GRO.computeMatrix.gz -o Hela.srt.sep.FengDEG.GRO.keep.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs1 $labs2 --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.srt.sep.FengDEG.GRO.computeMatrix.gz -o Hela.srt.sep.FengDEG.GRO.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs1 $labs2 --plotWidth 20 --plotHeight 10 --perGroup

computeMatrix scale-regions --missingDataAsZero \
-R Increase_genes_siTIP60.bed Decrease_genes_siTIP60.bed \
-S Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig \
-out Hela.FengDEG.GRO.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 6 --regionBodyLength 5000
plotHeatmap -m Hela.FengDEG.GRO.computeMatrix.gz -o Hela.FengDEG.GRO.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10
plotHeatmap -m Hela.FengDEG.GRO.computeMatrix.gz -o Hela.FengDEG.GRO.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.FengDEG.GRO.computeMatrix.gz -o Hela.FengDEG.GRO.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup

labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix scale-regions --missingDataAsZero \
-R Increase_genes_siTIP60.bed Decrease_genes_siTIP60.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig \
-out Hela.FengDEG.GRO.comb.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 6 --regionBodyLength 5000
plotHeatmap -m Hela.FengDEG.GRO.comb.computeMatrix.gz -o Hela.FengDEG.GRO.comb.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10
plotHeatmap -m Hela.FengDEG.GRO.comb.computeMatrix.gz -o Hela.FengDEG.GRO.comb.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.FengDEG.GRO.comb.computeMatrix.gz -o Hela.FengDEG.GRO.comb.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup

awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3}print $1,a,b,NR,$5,$6}' Increase_genes_siTIP60.bed > Feng.up.tss.bed
awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3}print $1,a,b,NR,$5,$6}' Decrease_genes_siTIP60.bed > Feng.down.tss.bed
awk -vOFS="\t" '{a=$3-1;b=$3;if($6=="-"){a=$2-1;b=$2}print $1,a,b,NR,$5,$6}' Increase_genes_siTIP60.bed > Feng.up.tes.bed
awk -vOFS="\t" '{a=$3-1;b=$3;if($6=="-"){a=$2-1;b=$2}print $1,a,b,NR,$5,$6}' Decrease_genes_siTIP60.bed > Feng.down.tes.bed
labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R Feng.up.tss.bed Feng.down.tss.bed Feng.up.tes.bed Feng.down.tes.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt{1,2}-GROseq.ucsc.bigWig \
-out Hela.FengDEG.GRO.tsstes.comb.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m Hela.FengDEG.GRO.tsstes.comb.computeMatrix.gz -o Hela.FengDEG.GRO.tsstes.comb.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10
plotHeatmap -m Hela.FengDEG.GRO.tsstes.comb.computeMatrix.gz -o Hela.FengDEG.GRO.tsstes.comb.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.FengDEG.GRO.tsstes.comb.computeMatrix.gz -o Hela.FengDEG.GRO.tsstes.comb.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup

labs=`ls Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R Feng.up.tss.bed Feng.down.tss.bed Feng.up.tes.bed Feng.down.tes.bed \
-S Hela-si{CTL,TIP60}-rpt{1,2}-*.ucsc.bigWig \
-out Hela.FengDEG.GRO.tsstes.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m Hela.FengDEG.GRO.tsstes.computeMatrix.gz -o Hela.FengDEG.GRO.tsstes.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10
plotHeatmap -m Hela.FengDEG.GRO.tsstes.computeMatrix.gz -o Hela.FengDEG.GRO.tsstes.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --sortRegions keep
plotProfile -m Hela.FengDEG.GRO.tsstes.computeMatrix.gz -o Hela.FengDEG.GRO.tsstes.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --perGroup
##### Feng's DEG


###########################################################
# do annatation
# For GSEA: Signal2Noise (default) uses the difference of means scaled by the standard deviation. Note: You must have at least three samples for each phenotype to use this metric.
# cut -f 1,8 Hela.gene.raw.txt | tail -n +2 | cut -f 1 -d"|" > homer.id2gene.txt
# awk '$23 < 0.01 || NR == 1' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv | cut -f 1,10-13| grep -v "genic_" > rld.q0.01.txt
# awk '$19!="NA"' rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv | cut -f 1,10-13 | grep -v "genic_" > rld.Hela.gene.txt
head -1 ../rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv | cut -f 1,10-13 > rld.Hela.gene.up.txt
head -1 ../rep12.Hela.gene.eRNA.counts.DESeq2.out.tsv | cut -f 1,10-13 > rld.Hela.gene.down.txt
cut -f 1,10-13 rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv | grep -v "genic_" >> rld.Hela.gene.up.txt
cut -f 1,10-13 rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv | grep -v "genic_" >> rld.Hela.gene.down.txt
# cut -f 1,10-13 rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv | grep -v "genic_" | perl $myperl ../homer.id2gene.txt /dev/stdin 0 0 | awk -vOFS="\t" '{print $7,$2,$3,$4,$5}' >> rld.Hela.gene.up.txt
# cut -f 1,10-13 rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv | grep -v "genic_" | perl $myperl ../homer.id2gene.txt /dev/stdin 0 0 | awk -vOFS="\t" '{print $7,$2,$3,$4,$5}' >> rld.Hela.gene.down.txt
# cut -f 1,10-13 rep12.Hela.gene.eRNA.counts.DESeq2.up.tsv | perl $myperl <(sort refGene.sim.txt | uniq) /dev/stdin 0 0 | awk -vOFS="\t" '{print $7,$2,$3,$4,$5}' >> rld.Hela.gene.up.txt
# cut -f 1,10-13 rep12.Hela.gene.eRNA.counts.DESeq2.down.tsv | perl $myperl <(sort refGene.sim.txt | uniq) /dev/stdin 0 0 | awk -vOFS="\t" '{print $7,$2,$3,$4,$5}' >> rld.Hela.gene.down.txt

cut -f 4 up.gene.bed > deseq2.eRNA.gene.q0.01.fc2.up.gene.txt
cut -f 4 down.gene.bed > deseq2.eRNA.gene.q0.01.fc2.down.gene.txt
for f in *.gene.txt
do
	perl $myperl <(sort refGene.sim.txt | uniq) $f 0 0 | cut -f 1,3 | perl $myperl Cosmic.CancerGeneCensus.all.gene.anno /dev/stdin 0 1 > $f.CancerGeneCensus
done

for f in *.gene.txt.CancerGeneCensus
do
	pre=`echo $f | cut -d"." -f7`
	all=`cat $pre.gene.bed | wc -l`
	oncogene=`grep -c "oncogene" $f`
	tsg=`grep -c "TSG" $f`
	echo $pre" "$all" "$oncogene" "$tsg
done
# down 671 16 23
# up 534 11 7