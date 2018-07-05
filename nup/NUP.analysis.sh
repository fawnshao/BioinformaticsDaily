#!/bin/bash
# rclone sync -L /home1/04935/shaojf/scratch/NUP_related_ChIP/201807.analysis mygoogle:NUP_project/201807.analysis

### deeptools for all ChIP-Seq
# enhancer region: MCF-7_p300.andMACS2.nopromoter.bed
# TSS region: hg19.refGene.tss.uniq.srt.bed
# NUP53 peak: Nup53blrp_E2_peaks.narrowPeak
# boundary region
# p300 peak region 
workdir=/home1/04935/shaojf/scratch/NUP_related_ChIP/201807.analysis
visdir=$workdir/all.bigwigs
name=all.nup
# makeMultiWigHub.pl $name hg19 -url $visdir -webdir $visdir -d *.mTD &

labels=`ls $visdir/all.nup/hg19/*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
multiBigwigSummary bins --binSize 1000 --bwfiles $visdir/all.nup/hg19/*mTD.ucsc.bigWig --outFileName MCF7-NUP.1k.out.npz --outRawCounts MCF7-NUP.1k.out.coverage --labels $labels
plotCorrelation -in MCF7-NUP.1k.out.npz -o MCF7-NUP.1k.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --colorMap bwr

# plotCorrelation -in MCF7-NUP.E2.1k.out.npz -o MCF7-NUP.E2.1k.bwr.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --colorMap bwr &
# plotCorrelation -in MCF7-NUP.E2.5k.out.npz -o MCF7-NUP.E2.5k.bwr.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --colorMap bwr &
# --colorMap RdYlBu

computeMatrix reference-point -R MCF-7_p300.andMACS2.nopromoter.bed -S $visdir/all.nup/hg19/*mTD.ucsc.bigWig -out MCF-7_p300.andMACS2.nopromoter.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 
plotHeatmap -m MCF-7_p300.andMACS2.nopromoter.out.gz -o MCF-7_p300.andMACS2.nopromoter.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labels --heatmapWidth 10
plotProfile -m MCF-7_p300.andMACS2.nopromoter.out.gz -o MCF-7_p300.andMACS2.nopromoter.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel $labels --numPlotsPerRow 3

computeMatrix reference-point -R hg19.refGene.tss.uniq.srt.bed -S $visdir/all.nup/hg19/*mTD.ucsc.bigWig -out hg19.refGene.tss.uniq.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 
plotHeatmap -m hg19.refGene.tss.uniq.out.gz -o hg19.refGene.tss.uniq.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labels --heatmapWidth 10
plotProfile -m hg19.refGene.tss.uniq.out.gz -o hg19.refGene.tss.uniq.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel $labels --numPlotsPerRow 3

cut -f 1-4 Nup53blrp_E2_peaks.narrowPeak > Nup53blrp_E2.bed
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/*mTD.ucsc.bigWig -out Nup53blrp_E2.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 
plotHeatmap -m Nup53blrp_E2.out.gz -o Nup53blrp_E2.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labels --heatmapWidth 10 --outFileSortedRegions Nup53blrp_E2.k3.sortedRegions.bed
plotProfile -m Nup53blrp_E2.out.gz -o Nup53blrp_E2.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel $labels --numPlotsPerRow 3
plotHeatmap -m Nup53blrp_E2.out.gz -o Nup53blrp_E2.k5.plotHeatmap.pdf --kmeans 5 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labels --heatmapWidth 10 --outFileSortedRegions Nup53blrp_E2.k5.sortedRegions.bed
# awk -vOFS="\t" '{print $1,int(($2+$3)/2-1),int(($2+$3)/2),$4}' Nup53blrp_E2.bed > Nup53blrp_E2.center.bed
labs=`ls $visdir/all.nup/hg19/H3*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/ChIP217*mTD.ucsc.bigWig $visdir/all.nup/hg19/ChIP216*mTD.ucsc.bigWig $visdir/all.nup/hg19/H3*mTD.ucsc.bigWig -out Nup53blrp_E2.H3.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 -p 10
plotHeatmap -m Nup53blrp_E2.H3.out.gz -o Nup53blrp_E2.H3.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --heatmapWidth 10 --refPointLabel "peak center" --outFileSortedRegions Nup53blrp_E2.H3.k3.sortedRegions.bed &
plotProfile -m Nup53blrp_E2.H3.out.gz -o Nup53blrp_E2.H3.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --numPlotsPerRow 3 &

labs=`ls $visdir/all.nup/hg19/{HECTD1,NCAPG}*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/ChIP217*mTD.ucsc.bigWig $visdir/all.nup/hg19/ChIP216*mTD.ucsc.bigWig $visdir/all.nup/hg19/{HECTD1,NCAPG}*mTD.ucsc.bigWig -out Nup53blrp_E2.HECTD1_NCAPG.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 -p 10
plotHeatmap -m Nup53blrp_E2.HECTD1_NCAPG.out.gz -o Nup53blrp_E2.HECTD1_NCAPG.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --heatmapWidth 10 --refPointLabel "peak center" --outFileSortedRegions Nup53blrp_E2.HECTD1_NCAPG.k3.sortedRegions.bed &
plotProfile -m Nup53blrp_E2.HECTD1_NCAPG.out.gz -o Nup53blrp_E2.HECTD1_NCAPG.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --numPlotsPerRow 3 &

labs=`ls $visdir/all.nup/hg19/{p300,POLII}*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/ChIP217*mTD.ucsc.bigWig $visdir/all.nup/hg19/ChIP216*mTD.ucsc.bigWig $visdir/all.nup/hg19/{p300,POLII}*mTD.ucsc.bigWig -out Nup53blrp_E2.p300_POLII.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 -p 10
plotHeatmap -m Nup53blrp_E2.p300_POLII.out.gz -o Nup53blrp_E2.p300_POLII.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --heatmapWidth 10 --refPointLabel "peak center" --outFileSortedRegions Nup53blrp_E2.p300_POLII.k3.sortedRegions.bed &
plotProfile -m Nup53blrp_E2.p300_POLII.out.gz -o Nup53blrp_E2.p300_POLII.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --numPlotsPerRow 3 &

labs=`ls $visdir/all.nup/hg19/{CTCF,RAD21}*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/ChIP217*mTD.ucsc.bigWig $visdir/all.nup/hg19/ChIP216*mTD.ucsc.bigWig $visdir/all.nup/hg19/{CTCF,RAD21}*mTD.ucsc.bigWig -out Nup53blrp_E2.CTCF_RAD21.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 -p 10
plotHeatmap -m Nup53blrp_E2.CTCF_RAD21.out.gz -o Nup53blrp_E2.CTCF_RAD21.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --heatmapWidth 10 --refPointLabel "peak center" --outFileSortedRegions Nup53blrp_E2.CTCF_RAD21.k3.sortedRegions.bed &
plotProfile -m Nup53blrp_E2.CTCF_RAD21.out.gz -o Nup53blrp_E2.CTCF_RAD21.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --numPlotsPerRow 3 &

labs=`ls $visdir/all.nup/hg19/{ERa,FoxA1,GATA3}*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/ChIP217*mTD.ucsc.bigWig $visdir/all.nup/hg19/ChIP216*mTD.ucsc.bigWig $visdir/all.nup/hg19/{ERa,FoxA1,GATA3}*mTD.ucsc.bigWig -out Nup53blrp_E2.ERa_FoxA1_GATA3.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 -p 10
plotHeatmap -m Nup53blrp_E2.ERa_FoxA1_GATA3.out.gz -o Nup53blrp_E2.ERa_FoxA1_GATA3.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --heatmapWidth 10 --refPointLabel "peak center" --outFileSortedRegions Nup53blrp_E2.ERa_FoxA1_GATA3.k3.sortedRegions.bed &
plotProfile -m Nup53blrp_E2.ERa_FoxA1_GATA3.out.gz -o Nup53blrp_E2.ERa_FoxA1_GATA3.plotProfile.pdf --plotType=fill --perGroup --kmeans 3 --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs --numPlotsPerRow 3 &


multiBigwigSummary BED-file --BED Nup53blrp_E2.bed --bwfiles $visdir/all.nup/hg19/*mTD.ucsc.bigWig --outFileName MCF7-NUP.Nup53blrp_E2.out.npz --outRawCounts MCF7-NUP.Nup53blrp_E2.out.coverage --labels $labels -p 10
plotCorrelation -in MCF7-NUP.Nup53blrp_E2.out.npz -o MCF7-NUP.Nup53blrp_E2.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf

###############################################################
# NUPChIP p300 POLII ERa K27ac NUPATAC
# mergePeaks cutsites.shNUP*_peaks.narrowPeak > MCF7_shNUP.ATAC.merged.peaks
# cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak	cutsites.shNUP53-2_E2_peaks.narrowPeak	cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak	cutsites.shNUP93-1_E2_peaks.narrowPeak	Total	Name
# 			X	4423	cutsites.shNUP93-1_E2_peaks.narrowPeak
# 		X		717	cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak
# 		X	X	501	cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak|cutsites.shNUP93-1_E2_peaks.narrowPeak
# 	X			6992	cutsites.shNUP53-2_E2_peaks.narrowPeak
# 	X		X	5018	cutsites.shNUP53-2_E2_peaks.narrowPeak|cutsites.shNUP93-1_E2_peaks.narrowPeak
# 	X	X		550	cutsites.shNUP53-2_E2_peaks.narrowPeak|cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak
# 	X	X	X	1847	cutsites.shNUP53-2_E2_peaks.narrowPeak|cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak|cutsites.shNUP93-1_E2_peaks.narrowPeak
# X				1647	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak
# X			X	1059	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak|cutsites.shNUP93-1_E2_peaks.narrowPeak
# X		X		151	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak|cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak
# X		X	X	396	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak|cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak|cutsites.shNUP93-1_E2_peaks.narrowPeak
# X	X			1805	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak|cutsites.shNUP53-2_E2_peaks.narrowPeak
# X	X		X	6812	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak|cutsites.shNUP53-2_E2_peaks.narrowPeak|cutsites.shNUP93-1_E2_peaks.narrowPeak
# X	X	X		527	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak|cutsites.shNUP53-2_E2_peaks.narrowPeak|cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak
# X	X	X	X	37463	cutsites.shNUP53-2_Dox_E2_peaks.narrowPeak|cutsites.shNUP53-2_E2_peaks.narrowPeak|cutsites.shNUP93-1_Dox_E2_peaks.narrowPeak|cutsites.shNUP93-1_E2_peaks.narrowPeak
# visdir=./all.bigwigs
labs1=`ls $visdir/all.nup/hg19/{ERa,p300,POLII}*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
labs2=`ls $visdir/all.nup/hg19/H3K27ac*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
labs3=`ls $visdir/bigwigs/*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/NS9-Xiaoyu-MCF7-//;s/_ATAC[1-4]_S[1-4]/_ATAC/;s/.mTD.ucsc.bigWig//' | tr "\n" " "`
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/ChIP217*mTD.ucsc.bigWig $visdir/all.nup/hg19/ChIP216*mTD.ucsc.bigWig $visdir/all.nup/hg19/{ERa,p300,POLII}*mTD.ucsc.bigWig $visdir/all.nup/hg19/H3K27ac*mTD.ucsc.bigWig $visdir/bigwigs/*mTD.ucsc.bigWig -out Nup53blrp_E2.ERa_p300_POLII_H3K27ac_ATAC.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 -p 48
plotHeatmap -m Nup53blrp_E2.ERa_p300_POLII_H3K27ac_ATAC.out.gz -o Nup53blrp_E2.ERa_p300_POLII_H3K27ac_ATAC.plotHeatmap.pdf --kmeans 5 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs1 $labs2 $labs3 --heatmapWidth 10 --refPointLabel "peak center" --outFileSortedRegions Nup53blrp_E2.ERa_p300_POLII_H3K27ac_ATAC.k5.sortedRegions.bed &

labs1=`ls $visdir/all.nup/hg19/{p300,POLII}*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
labs2=`ls $visdir/all.nup/hg19/H3K27ac*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/.mTD.ucsc.bigWig//' | tr "\n" " "`
labs3=`ls $visdir/bigwigs/*mTD.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/NS9-Xiaoyu-MCF7-//;s/_ATAC[1-4]_S[1-4]/_ATAC/;s/.mTD.ucsc.bigWig//' | tr "\n" " "`
computeMatrix reference-point -R Nup53blrp_E2.bed -S $visdir/all.nup/hg19/ChIP217*mTD.ucsc.bigWig $visdir/all.nup/hg19/ChIP216*mTD.ucsc.bigWig $visdir/all.nup/hg19/{p300,POLII}*mTD.ucsc.bigWig $visdir/all.nup/hg19/H3K27ac*mTD.ucsc.bigWig $visdir/bigwigs/*mTD.ucsc.bigWig -out Nup53blrp_E2.p300_POLII_H3K27ac_ATAC.out.gz --referencePoint center -b 1000 -a 1000 -bs 10 -p 48
plotHeatmap -m Nup53blrp_E2.p300_POLII_H3K27ac_ATAC.out.gz -o Nup53blrp_E2.p300_POLII_H3K27ac_ATAC.plotHeatmap.pdf --kmeans 4 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel Nup53blrp_E2 Nup53blrp_veh $labs1 $labs2 $labs3 --heatmapWidth 10 --refPointLabel "peak center" --outFileSortedRegions Nup53blrp_E2.p300_POLII_H3K27ac_ATAC.k4.sortedRegions.bed &

###############################################################
awk '$NF=="cluster_1"' Nup53blrp_E2.p300_POLII.k2.sortedRegions.bed | cut -f 1-4 > Nup53blrp_E2.p300_POLII.k2.cluster_1.bed
awk '$NF=="cluster_2"' Nup53blrp_E2.p300_POLII.k2.sortedRegions.bed | cut -f 1-4 > Nup53blrp_E2.p300_POLII.k2.cluster_2.bed
 

awk '$NF=="cluster_1"' Nup53blrp_E2.k3.sortedRegions.bed | cut -f 1-4 > Nup53blrp_E2.k3.cluster_1.bed
awk '$NF=="cluster_2"' Nup53blrp_E2.k3.sortedRegions.bed | cut -f 1-4 > Nup53blrp_E2.k3.cluster_2.bed
awk '$NF=="cluster_3"' Nup53blrp_E2.k3.sortedRegions.bed | cut -f 1-4 > Nup53blrp_E2.k3.cluster_3.bed
