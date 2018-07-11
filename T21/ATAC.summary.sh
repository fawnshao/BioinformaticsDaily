#!/bash/bash
# rclone sync -L /home1/04935/shaojf/scratch/T21/ mygoogle:T21/
mydeseq=/home1/04935/shaojf/myTools/BioinformaticsDaily/referenceCode/runDESeq2.R
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

mergePeaks *T21C1*peaks.narrowPeak > NSC.T21C1_peaks.merged
# 	Max distance to merge: direct overlap required (-d given)
# 	Merging peaks... 
# 	Comparing NS13-Xiaoyu-NSC-T21C1_P9-350k_peaks.narrowPeak (55057 total) and NS13-Xiaoyu-NSC-T21C1_P9-350k_peaks.narrowPeak (55057 total)
# 	Comparing NS13-Xiaoyu-NSC-T21C1_P9-350k_peaks.narrowPeak (55057 total) and NS14-Xiaoyu-NSC-T21C1_P8-50k_peaks.narrowPeak (82683 total)
# 	Comparing NS14-Xiaoyu-NSC-T21C1_P8-50k_peaks.narrowPeak (82683 total) and NS13-Xiaoyu-NSC-T21C1_P9-350k_peaks.narrowPeak (55057 total)
# 	Comparing NS14-Xiaoyu-NSC-T21C1_P8-50k_peaks.narrowPeak (82683 total) and NS14-Xiaoyu-NSC-T21C1_P8-50k_peaks.narrowPeak (82683 total)

# NS13-Xiaoyu-NSC-T21C1_P9-350k_peaks.narrowPeak	NS14-Xiaoyu-NSC-T21C1_P8-50k_peaks.narrowPeak	Total	Name
# 	X	32506	NS14-Xiaoyu-NSC-T21C1_P8-50k_peaks.narrowPeak
# X		3097	NS13-Xiaoyu-NSC-T21C1_P9-350k_peaks.narrowPeak
# X	X	50069	NS13-Xiaoyu-NSC-T21C1_P9-350k_peaks.narrowPeak|NS14-Xiaoyu-NSC-T21C1_P8-50k_peaks.narrowPeak
mergePeaks *T21C5*peaks.narrowPeak > NSC.T21C5_peaks.merged
# 	Max distance to merge: direct overlap required (-d given)
# 	Merging peaks... 
# 	Comparing NS13-Xiaoyu-NSC-T21C5_P4-500k_peaks.narrowPeak (45483 total) and NS13-Xiaoyu-NSC-T21C5_P4-500k_peaks.narrowPeak (45483 total)
# 	Comparing NS13-Xiaoyu-NSC-T21C5_P4-500k_peaks.narrowPeak (45483 total) and NS14-Xiaoyu-NSC-T21C5_P4-50k_peaks.narrowPeak (11156 total)
# 	Comparing NS14-Xiaoyu-NSC-T21C5_P4-50k_peaks.narrowPeak (11156 total) and NS13-Xiaoyu-NSC-T21C5_P4-500k_peaks.narrowPeak (45483 total)
# 	Comparing NS14-Xiaoyu-NSC-T21C5_P4-50k_peaks.narrowPeak (11156 total) and NS14-Xiaoyu-NSC-T21C5_P4-50k_peaks.narrowPeak (11156 total)

# NS13-Xiaoyu-NSC-T21C5_P4-500k_peaks.narrowPeak	NS14-Xiaoyu-NSC-T21C5_P4-50k_peaks.narrowPeak	Total	Name
# 	X	993	NS14-Xiaoyu-NSC-T21C5_P4-50k_peaks.narrowPeak
# X		35666	NS13-Xiaoyu-NSC-T21C5_P4-500k_peaks.narrowPeak
# X	X	9778	NS13-Xiaoyu-NSC-T21C5_P4-500k_peaks.narrowPeak|NS14-Xiaoyu-NSC-T21C5_P4-50k_peaks.narrowPeak

awk '$8==2{print $2,$3,$4,$5,$6}' NSC.T21C1_peaks.merged | awk -vOFS="\t" '{print $1,$2,$3,"NSC.T21C1."NR,$5,$4}' > common.NSC.T21C1_peaks
awk '$8==2{print $2,$3,$4,$5,$6}' NSC.T21C5_peaks.merged | awk -vOFS="\t" '{print $1,$2,$3,"NSC.T21C4."NR,$5,$4}' > common.NSC.T21C5_peaks
mergePeaks common.NSC.T21C1_peaks common.NSC.T21C5_peaks > common.NSC.T21C1-5_peaks.merged
# 	Max distance to merge: direct overlap required (-d given)
# 	Merging peaks... 
# 	Comparing common.NSC.T21C1_peaks (48143 total) and common.NSC.T21C1_peaks (48143 total)
# 	Comparing common.NSC.T21C1_peaks (48143 total) and common.NSC.T21C5_peaks (9368 total)
# 	Comparing common.NSC.T21C5_peaks (9368 total) and common.NSC.T21C1_peaks (48143 total)
# 	Comparing common.NSC.T21C5_peaks (9368 total) and common.NSC.T21C5_peaks (9368 total)

# common.NSC.T21C1_peaks	common.NSC.T21C5_peaks	Total	Name
# 	X	1041	common.NSC.T21C5_peaks
# X		39812	common.NSC.T21C1_peaks
# X	X	8315	common.NSC.T21C1_peaks|common.NSC.T21C5_peaks
Rscript ~/stampede2/myTools/BioinformaticsDaily/RVisualization/vennplot_from_a_stdin.R NSC.T21C5 NSC.T21C1 1041 39812 8315

mergePeaks NSC.T21C1_peaks.merged NSC.T21C5_peaks.merged > all.NSC.T21C1-5_peaks.merged
# 	Max distance to merge: direct overlap required (-d given)
# 	Merging peaks... 
# 	Comparing NSC.T21C1_peaks.merged (85672 total) and NSC.T21C1_peaks.merged (85672 total)
# 	Comparing NSC.T21C1_peaks.merged (85672 total) and NSC.T21C5_peaks.merged (46437 total)
# 	Comparing NSC.T21C5_peaks.merged (46437 total) and NSC.T21C1_peaks.merged (85672 total)
# 	Comparing NSC.T21C5_peaks.merged (46437 total) and NSC.T21C5_peaks.merged (46437 total)

# NSC.T21C1_peaks.merged	NSC.T21C5_peaks.merged	Total	Name
# 	X	3860	NSC.T21C5_peaks.merged
# X		44553	NSC.T21C1_peaks.merged
# X	X	40871	NSC.T21C1_peaks.merged|NSC.T21C5_peaks.merged
Rscript ~/stampede2/myTools/BioinformaticsDaily/RVisualization/vennplot_from_a_stdin.R NSC.T21C5.all NSC.T21C1.all 3860 44553 40871

analyzeRepeats.pl all.NSC.T21C1-5_peaks.merged hg19 -d nochrM.*.mTD -strand both -raw > all.NSC.T21C1-5_peaks.raw.txt
echo "ID NS13-Xiaoyu-NSC-T21C1_P9-350k NS13-Xiaoyu-NSC-T21C5_P4-500k NS14-Xiaoyu-NSC-T21C1_P8-50k NS14-Xiaoyu-NSC-T21C5_P4-50k" | tr " " "\t" > all.NSC.T21C1-5_peaks.sim.count.txt
cut -f 1,9-12 all.NSC.T21C1-5_peaks.raw.txt | tail -n +2 >> all.NSC.T21C1-5_peaks.sim.count.txt
Rscript $mydeseq all.NSC.T21C1-5_peaks.sim.count.txt T21C1 T21C5 T21C1 T21C5

awk '$23 < 0.01 && ($19 > 1 || $19 < -1)' all.NSC.T21C1-5_peaks.sim.count.txt.DESeq2.out.tsv > all.NSC.T21C1-5_peaks.de.txt
perl $myperl <(cut -f 1-4 all.NSC.T21C1-5_peaks.merged | tail -n +2) <(cut -f 1 all.NSC.T21C1-5_peaks.de.txt) 0 0 | awk -vOFS="\t" '{print $3,$4,$5,$1}' > all.NSC.T21C1-5_peaks.de.bed

labs="NSC.T21C1 NSC.T21C5 NSC.T21C1 NSC.T21C5"
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R all.NSC.T21C1-5_peaks.de.bed \
-S *.bigWig \
-out all.NSC.T21C1-5_peaks.de.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotHeatmap -m all.NSC.T21C1-5_peaks.de.computeMatrix.gz -o all.NSC.T21C1-5_peaks.de.srt.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10
plotProfile -m all.NSC.T21C1-5_peaks.de.computeMatrix.gz -o all.NSC.T21C1-5_peaks.de.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10
plotHeatmap -m all.NSC.T21C1-5_peaks.de.computeMatrix.gz -o all.NSC.T21C1-5_peaks.de.srt.k3.plotHeatmap.pdf --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --kmeans 3 --outFileSortedRegions all.NSC.T21C1-5_peaks.de.srt.k3.bed
plotProfile -m all.NSC.T21C1-5_peaks.de.computeMatrix.gz -o all.NSC.T21C1-5_peaks.de.k3.plotProfile.pdf --plotFileFormat pdf --samplesLabel $labs --plotWidth 20 --plotHeight 10 --kmeans 3
