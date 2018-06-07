#!/bin/sh
#### NUP tRNA
# rclone sync /home1/04935/shaojf/stampede2/NUP_project/NUP_integrate/ mygoogle:NUP_project/NUP_integrate/
# multiBigwigSummary bins --binSize 1000 --bwfiles bigwigs/*ucsc.bigWig other.bw/*ucsc.bigWig --outFileName MCF7-NUP.bin1000.out.npz --outRawCounts MCF7-NUP.bin1000.out.coverage
# plotCorrelation -in MCF7-NUP.bin1000.out.npz -o MCF7-NUP.bin1000.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf --labels ATAC.shNUP53-2_Dox_E2 ATAC.shNUP53-2_E2 ATAC.shNUP93-1_Dox_E2 ATAC.shNUP93-1_E2 ChIP.nup53blrp_EtOH ChIP.nup53blrp_E2 ChIP.Nup93scbt_EtOH ChIP.Nup93scbt_E2 ChIP.Nup93bethyl_EtOH ChIP.Nup93bethyl_E2 ChIP.GRO.siCTL-E2-GROneg GRO.siCTL-E2-GROpos GRO.siNUP53-E2-GROneg GRO.siNUP53-E2-GROpos GRO.siNUP93-E2-GROneg GRO.siNUP93-E2-GROpos

# multiBigwigSummary bins --binSize 1000 --bwfiles bigwigs/*ucsc.bigWig --outFileName NS9-Xiaoyu-MCF7-NUP.ATAC.bin1000.out.npz --outRawCounts NS9-Xiaoyu-MCF7-NUP.ATAC.bin1000.out.coverage
# plotCorrelation -in NS9-Xiaoyu-MCF7-NUP.ATAC.bin1000.out.npz -o NS9-Xiaoyu-MCF7-NUP.ATAC.bin1000.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf

awk -vOFS="\t" '{print $1,$2-1000,$3+1000,$4,$5,$6}' hg19-tRNAs.bed > hg19-tRNAs.flank.bed 

multiBigwigSummary BED-file --BED hg19-tRNAs.flank.bed --bwfiles bigwigs/*ucsc.bigWig other.bw/*ucsc.bigWig --outFileName Xiaoyu-MCF7-NUP.tRNA.out.npz --outRawCounts Xiaoyu-MCF7-NUP.tRNA.out.coverage
plotCorrelation -in Xiaoyu-MCF7-NUP.tRNA.out.npz -o Xiaoyu-MCF7-NUP.tRNA.out.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotFileFormat pdf

# computeMatrix reference-point -R tmp.bed -S bigwigs/*ucsc.bigWig -out NS9-Xiaoyu-MCF7-NUP.ATAC.Nup53blrp_E2_peaks.computeMatrix.gz --referencePoint center -b 3000 -a 3000 -bs 10 
# plotHeatmap -m NS9-Xiaoyu-MCF7-NUP.ATAC.Nup53blrp_E2_peaks.computeMatrix.gz -o NS9-Xiaoyu-MCF7-NUP.ATAC.Nup53blrp_E2_peaks.computeMatrix.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel shNUP53-2_Dox_E2 shNUP53-2_E2 shNUP93-1_Dox_E2 shNUP93-1_E2 --heatmapWidth 10

computeMatrix reference-point -R hg19-tRNAs.bed -S bigwigs/*ucsc.bigWig other.bw/*ucsc.bigWig -out Xiaoyu-MCF7-NUP.tRNA.out.1k.gz --referencePoint center -b 1000 -a 1000 -bs 10 
plotHeatmap -m Xiaoyu-MCF7-NUP.tRNA.out.1k.gz -o Xiaoyu-MCF7-NUP.tRNA.out.1k.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel ATAC.shNUP53-2_Dox_E2 ATAC.shNUP53-2_E2 ATAC.shNUP93-1_Dox_E2 ATAC.shNUP93-1_E2 ChIP.NUP53.EtOH ChIP.NUP53.E2 ChIP.NUP93scbt.EtOH ChIP.NUP93scbt.E2 ChIP.NUP93bethyl.EtOH ChIP.NUP93bethyl.E2 shCTL-GRO-E2neg shCTL-GRO-E2pos shNUP53-2-GRO-E2neg shNUP53-2-GRO-E2pos shNUP53-3-GRO-E2neg shNUP53-3-GRO-E2pos siCTL-E2-GROneg siCTL-E2-GROpos siCTL-ETOH-GROneg siCTL-ETOH-GROpos siNUP53-E2-GROneg siNUP53-E2-GROpos siNUP53-ETOH-GROneg siNUP53-ETOH-GROpos siNUP93-E2-GROneg siNUP93-E2-GROpos siNUP93-ETOH-GROneg siNUP93-ETOH-GROpos --heatmapWidth 10

computeMatrix reference-point -R hg19-tRNAs.bed -S bigwigs/*ucsc.bigWig other.bw/*ucsc.bigWig -out Xiaoyu-MCF7-NUP.tRNA.out.3k.gz --referencePoint center -b 3000 -a 3000 -bs 10 
plotHeatmap -m Xiaoyu-MCF7-NUP.tRNA.out.3k.gz -o Xiaoyu-MCF7-NUP.tRNA.out.3k.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel ATAC.shNUP53-2_Dox_E2 ATAC.shNUP53-2_E2 ATAC.shNUP93-1_Dox_E2 ATAC.shNUP93-1_E2 ChIP.NUP53.EtOH ChIP.NUP53.E2 ChIP.NUP93scbt.EtOH ChIP.NUP93scbt.E2 ChIP.NUP93bethyl.EtOH ChIP.NUP93bethyl.E2 shCTL-GRO-E2neg shCTL-GRO-E2pos shNUP53-2-GRO-E2neg shNUP53-2-GRO-E2pos shNUP53-3-GRO-E2neg shNUP53-3-GRO-E2pos siCTL-E2-GROneg siCTL-E2-GROpos siCTL-ETOH-GROneg siCTL-ETOH-GROpos siNUP53-E2-GROneg siNUP53-E2-GROpos siNUP53-ETOH-GROneg siNUP53-ETOH-GROpos siNUP93-E2-GROneg siNUP93-E2-GROpos siNUP93-ETOH-GROneg siNUP93-ETOH-GROpos --heatmapWidth 10

cut -f 2-4 promoters.ann.txt > hg19.promoters.bed
computeMatrix reference-point -R hg19.promoters.bed -S bigwigs/*ucsc.bigWig other.bw/*ucsc.bigWig -out Xiaoyu-MCF7-NUP.promoters.out.3k.gz --referencePoint center -b 3000 -a 3000 -bs 10 
plotHeatmap -m Xiaoyu-MCF7-NUP.promoters.out.3k.gz -o Xiaoyu-MCF7-NUP.promoters.out.3k.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel ATAC.shNUP53-2_Dox_E2 ATAC.shNUP53-2_E2 ATAC.shNUP93-1_Dox_E2 ATAC.shNUP93-1_E2 ChIP.NUP53.EtOH ChIP.NUP53.E2 ChIP.NUP93scbt.EtOH ChIP.NUP93scbt.E2 ChIP.NUP93bethyl.EtOH ChIP.NUP93bethyl.E2 shCTL-GRO-E2neg shCTL-GRO-E2pos shNUP53-2-GRO-E2neg shNUP53-2-GRO-E2pos shNUP53-3-GRO-E2neg shNUP53-3-GRO-E2pos siCTL-E2-GROneg siCTL-E2-GROpos siCTL-ETOH-GROneg siCTL-ETOH-GROpos siNUP53-E2-GROneg siNUP53-E2-GROpos siNUP53-ETOH-GROneg siNUP53-ETOH-GROpos siNUP93-E2-GROneg siNUP93-E2-GROpos siNUP93-ETOH-GROneg siNUP93-ETOH-GROpos --heatmapWidth 10
# Traceback (most recent call last):
#   File "/home1/04935/shaojf/.local/bin/plotHeatmap", line 11, in <module>
#     main(args)
#   File "/home1/04935/shaojf/.local/lib/python2.7/site-packages/deeptools/plotHeatmap.py", line 605, in main
#     hm.read_matrix_file(matrix_file)
#   File "/home1/04935/shaojf/.local/lib/python2.7/site-packages/deeptools/heatmapper.py", line 747, in read_matrix_file
#     for line in fh:
#   File "/opt/apps/intel17/python/2.7.13/lib/python2.7/gzip.py", line 464, in readline
#     c = self.read(readsize)
#   File "/opt/apps/intel17/python/2.7.13/lib/python2.7/gzip.py", line 268, in read
#     self._read(readsize)
#   File "/opt/apps/intel17/python/2.7.13/lib/python2.7/gzip.py", line 315, in _read
#     self._read_eof()
#   File "/opt/apps/intel17/python/2.7.13/lib/python2.7/gzip.py", line 354, in _read_eof
#     hex(self.crc)))
# IOError: CRC check failed 0x5cb7814 != 0xf166d567L

computeMatrix reference-point -R hg19.promoters.bed -S other.bw/MCF7-sh* -out Xiaoyu-MCF7-NUP.sh.promoters.out.3k.gz --referencePoint center -b 3000 -a 3000 -bs 10 
plotHeatmap -m Xiaoyu-MCF7-NUP.sh.promoters.out.3k.gz -o Xiaoyu-MCF7-NUP.sh.promoters.out.3k.plotHeatmap.pdf --kmeans 3 --colorMap Greens --whatToShow 'heatmap and colorbar' --plotFileFormat pdf --samplesLabel shCTL-GRO-E2neg shCTL-GRO-E2pos shNUP53-2-GRO-E2neg shNUP53-2-GRO-E2pos shNUP53-3-GRO-E2neg shNUP53-3-GRO-E2pos --heatmapWidth 10

# plotProfile -m NS9-Xiaoyu-MCF7-NUP.ATAC.Nup53blrp_E2_peaks.computeMatrix.gz -o NS9-Xiaoyu-MCF7-NUP.ATAC.Nup53blrp_E2_peaks.computeMatrix.plotProfile.pdf --kmeans 1 --plotFileFormat pdf --samplesLabel shNUP53-2_Dox_E2 shNUP53-2_E2 shNUP93-1_Dox_E2 shNUP93-1_E2 --plotWidth 10 --plotHeight 10 --refPointLabel PeakCenter --plotType fill
# plotProfile -m NS9-Xiaoyu-MCF7-NUP.ATAC.Nup53blrp_E2_peaks.computeMatrix.1k.gz -o NS9-Xiaoyu-MCF7-NUP.ATAC.Nup53blrp_E2_peaks.computeMatrix.1k.plotProfile.pdf --kmeans 1 --plotFileFormat pdf --samplesLabel shNUP53-2_Dox_E2 shNUP53-2_E2 shNUP93-1_Dox_E2 shNUP93-1_E2 --plotWidth 10 --plotHeight 10 --refPointLabel PeakCenter --plotType fill

annotatePeaks.pl <(awk -vOFS="\t" '{print $1,$2}')