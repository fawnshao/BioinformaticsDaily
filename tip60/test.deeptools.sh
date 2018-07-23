pre=all.enhancer
labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R enhancer.p300.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotProfile -m $pre.computeMatrix.gz -o $pre.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 10 --plotHeight 20
plotProfile -m $pre.computeMatrix.gz -o $pre.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --perGroup 
plotHeatmap -m $pre.computeMatrix.gz -o $pre.plotHeatmap.srtbysiTIP.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.srt.bed --sortUsingSamples 5 --zMax 20 40 20 3 3
plotHeatmap -m $pre.computeMatrix.gz -o $pre.plotHeatmap.srtbyp300.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.srt.bed --sortUsingSamples 1 --zMax 20 40 20 3 3
plotHeatmap -m $pre.computeMatrix.gz -o $pre.plotHeatmap.k3.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.k3.bed --kmeans 3 --zMax 20 40 20 3 3



pre=all.enhancer
labs=`ls Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R enhancer.p300.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig \
-out $pre.GROsep.computeMatrix.gz \
-b 3000 -a 3000 -bs 10 -p 272
plotProfile -m $pre.GROsep.computeMatrix.gz -o $pre.GROsep.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 10 --plotHeight 20
plotProfile -m $pre.GROsep.computeMatrix.gz -o $pre.GROsep.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --perGroup 
plotHeatmap -m $pre.GROsep.computeMatrix.gz -o $pre.GROsep.plotHeatmap.srtbyGRO.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.srt4567.bed --sortUsingSamples 4 5 6 7 --zMax 20 40 20 2 2 2 2
plotHeatmap -m $pre.GROsep.computeMatrix.gz -o $pre.GROsep.plotHeatmap.srtbyp300.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.srt.bed --sortUsingSamples 1 --zMax 20 40 20 3 3 3 3
plotHeatmap -m $pre.GROsep.computeMatrix.gz -o $pre.GROsep.plotHeatmap.k3.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.k3.bed --kmeans 3 --zMax 20 40 20 3 3 3 3

pre=all.enhancer
labs=`ls Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R enhancer.p300.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig \
-out $pre.GROsep.bs50.computeMatrix.gz \
-b 3000 -a 3000 -bs 50 -p 272
plotProfile -m $pre.GROsep.bs50.computeMatrix.gz -o $pre.GROsep.bs50.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 10 --plotHeight 20
plotProfile -m $pre.GROsep.bs50.computeMatrix.gz -o $pre.GROsep.bs50.k3.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 10 --plotHeight 20 --kmeans 3 
plotProfile -m $pre.GROsep.bs50.computeMatrix.gz -o $pre.GROsep.bs50.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --perGroup 
plotHeatmap -m $pre.GROsep.bs50.computeMatrix.gz -o $pre.GROsep.bs50.plotHeatmap.srtbyGRO.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs50.srt4567.bed --sortUsingSamples 4 5 6 7 --zMax 20 40 20 2 2 2 2
plotHeatmap -m $pre.GROsep.bs50.computeMatrix.gz -o $pre.GROsep.bs50.plotHeatmap.srtbyp300.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs50.srt.bed --sortUsingSamples 1 --zMax 20 40 20 2 2 2 2
plotHeatmap -m $pre.GROsep.bs50.computeMatrix.gz -o $pre.GROsep.bs50.plotHeatmap.k3.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs50.k3.bed --kmeans 3 --zMax 20 40 20 2 2 2 2


pre=all.enhancer
labs=`ls Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R enhancer.p300.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig \
-out $pre.GROsep.bs400.computeMatrix.gz \
-b 10000 -a 10000 -bs 400 -p 272
plotProfile -m $pre.GROsep.bs400.computeMatrix.gz -o $pre.GROsep.bs400.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 10 --plotHeight 20
plotProfile -m $pre.GROsep.bs400.computeMatrix.gz -o $pre.GROsep.bs400.k3.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 10 --plotHeight 20 --kmeans 3 
plotProfile -m $pre.GROsep.bs400.computeMatrix.gz -o $pre.GROsep.bs400.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --perGroup 
plotHeatmap -m $pre.GROsep.bs400.computeMatrix.gz -o $pre.GROsep.bs400.plotHeatmap.srtbyGRO.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs400.srt4567.bed --sortUsingSamples 4 5 6 7 --zMax 20 40 20 2 2 2 2
plotHeatmap -m $pre.GROsep.bs400.computeMatrix.gz -o $pre.GROsep.bs400.plotHeatmap.srtbyp300.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs400.srt.bed --sortUsingSamples 1 --zMax 20 40 20 2 2 2 2
plotHeatmap -m $pre.GROsep.bs400.computeMatrix.gz -o $pre.GROsep.bs400.plotHeatmap.k3.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs400.k3.bed --kmeans 3 --zMax 20 40 20 2 2 2 2

pre=all.enhancer
labs=`ls Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R enhancer.p300.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
Hela-si{CTL,TIP60}-rpt1-GROseq{pos,neg}.ucsc.bigWig \
-out $pre.GROsep.bs5.computeMatrix.gz \
-b 5000 -a 5000 -bs 5 -p 270
plotProfile -m $pre.GROsep.bs5.computeMatrix.gz -o $pre.GROsep.bs5.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --yMax 15 20 10 1.5 1.5 1.5 1.5
plotProfile -m $pre.GROsep.bs5.computeMatrix.gz -o $pre.GROsep.bs5.k3.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --kmeans 3 --yMax 20 70 20 3 3 3 3
plotProfile -m $pre.GROsep.bs5.computeMatrix.gz -o $pre.GROsep.bs5.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --perGroup
plotHeatmap -m $pre.GROsep.bs5.computeMatrix.gz -o $pre.GROsep.bs5.plotHeatmap.srtbyGRO.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs5.srt4567.bed --sortUsingSamples 4 5 6 7 --zMax 10 20 10 1 1 1 1 --yMax 10 20 10 1 1 1 1
plotHeatmap -m $pre.GROsep.bs5.computeMatrix.gz -o $pre.GROsep.bs5.plotHeatmap.srtbyp300.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs5.srt.bed --sortUsingSamples 1 --zMax 10 20 10 1 1 1 1 --yMax 10 20 10 1 1 1 1
plotHeatmap -m $pre.GROsep.bs5.computeMatrix.gz -o $pre.GROsep.bs5.plotHeatmap.k3.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROsep.bs5.k3.bed --kmeans 3 --zMax 10 20 10 3 3 3 3 --yMax 20 70 20 3 3 3 3


pre=all.enhancer
labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R enhancer.p300.bed \
-S hg19.MichaelSnyder.p300.rep12.ENCFF289TVJ.bigWig \
hg19.BradleyBernstein.H3K27ac.rep12.ENCFF388WMD.bigWig \
hg19.RichardMyers.pol2.rep12.ENCFF959MZN.bigWig \
../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.GROcomb.bs5.computeMatrix.gz \
-b 5000 -a 5000 -bs 5 -p 270
plotProfile -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --yMax 15 20 10 1.5 1.5
plotProfile -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.k3.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --kmeans 3 --yMax 20 70 20 3 3
plotProfile -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.pG.plotProfile.pdf --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --plotWidth 20 --plotHeight 10 --perGroup
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.srtbyGRO.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROcomb.bs5.srt45.bed --sortUsingSamples 4 5 --zMax 10 20 10 1 1 --yMax 10 20 10 1 1
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.srtbyp300.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROcomb.bs5.srt.bed --sortUsingSamples 1 --zMax 10 20 10 1 1 --yMax 10 20 10 1 1
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.k3.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel p300 H3K27ac Pol2 $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROcomb.bs5.k3.bed --kmeans 3 --zMax 10 20 10 1 1 --yMax 20 70 20 3 3






awk -vOFS="\t" -F"\t" '$9+$15>1{print $2,$3,$4,$1,$6,$5}' Hela.gene.fpkm.txt | tail -n +2 > Hela.expressed.gene.bed
awk -vOFS="\t" -F"\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3}print $1,a,b,$4,$5,$6}' Hela.expressed.gene.bed > Hela.expressed.tss.bed
awk -vOFS="\t" -F"\t" '{a=$2-1;b=$2;if($6=="+"){a=$3-1;b=$3}print $1,a,b,$4,$5,$6}' Hela.expressed.gene.bed > Hela.expressed.tes.bed

pre=all.tss
labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R Hela.expressed.tss.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.GROcomb.bs5.computeMatrix.gz \
-b 5000 -a 5000 -bs 5 -p 272
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.k5.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROcomb.bs5.k5.bed --kmeans 5 --zMax 10 10 &
pre=all.tes
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R Hela.expressed.tes.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.GROcomb.bs5.computeMatrix.gz \
-b 5000 -a 5000 -bs 5 -p 272
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.k5.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROcomb.bs5.k5.bed --kmeans 5 --zMax 10 10 &

perl $myperl <(cat up.tss.bed down.tss.bed | cut -f 4) Hela.expressed.tss.bed 0 3 | grep "/" | cut -f 1-6 > other.tss.bed
perl $myperl <(cat up.tes.bed down.tes.bed | cut -f 4) Hela.expressed.tes.bed 0 3 | grep "/" | cut -f 1-6 > other.tes.bed
labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
pre=up.down.other.tss
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.tss.bed down.tss.bed other.tss.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.GROcomb.bs5.computeMatrix.gz \
-b 5000 -a 5000 -bs 5 -p 272
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROcomb.bs5.bed --zMax 10 10 &
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.1.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 30 --outFileSortedRegions $pre.GROcomb.bs5.bed --zMax 5 5 &

pre=up.down.other.tes
computeMatrix reference-point --referencePoint center --missingDataAsZero \
-R up.tss.bed down.tss.bed other.tss.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.GROcomb.bs5.computeMatrix.gz \
-b 5000 -a 5000 -bs 5 -p 272
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 20 --outFileSortedRegions $pre.GROcomb.bs5.bed --zMax 10 10 &
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.1.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 10 --heatmapHeight 30 --outFileSortedRegions $pre.GROcomb.bs5.bed --zMax 5 5 &

pre=up.down.other.tss
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.1.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 --heatmapHeight 30 --outFileSortedRegions $pre.GROcomb.bs5.bed --zMax 5 5 &
pre=up.down.other.tes
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.1.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 --heatmapHeight 30 --outFileSortedRegions $pre.GROcomb.bs5.bed --zMax 5 5 &


perl $myperl <(cat up.gene.bed down.gene.bed | cut -f 4) Hela.expressed.gene.bed 0 3 | grep "/" | cut -f 1-6 > other.gene.bed
labs=`ls ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig | awk -F"/" '{print $NF}' | sed 's/Hela-//g;s/.ucsc.bigWig//g'`
pre=up.down.other.gene
computeMatrix scale-regions --regionBodyLength 3000 --missingDataAsZero \
-R up.gene.bed down.gene.bed other.gene.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.GROcomb.bs5.computeMatrix.gz \
-b 5000 -a 5000 -bs 5 -p 272 
plotHeatmap -m $pre.GROcomb.bs5.computeMatrix.gz -o $pre.GROcomb.bs5.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 --heatmapHeight 30 --outFileSortedRegions $pre.GROcomb.bs5.bed --zMax 5 5 &
computeMatrix scale-regions --regionBodyLength 5000 --missingDataAsZero \
-R up.gene.bed down.gene.bed other.gene.bed \
-S ../Feng.data/bigwig/Hela-si{CTL,TIP60}-rpt1-GROseq.ucsc.bigWig \
-out $pre.GROcomb.bs5.5k.computeMatrix.gz \
-b 1000 -a 1000 -bs 5 -p 272 
plotHeatmap -m $pre.GROcomb.bs5.5k.computeMatrix.gz -o $pre.GROcomb.bs5.5k.plotHeatmap.pdf --colorMap Greens --whatToShow 'plot, heatmap and colorbar' --plotFileFormat pdf --samplesLabel $labs --heatmapWidth 20 --heatmapHeight 30 --outFileSortedRegions $pre.GROcomb.bs5.5k.bed --zMax 5 5 &
