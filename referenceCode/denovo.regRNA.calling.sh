#!/bin/bash
# common files
workdir=/public/workspace/shaojf/shaojf/mydatasets/Lingqian/test.pipeline
hisat2ref=/public/workspace/shaojf/shaojf/mydatasets/Lingqian/reference/hg19
genegtf=/public/workspace/shaojf/shaojf/mydatasets/Lingqian/reference/gencode.v32lift37.annotation.sorted.gtf
dhsfile=/public/workspace/shaojf/shaojf/mydatasets/Lingqian/ENCODE/K562/hg19.K562.DNase.hotspots.ENCFF635VMY.bed
blackListFileName=/public/workspace/shaojf/shaojf/mydatasets/Lingqian/ENCODE/hg19.ENCFF001TDO.merged.bed

# input fastq
prefix=SRR5243202
sample=K562.TT

# common parameters
Threads=16

# enviroment
bytlib load bedtools-2.25.0

#######################################################
# skipped hisat2 steps
# samtools view -1 -q 10 -bo $BamDir/$i.q10.bam --threads $Threads $BamDir/$i.Aligned.out.bam
# samtools sort --threads 10 $BamDir/$i.q10.bam -o $BamDir/$i.sorted.bam
# samtools index $BamDir/$i.sorted.bam
#######################################################


#######################################################
# stringtie step
stringtie $prefix.hisat2.bam --rf -G $genegtf \
	-o $prefix.stringtie.gtf -A $prefix.stringtie.abundance.txt \
	-l $sample -m 500 -g 200 -p $Threads
bedtools intersect -v -s \
	-a <(awk '$3=="transcript"' $prefix.stringtie.gtf) \
	-b <(awk '$1 !~ /GL/ && $3=="transcript"' $genegtf) > $prefix.stringtie.novel.gtf
Rscript ????
#######################################################
# K562.TT.8 novel eRNA?
# awk '$3=="transcript" && $0 !~ /reference_id/' $prefix.stringtie.gtf
# filter by quantile?


#######################################################
# could manually rename the files
# --filterRNAstrand {forward,reverse} This option assumes a standard dUTP-based library 
# 	preparation (that is, --filterRNAstrand=forward keeps minus-strand reads, 
# 		which originally came from genes on the forward strand using a dUTP-based method).
# deeptools step
bamCoverage -b $prefix.hisat2.bam -o CPM.$prefix.plus.bw -of bigwig -bs 1 -p $Threads \
	--blackListFileName $blackListFileName --filterRNAstrand forward \
	--effectiveGenomeSize 2864785220 --normalizeUsing CPM
bamCoverage -b $prefix.hisat2.bam -o CPM.$prefix.minus.bw -of bigwig -bs 1 -p Threads \
	--blackListFileName $blackListFileName --filterRNAstrand reverse \
	--effectiveGenomeSize 2864785220 --normalizeUsing CPM
#######################################################


# bigWig.jar  compbio.jar  compbioLib.jar  log4j-1.2.15.jar
java -Xmx48000m -cp compbio.jar:log4j-1.2.15.jar:compbioLib.jar:bigWig.jar scripts.lincs.AnalyzeEnhancerRNAs analyze_dhs CPM.K562.GRO NONE hg19.K562.DNase.hotspots.ENCFF635VMY.nonPromoter.DHS.bed -minSum 2 hg19.chrom.sizes CPM.K562.GRO



########################################################
bytlib load bedtools-2.25.0
bedtools merge -s -d 300 -c 4,4 -o count_distinct,distinct -i <(cat K562.*HStools.bed | bedtools sort -i) > merged.K562.HStools.bed 
awk -v OFS="\t" '$6~/TT/ && $6!~/PRO/ && $6!~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.HStools.bed > K562.HStools.TT.u.bed
awk -v OFS="\t" '$6!~/TT/ && $6~/PRO/ && $6!~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.HStools.bed > K562.HStools.PRO.u.bed
awk -v OFS="\t" '$6!~/TT/ && $6!~/PRO/ && $6~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.HStools.bed > K562.HStools.GRO.u.bed
awk -v OFS="\t" '$6~/TT/ && $6~/PRO/ && $6~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.HStools.bed > K562.HStools.all.bed
# wc -l K562.HStools.*
  #   11 K562.HStools.all.bed
  #  944 K562.HStools.GRO.u.bed
  #  686 K562.HStools.PRO.u.bed
  # 6459 K562.HStools.TT.u.bed
  # 8100 total

# awk -v OFS="\t" '$4=="+"{print $1,$2-1,$2,$6,$5,$4}' merged.K562.HStools.bed > fwd.merged.K562.HStools.bed
# awk -v OFS="\t" '$4=="-"{print $1,$3-1,$3,$6,$5,$4}' merged.K562.HStools.bed > rev.merged.K562.HStools.bed

# for f in fwd.merged.K562.HStools.bed
# do
#     computeMatrix reference-point \
#         -R $f \
#         -S CPM.K562.*.plus.bw \
#         --blackListFileName hg19.ENCFF001TDO.merged.bed \
#         -out $f.computeMatrix.gz \
#         --missingDataAsZero \
#         -b 3000 -a 3000 -bs 5 -p 10 &
# done
# for f in rev.merged.K562.HStools.bed
# do
#     computeMatrix reference-point \
#         -R $f \
#         -S CPM.K562.*.minus.bw \
#         --blackListFileName hg19.ENCFF001TDO.merged.bed \
#         -out $f.computeMatrix.gz \
#         --missingDataAsZero \
#         -b 3000 -a 3000 -bs 5 -p 10 &
# done

# for f in fwd.merged.K562.HStools.bed
# do
#     computeMatrix reference-point \
#         -R $f \
#         -S CPM.K562.*.minus.bw \
#         --blackListFileName hg19.ENCFF001TDO.merged.bed \
#         -out antisense.$f.computeMatrix.gz \
#         --missingDataAsZero \
#         -b 3000 -a 3000 -bs 5 -p 10 &
# done
# for f in rev.merged.K562.HStools.bed
# do
#     computeMatrix reference-point \
#         -R $f \
#         -S CPM.K562.*.plus.bw \
#         --blackListFileName hg19.ENCFF001TDO.merged.bed \
#         -out antisense.$f.computeMatrix.gz \
#         --missingDataAsZero \
#         -b 3000 -a 3000 -bs 5 -p 10 &
# done


# cat <(zcat fwd.merged.K562.HStools.bed.computeMatrix.gz) \
#     <(zcat rev.merged.K562.HStools.bed.computeMatrix.gz | tail -n +2) > sense.merged.K562.HStools.computeMatrix
# cat <(zcat antisense.fwd.merged.K562.HStools.bed.computeMatrix.gz) \
#     <(zcat antisense.rev.merged.K562.HStools.bed.computeMatrix.gz | tail -n +2) > antisense.merged.K562.HStools.computeMatrix

#      # 15596 antisense.merged.K562.HStools.computeMatrix
#      # 15596 sense.merged.K562.HStools.computeMatrix
# # vim computeMatrix

# gzip antisense.merged.K562.HStools.computeMatrix 
# gzip sense.merged.K562.HStools.computeMatrix 

# plotHeatmap -m sense.merged.K562.HStools.computeMatrix.gz \
# 	-o plotHeatmap.sense.merged.K562.HStools.pdf \
# 	--outFileSortedRegions plotHeatmap.srt.sense.merged.K562.HStools.bed \
# 	--colorMap Oranges --heatmapHeight 16 --heatmapWidth 8  \
# 	--yMax 0.3 0.3 0.3 0.2 0.2 --zMax 0.15 0.15 0.15 0.1 0.1 \
# 	-x "HStools ncRNA" --samplesLabel GRO PRO.1 PRO.2 TT.1 TT.2
# zcat antisense.merged.K562.HStools.computeMatrix.gz | head -1 > srt.antisense.merged.K562.HStools.computeMatrix
# perl /public/workspace/shaojf/shaojf/mytools/add_any_2files_together.pl <(zcat antisense.merged.K562.HStools.computeMatrix.gz | tail -n +2) <(tail -n +2 plotHeatmap.srt.sense.merged.K562.HStools.bed | cut -f 4) 3 0 | awk '$2!="-"' | cut -f 2- >> srt.antisense.merged.K562.HStools.computeMatrix
# # 15597
# gzip srt.antisense.merged.K562.HStools.computeMatrix 
# plotHeatmap -m srt.antisense.merged.K562.HStools.computeMatrix.gz \
# 	-o plotHeatmap.antisense.merged.K562.HStools.pdf \
# 	--sortRegions no \
# 	--colorMap Oranges --heatmapHeight 16 --heatmapWidth 8  \
# 	--yMax 1.5 1 1 0.2 0.2 --zMax 1.5 1 1 0.1 0.1 \
# 	-x "HStools ncRNA" --samplesLabel GRO PRO.1 PRO.2 TT.1 TT.2

# tail -n +2 plotHeatmap.srt.sense.merged.K562.HStools.bed | cut -f 1-6 > plotHeatmap.srt.sense.merged.K562.HStools.nohead.bed
# for f in plotHeatmap.srt.sense.merged.K562.HStools.nohead.bed
# do
#     computeMatrix reference-point \
#         -R $f \
#         -S K562/Histone/hg19.K562.H*.bigWig \
#         --blackListFileName hg19.ENCFF001TDO.merged.bed \
#         -out hist.$f.computeMatrix.gz \
#         --missingDataAsZero \
#         -b 3000 -a 3000 -bs 5 -p 12 &
# done
# labs=`ll K562/Histone/hg19.K562.H*.bigWig | awk -F"." '{print $3}'`
# plotHeatmap -m hist.plotHeatmap.srt.sense.merged.K562.HStools.nohead.bed.computeMatrix.gz \
# 	-o plotHeatmap.hist.merged.K562.HStools.pdf \
# 	--sortRegions no \
# 	--colorMap Blues --heatmapHeight 16 --heatmapWidth 8  \
# 	-x "HStools ncRNA" --samplesLabel $labs


####################
awk -vOFS="\t" '{print $1,$2,$3,$6,$5,$4}' merged.K562.HStools.bed > merged.K562.HStools.arranged.bed
awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3;}print $1,a,b,$4,$5,$6}' merged.K562.HStools.arranged.bed > merged.K562.HStools.arranged.tss.bed

bedtools intersect -u -a merged.K562.HStools.arranged.bed -b <(cat K562/hg19.K562.*hotspots*.bed) > merged.K562.HStools.DHSintersect.full.bed
awk -vOFS="\t" '{a=$2-1;b=$2;if($6=="-"){a=$3-1;b=$3;}print $1,a,b,$4,$5,$6}' merged.K562.HStools.DHSintersect.full.bed > merged.K562.HStools.DHSintersect.bed
# wc -l merged.K562.HStools.arranged.bed
# 8344 merged.K562.HStools.arranged.bed
# wc -l merged.K562.HStools.DHSintersect.bed 
# 5193 merged.K562.HStools.DHSintersect.bed
# bedtools intersect -u -a K562/hg19.K562.DNase.hotspots.ENCFF162UKK.bed \
# 	-b <(awk -vOFS="\t" '{print $1,$2,$2+1"\n"$1,$3,$3+1}' merged.K562.HStools.arranged.bed) > merged.K562.HStools.DHSintersect.bed
# wc -l merged.K562.HStools.DHSintersect.bed 
# 5485 merged.K562.HStools.DHSintersect.bed

# (14:35 stu16230119@server Summary)$ more merged.K562.HStools.DHSintersect.bed | wc -l
# 10049
ln -s merged.K562.HStools.arranged.tss.bed merged.K562.HStools.DHSintersect.bed
for f in merged.K562.HStools.DHSintersect.bed
do
    computeMatrix reference-point \
        -R $f \
        -S K562/hg19.K562.*DNase.*.bigWig \
        --blackListFileName hg19.ENCFF001TDO.merged.bed \
        -out DNase.$f.computeMatrix.gz \
        --missingDataAsZero \
        -b 5000 -a 5000 -bs 5 -p 24
done
plotHeatmap -m DNase.merged.K562.HStools.DHSintersect.bed.computeMatrix.gz \
	-o plotHeatmap.DNase.merged.K562.HStools.DHSintersect.pdf \
	--outFileSortedRegions plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.bed \
	--colorMap Greens --heatmapHeight 16 --heatmapWidth 8  \
	-z "HStools ncRNA" --samplesLabel DNase G1.DNase G2.DNase

tail -n +2 plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.bed | cut -f 1-6 > plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.nohead.bed

# for f in GRO PRO.1 PRO TT.1 TT
# do
#     ./bigWigMerge CPM.K562.$f.minus.bw CPM.K562.$f.plus.bw CPM.K562.$f.bdg &
# done
# for f in GRO PRO.1 PRO TT.1 TT
# do
#     bedtools sort -i CPM.K562.$f.bdg > CPM.K562.$f.srt.bdg 
#     ./bedGraphToBigWig CPM.K562.$f.srt.bdg hg19.chrom.sizes CPM.K562.$f.merged.bw &
# done
# for f in GRO PRO.1 PRO TT.1 TT
# do
#     computeMatrix reference-point \
#         -R plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.nohead.bed \
#         -S CPM.K562.$f.merged.bw \
#         --blackListFileName hg19.ENCFF001TDO.merged.bed \
#         -out DNase.srt.$f.merged.computeMatrix.gz \
#         --missingDataAsZero \
#         -b 3000 -a 3000 -bs 5 -p 20
# done
# for f in GRO PRO.1 PRO TT.1 TT
# do
#     plotHeatmap -m DNase.srt.$f.merged.computeMatrix.gz \
# 		-o plotHeatmap.DNase.srt.$f.merged.pdf \
# 		--sortRegions no \
# 		--colorMap Oranges --heatmapHeight 16 --heatmapWidth 8 \
# 		-z "HStools ncRNA" --samplesLabel $f
# done

# for f in K562/Histone/hg19.K562.H*.bigWig
# do
# 	name=`echo $f | awk -F"." '{print $3"."$4}'`
# 	computeMatrix reference-point \
#         -R plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.nohead.bed \
#         -S $f \
#         --blackListFileName hg19.ENCFF001TDO.merged.bed \
#         -out DNase.srt.$name.computeMatrix.gz \
#         --missingDataAsZero \
#         -b 3000 -a 3000 -bs 5 -p 20
# done
# for f in K562/Histone/hg19.K562.H*.bigWig
# do
# 	name=`echo $f | awk -F"." '{print $3"."$4}'`
# 	plotHeatmap -m DNase.srt.$name.computeMatrix.gz \
# 		-o plotHeatmap.DNase.srt.$name.pdf \
# 		--sortRegions no \
# 		--colorMap Blues --heatmapHeight 16 --heatmapWidth 8 \
# 		-x "HStools ncRNA" --samplesLabel $name
# done


#######
f=plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.nohead.bed
computeMatrix reference-point \
    -R $f \
    -S CPM.K562.*.merged.bw \
    --blackListFileName hg19.ENCFF001TDO.merged.bed \
    -out DNase.srt.K562nascentRNA.computeMatrix.gz \
    --missingDataAsZero \
    -b 5000 -a 5000 -bs 5 -p 24
plotHeatmap -m DNase.srt.K562nascentRNA.computeMatrix.gz \
	-o plotHeatmap.DNase.srt.K562nascentRNA.pdf \
	--sortRegions no \
	--colorMap Oranges --heatmapHeight 16 --heatmapWidth 8  \
	--yMax 4.5 2.1 2.1 0.45 0.45 --zMax 9 4.2 4.2 0.5 0.5 \
	-z "HStools novel ncRNA" --samplesLabel GRO PRO.1 PRO TT.1 TT
plotHeatmap -m DNase.srt.K562nascentRNA.computeMatrix.gz \
	-o plotHeatmap.only.DNase.srt.K562nascentRNA.pdf \
	--sortRegions no \
	--whatToShow "heatmap only" \
	--colorMap Oranges --heatmapHeight 16 --heatmapWidth 8  \
	--yMax 4.5 2.1 2.1 0.45 0.45 --zMax 9 4.2 4.2 0.5 0.5 \
	-z "HStools novel ncRNA" --samplesLabel GRO PRO.1 PRO TT.1 TT

computeMatrix reference-point \
    -R $f \
    -S K562/Histone/hg19.K562.{H3K4me1-human.ENCFF444SGK,H3K4me2-human.ENCFF118MMT,H3K4me3-human.ENCFF879AFU,H3K9ac-human.ENCFF866KTJ,H3K27ac-human.ENCFF010PHG,H3K79me2-human.ENCFF003CLZ}.bigWig \
    --blackListFileName hg19.ENCFF001TDO.merged.bed \
    -out DNase.srt.K562hist.computeMatrix.gz \
    --missingDataAsZero \
    -b 5000 -a 5000 -bs 5 -p 24
plotHeatmap -m DNase.srt.K562hist.computeMatrix.gz \
	-o plotHeatmap.DNase.srt.K562hist.pdf \
	--sortRegions no \
	--colorMap Blues --heatmapHeight 16 --heatmapWidth 8  \
	--yMax 2.5 6 10 4.5 6.5 3 --zMax 5 12 20 9 13 6 \
	-z "HStools novel ncRNA" \
	--samplesLabel H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K27ac H3K79me2

#####
bedtools genomecov -bg -i hg19.mut.srt.bed -g hg19.chrom.sizes > hg19.mut.srt.bdg
bedtools genomecov -bg -i hg19.cnv.srt.bed -g hg19.chrom.sizes > hg19.cnv.srt.bdg
for f in cnv mut
do
    ./bedGraphToBigWig hg19.$f.srt.bdg hg19.chrom.sizes hg19.$f.bw &
done

f=plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.nohead.bed
for bw in hg19.mut.bw hg19.cnv.bw
do
	computeMatrix reference-point \
	    -R $f \
	    -S $bw \
	    --blackListFileName hg19.ENCFF001TDO.merged.bed \
	    -out DNase.srt.$bw.computeMatrix.gz \
	    --missingDataAsZero \
	    -b 5000 -a 5000 -bs 100 -p 24 &
done
for bw in hg19.mut.bw hg19.cnv.bw
do
	plotHeatmap -m DNase.srt.$bw.computeMatrix.gz \
		-o plotHeatmap.DNase.srt.$bw.pdf \
		--outFileSortedRegions plotHeatmap.srt.$bw.bed \
		--colorMap Reds --heatmapHeight 16 --heatmapWidth 8 \
		-y "HStools novel ncRNA" \
		--samplesLabel $bw
done

# bedtools intersect -u -a merged.K562.HStools.bed -b ../ICGC.release28.leukemia/hg19.mut.srt.bed > merged.K562.HStools.ICGCmutintersect.bed
# bedtools intersect -u -a merged.K562.HStools.bed -b ../ICGC.release28.leukemia/hg19.cnv.srt.bed > merged.K562.HStools.ICGCcnvintersect.bed
# wc -l merged.K562.HStools*.bed
#   15663 merged.K562.HStools.bed
#   10049 merged.K562.HStools.DHSintersect.bed
#   15427 merged.K562.HStools.ICGCcnvintersect.bed
#   11650 merged.K562.HStools.ICGCmutintersect.bed
bedtools intersect -wao -a merged.K562.HStools.arranged.bed -b ../ICGC.release28.leukemia/hg19.mut.srt.bed > merged.K562.HStools.ICGCmutintersect.detail.bed
bedtools intersect -wao -a merged.K562.HStools.arranged.bed -b ../ICGC.release28.leukemia/hg19.cnv.srt.bed > merged.K562.HStools.ICGCcnvintersect.detail.bed
awk '$7!="."{print $1,$2,$3,$4,$10}' merged.K562.HStools.ICGCmutintersect.detail.bed | sort -u | awk '{print $1,$2,$3,$4}' | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$4,$5,$1}' | sort -k5,5nr >  merged.K562.HStools.ICGCmutintersect.counts
awk '$7!="."{print $1,$2,$3,$4,$10}' merged.K562.HStools.ICGCcnvintersect.detail.bed | sort -u | awk '{print $1,$2,$3,$4}' | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$4,$5,$1}' | sort -k5,5nr >  merged.K562.HStools.ICGCcnvintersect.counts

head -500 merged.K562.HStools.ICGCmutintersect.counts | cut -f 1-4 > merged.K562.HStools.ICGCmutintersect.top.bed
head -500 merged.K562.HStools.ICGCcnvintersect.counts | cut -f 1-4 > merged.K562.HStools.ICGCcnvintersect.top.bed
head -3000 plotHeatmap.srt.DNase.merged.K562.HStools.DHSintersect.nohead.bed | cut -f 1-3 > top.srt.DNase.merged.K562.HStools.DHSintersect.nohead.bed


#####
bytlib load bedtools-2.25.0
bedtools merge -s -d 300 -c 4,4 -o count_distinct,distinct -i <(cat K562.*EPCtools.bed | bedtools sort -i) > merged.K562.EPCtools.bed 
awk -v OFS="\t" '$6~/TT/ && $6!~/PRO/ && $6!~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.EPCtools.bed > K562.EPCtools.TT.u.bed
awk -v OFS="\t" '$6!~/TT/ && $6~/PRO/ && $6!~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.EPCtools.bed > K562.EPCtools.PRO.u.bed
awk -v OFS="\t" '$6!~/TT/ && $6!~/PRO/ && $6~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.EPCtools.bed > K562.EPCtools.GRO.u.bed
awk -v OFS="\t" '$6~/TT/ && $6~/PRO/ && $6~/GRO/{print $1,$2,$3,$6,$5,$4}' merged.K562.EPCtools.bed > K562.EPCtools.all.bed
# wc -l *K562.*EPCtools.bed
#    7384 K562.GRO.EPCtools.bed
#   11257 K562.PRO.1.EPCtools.bed
#   12249 K562.PRO.EPCtools.bed
#   14474 K562.TT.1.EPCtools.bed
#   13714 K562.TT.EPCtools.bed
#   31266 merged.K562.EPCtools.bed
awk -vOFS="\t" '{print $1,$2,$3,$6,$5,$4}' merged.K562.EPCtools.bed > merged.K562.EPCtools.arranged.bed
mergePeaks -strand -prefix HStools.vs.EPCtools merged.K562.HStools.arranged.bed merged.K562.EPCtools.arranged.bed

