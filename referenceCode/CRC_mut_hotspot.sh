#!/bin/sh
# the genomic feature of mutation hotspot
### add chr or not?
# high expression more mutation?
# mutation in a TAD?

for f in WGS.COAD-US.bed WGS.COCA-CN.bed
do
	pre=`echo $f | sed 's/.bed//'`
	cut -f 4 $f | cut -f 1 -d"|" | sort | uniq -c | awk -vOFS="\t" '{print $2"\t"$1}' | sort -k2,2nr > samplecount.$pre
done

for f in WGS.COAD-US.bed WGS.COCA-CN.bed
do
	pre=`echo $f | sed 's/.bed//'`
	bedtools merge -d 1000 -i $f > bedtools.merge.$pre
done

for f in WGS.COAD-US.bed WGS.COCA-CN.bed
do
	pre=`echo $f | sed 's/.bed//'`
	bedtools cluster -d 1000 -i $f > bedtools.cluster.$pre
done

for f in WGS.COAD-US.bed WGS.COCA-CN.bed
do
	pre=`echo $f | sed 's/.bed//'`
	cut -f5 bedtools.cluster.$pre | uniq -c | awk -vOFS="\t" '{print $2"\t"$1}' | sort -k2,2nr > bedtools.cluster.$pre.counts
done
# cut -f 5 bedtools.cluster.WGS.COAD-US | uniq | wc -l
# 37998
# cut -f 5 bedtools.cluster.WGS.COCA-CN | uniq | wc -l
# 182145
#    56122 bedtools.cluster.WGS.COAD-US
#   521271 bedtools.cluster.WGS.COCA-CN
myperl=add_any_2files_together.pl 
# for f in WGS.COAD-US.bed WGS.COCA-CN.bed
# do
# 	pre=`echo $f | sed 's/.bed//'`
# 	perl $myperl bedtools.cluster.$pre <(awk '$2>5{print $1}' bedtools.cluster.$pre.counts) 4 0 | cut -f 2-5 > gt5.$f
# done
pre=WGS.COAD-US
perl $myperl bedtools.cluster.$pre <(awk '$2>5{print $1}' bedtools.cluster.$pre.counts) 4 0 | cut -f 2-5 > top.$pre.bed
pre=WGS.COCA-CN
perl $myperl bedtools.cluster.$pre <(awk '$2>100{print $1}' bedtools.cluster.$pre.counts) 4 0 | cut -f 2-5 > top.$pre.bed

for f in WGS.COAD-US.bed WGS.COCA-CN.bed
do
	pre=`echo $f | sed 's/.bed//'`
	annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{print "chr"$1,$2,$3,$4}' $f) hg19 -size given 1> $pre.homer.txt 2> $pre.homer.log &
done


##### DNA methylation
# bigWigAverageOverBed in.bw in.bed out.tab
# for f in WGS.COAD-US.bed WGS.COCA-CN.bed
# bedtools random -l 1000 -n 50000 -seed 123456 -g /home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes > hg19.seed123456.bed
for f in WGS.COAD-US.bed WGS.COCA-CN.bed
do
	pre=`echo $f | sed 's/.bed//'`
	for m in *FractionalMethylation.bigwig
	do
		tail=`echo $m | sed 's/_FractionalMethylation.bigwig//'`
		bigWigAverageOverBed $m <(awk -F"\t" -vOFS="\t" '{print "chr"$1,$2,$3,"chr"$1":"$2"-"$3}' $f | sort | uniq) mut.$pre.$tail.tab &
	done
done

for f in bedtools.merge.WGS.*
do
	pre=`echo $f | sed 's/bedtools.merge.//'`
	for m in *FractionalMethylation.bigwig
	do
		tail=`echo $m | sed 's/_FractionalMethylation.bigwig//'`
		bigWigAverageOverBed $m <(awk -F"\t" -vOFS="\t" '{print "chr"$1,$2,$3,"chr"$1":"$2"-"$3}' $f) $pre.$tail.tab &
	done
done

# for f in bedtools.merge.WGS.*
# do
# 	pre=`echo $f | sed 's/bedtools.merge.//'`
# 	bedtools flank -i <(awk -F"\t" -vOFS="\t" '{print "chr"$1,$2,$3,"chr"$1":"$2"-"$3}' $f) -g /home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes -b 1000 > $pre.flank.bed
# 	for m in *FractionalMethylation.bigwig
# 	do
# 		tail=`echo $m | sed 's/_FractionalMethylation.bigwig//'`
# 		bigWigAverageOverBed $m $pre.flank.bed flank.$pre.$tail.tab &
# 	done
# done


for reg in WGS.mut/d100/top.WGS.CO*
do
	pre=`echo $reg | awk -F"/" '{print $NF}' | awk -F"." '{print $(NF-1)}'`
	for f in `ls bigwig.files/ | awk -F"." '{print $3}' | sort | uniq`
	do
		computeMatrix reference-point \
			-R $reg \
			-S bigwig.files/*${f}*.bigwig \
			-p 48 \
			-b 100 -a 100 -bs 5 \
			--referencePoint center \
			--skipZeros -o d100.$pre.${f}_center.gz \
			--outFileNameMatrix d100.$pre.${f}_center.tab \
			--outFileSortedRegions d100.$pre.${f}_center.bed
		labs=`ls bigwig.files/*${f}*.bigwig | awk -F"/" '{print $NF}' | awk -F"." '{print $2"."$3}' | tr "\n" " "`
		plotHeatmap -m d100.$pre.${f}_center.gz \
			-o d100.$pre.${f}_center.plotHeatmap.pdf \
			--colorMap Greens --whatToShow 'heatmap and colorbar' \
			--plotFileFormat pdf --samplesLabel $labs \
			--heatmapWidth 10 \
			--kmeans 3 \
			--outFileSortedRegions d100.$pre.${f}_center.k3.sortedRegions.bed &
		computeMatrix reference-point \
			-R $reg \
			-S bigwig.files/*${f}*.bigwig \
			-p 48 \
			-b 5000 -a 5000 -bs 5 \
			--referencePoint center \
			--skipZeros -o d100.5k.$pre.${f}_center.gz \
			--outFileNameMatrix d100.5k.$pre.${f}_center.tab \
			--outFileSortedRegions d100.5k.$pre.${f}_center.bed
		plotHeatmap -m d100.5k.$pre.${f}_center.gz \
			-o d100.$pre.5k.${f}_center.plotHeatmap.pdf \
			--colorMap Greens --whatToShow 'heatmap and colorbar' \
			--plotFileFormat pdf --samplesLabel $labs \
			--heatmapWidth 10 \
			--kmeans 3 \
			--outFileSortedRegions d100.$pre.5k.${f}_center.k3.sortedRegions.bed &
	done	
done

for input in WGS.mut/d100/bedtools.merge.WGS.CO*
do
	pre=`echo $input | awk -F"/" '{print $NF}' | awk -F"." '{print $NF}'`
	awk -vOFS="\t" '{print $0, $3-$2}' $input | sort -k4,4nr | head -100 | awk -vOFS="\t" '{print $1,$2,$3,$1":"$2"-"$3"#"$4}' > d100.top100.$pre.bed
	for f in `ls bigwig.files/ | awk -F"." '{print $3}' | sort | uniq`
	do
		labs=`ls bigwig.files/*${f}*.bigwig | awk -F"/" '{print $NF}' | awk -F"." '{print $2"."$3}' | tr "\n" " "`
		computeMatrix reference-point \
			-R d100.top100.$pre.bed \
			-S bigwig.files/*${f}*.bigwig \
			-p 48 \
			-b 3000 -a 3000 -bs 5 \
			--referencePoint center \
			--skipZeros -o d100.top100.$pre.${f}_center.gz \
			--outFileNameMatrix d100.top100.$pre.${f}_center.tab \
			--outFileSortedRegions d100.top100.$pre.${f}_center.bed
		plotHeatmap -m d100.top100.$pre.${f}_center.gz \
			-o d100.top100.$pre.${f}_center.plotHeatmap.pdf \
			--colorMap Greens --whatToShow 'heatmap and colorbar' \
			--plotFileFormat pdf --samplesLabel $labs \
			--heatmapWidth 10 \
			--kmeans 3 \
			--outFileSortedRegions d100.top100.$pre.${f}_center.k3.sortedRegions.bed &
		computeMatrix scale-regions \
			-R d100.top100.$pre.bed \
			-S bigwig.files/*${f}*.bigwig \
			-p 48 \
			-b 1000 -a 1000 -bs 5 \
			--regionBodyLength 1000 \
			--skipZeros -o d100.top100.$pre.${f}_scaled.gz \
			--outFileNameMatrix d100.top100.$pre.${f}_scaled.tab \
			--outFileSortedRegions d100.top100.$pre.${f}_genes.bed
		plotHeatmap -m d100.top100.$pre.${f}_scaled.gz \
			-o d100.top100.$pre.${f}_scaled.plotHeatmap.pdf \
			--colorMap Greens --whatToShow 'heatmap and colorbar' \
			--plotFileFormat pdf --samplesLabel $labs \
			--heatmapWidth 10 \
			--kmeans 3 \
			--outFileSortedRegions d100.top100.$pre.${f}_scaled.k3.sortedRegions.bed &
	done	
done



for f in bedtools.merge.WGS.*
do
	awk '$3-$2<100' $f > single.$f
	awk '$3-$2>=100' $f > multip.$f
done

# WGS.mut/single.bedtools.merge.WGS.COAD-US  WGS.mut/single.bedtools.merge.WGS.COCA-CN
# WGS.mut/multip.bedtools.merge.WGS.COAD-US  WGS.mut/multip.bedtools.merge.WGS.COCA-CN
# for reg in WGS.mut/{single,multip}.bedtools.merge.WGS.CO*
for reg in WGS.mut/multip.bedtools.merge.WGS.CO*
do
	pre=`echo $reg | awk -F"/" '{print $NF}' | awk -F"." '{print $NF}'`
	for f in `ls bigwig.files/ | awk -F"." '{print $3}' | sort | uniq`
	do
		computeMatrix scale-regions \
			-R $reg \
			-S bigwig.files/*${f}*.bigwig \
			-p 48 \
			-b 1000 -a 1000 -bs 10 \
			--regionBodyLength 1000 \
			--skipZeros -o multip.$pre.${f}_scaled.gz \
			--outFileNameMatrix multip.$pre.${f}_scaled.tab \
			--outFileSortedRegions multip.$pre.${f}_genes.bed
		labs=`ls bigwig.files/*${f}*.bigwig | awk -F"/" '{print $NF}' | awk -F"." '{print $2"."$3}' | tr "\n" " "`
		plotHeatmap -m multip.$pre.${f}_scaled.gz \
			-o multip.$pre.${f}_scaled.plotHeatmap.pdf \
			--colorMap Greens --whatToShow 'heatmap and colorbar' \
			--plotFileFormat pdf --samplesLabel $labs \
			--heatmapWidth 10 \
			--kmeans 3 \
			--outFileSortedRegions multip.$pre.${f}.k3.sortedRegions.bed &
	done	
done
for reg in WGS.mut/single.bedtools.merge.WGS.CO*
do
	pre=`echo $reg | awk -F"/" '{print $NF}' | awk -F"." '{print $NF}'`
	for f in `ls bigwig.files/ | awk -F"." '{print $3}' | sort | uniq`
	do
		computeMatrix reference-point \
			-R $reg \
			-S bigwig.files/*${f}*.bigwig \
			-p 48 \
			-b 1000 -a 1000 -bs 10 \
			--referencePoint center \
			--skipZeros -o single.$pre.${f}_center.gz \
			--outFileNameMatrix single.$pre.${f}_center.tab \
			--outFileSortedRegions single.$pre.${f}_genes.bed
		labs=`ls bigwig.files/*${f}*.bigwig | awk -F"/" '{print $NF}' | awk -F"." '{print $2"."$3}' | tr "\n" " "`
		plotHeatmap -m single.$pre.${f}_center.gz \
			-o single.$pre.${f}center.plotHeatmap.pdf \
			--colorMap Greens --whatToShow 'heatmap and colorbar' \
			--plotFileFormat pdf --samplesLabel $labs \
			--heatmapWidth 10 \
			--kmeans 3 \
			--outFileSortedRegions single.$pre.${f}.k3.sortedRegions.bed &
	done	
done
wait




