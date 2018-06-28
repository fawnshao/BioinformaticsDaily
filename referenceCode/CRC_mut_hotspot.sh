#!/bin/sh
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
