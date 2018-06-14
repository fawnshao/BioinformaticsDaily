#!/bin/sh
for f in WGS.COAD-US.bed WGS.COCA-CN.bed
do
	pre=`echo $f | sed 's/.bed//'`
	cut -f 4 $f | cut -f 1 -d"|" | sort | uniq -c | awk -vOFS="\t" '{print $2"\t"$1}' | sort -k2,2nr > samplecount.$pre
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
