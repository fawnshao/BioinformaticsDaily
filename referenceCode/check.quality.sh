#!/bin/bash
for f in *.fastq.gz
do
	gunzip -c $f | grep -c "^+"
done

for f in *.Aligned.out.bam
do
	pre=`echo $f | sed 's?.Aligned.out.bam??'`
	samtools view $f | awk -v var=$pre '$9>0{print var"\t"$9}' >> STAR.mapstats.tsv
done
