#!/bin/bash
# for f in `ls ChIP/`
for f in `ls $1/`
do
	name=`echo $f | sed 's/.ucsc.bigWig//'`
	echo "track $name"
	echo "bigDataUrl $1/$f"
	echo "shortLabel $name"
	echo "longLabel $name"
	echo "type bigWig"
	a=`seq 0 255 | shuf | head -1`
	b=`seq 0 255 | shuf | head -1`
	c=`seq 0 255 | shuf | head -1`
	echo "color $a,$b,$c"
	echo "autoScale on"
	echo "visibility full"
	echo "maxHeightPixels 50:30:8"
	echo ""
done
