#!/bin/bash
WorkDir=/data/shaojf/LiLab.Seq/NS36
FastqDir=$WorkDir/FastqDir
FastqcDir=$WorkDir/FastqcDir
mkdir $FastqcDir
fastqc --threads 60 --outdir $FastqcDir $FastqDir/*.gz

cd $FastqcDir
for f in *.zip
do
	unzip $f &
done
wait

# grep "Total Sequences" *_1_fastqc/fastqc_data.txt | awk '{print $1"\t"$3}'
# grep -A 2 "Overrepresented sequences" *_1_fastqc/fastqc_data.txt | sed -n '3~4p' | awk '{print $1"\t"$3}'

grep "Total Sequences" *_R1_001_fastqc/fastqc_data.txt | awk '{print $1"\t"$3}'
grep -A 2 "Overrepresented sequences" *_R1_001_fastqc/fastqc_data.txt | sed -n '3~4p' | awk '{print $1"\t"$3}'

# grep "Total Sequences" NS36-Moliu-ChIP**_R1_001_fastqc/fastqc_data.txt | awk '{print $1"\t"$3}'
# grep -A 2 "Overrepresented sequences" NS36-Moliu-ChIP*_R1_001_fastqc/fastqc_data.txt | sed -n '3~4p' | awk '{print $1"\t"$3}'