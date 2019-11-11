#!/bin/bash
workdir=/home1/04935/shaojf/scratch/Xiaoyu.Seq/NovoGene.2-3
fastqDir=$workdir/FastqDir
fastqcDir=$workdir/FastQCDir
mkdir $fastqcDir

fastqc --threads 272 --outdir $fastqcDir $fastqDir/*.gz

cd $fastqcDir
for f in *.zip
do
	unzip $f
done

grep "Total Sequences" *_1_fastqc/fastqc_data.txt | awk '{print $1"\t"$3}'
grep -A 2 "Overrepresented sequences" *_1_fastqc/fastqc_data.txt | sed -n '3~4p' | awk '{print $1"\t"$3}'