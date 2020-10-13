#!/bin/bash
conda activate ngs
bismark_ref=/public/workspace/shaojf/myref/bismark.hg38/
WorkDir=/public/workspace/shaojf/wgbs.analysis
FastqDir=$WorkDir/FastqDir
prefix=$1

for fastq in $FastqDir/$prefix*.fq.bz2
do
	echo $fastq
	echo `date`
	name=`basename $fastq | sed 's/.fq.bz2//'`
	bismark --phred64-quals --local --bowtie2 -p 30 $bismark_ref $fastq
	bismark_methylation_extractor -s --gzip --parallel 30 --comprehensive --bedGraph \
		--genome_folder $bismark_ref ${fastq}_bismark_bt2.sam
	echo `date`
	echo "Finished!"
done

conda deactivate
# Perl module GD::Graph::lines is not installed, skipping drawing M-bias plots (only writing out M-bias plot table)
# z - C in CpG context - unmethylated
# Z - C in CpG context - methylated
# x - C in CHG context - unmethylated
# X - C in CHG context - methylated
# h - C in CHH context - unmethylated
# H - C in CHH context - methylated
# u - C in Unknown context (CN or CHN) - unmethylated
# U - C in Unknown context (CN or CHN) - methylated
# . - not a C or irrelevant position
# a=0;b=0;c=0;
# for(i=1;i<=length($0);i++){
# 	if(substr($0,i,1)=="Z"){a++;}
# 	else if(substr($0,i,1)=="X"){b++;}
# 	else if(substr($0,i,1)=="H"){c++;}
# }
# print a"\t"b"\t"c}
