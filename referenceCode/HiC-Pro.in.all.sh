#!/bin/bash
THREADS=30
MAPQ_Thr=30
GSIZE='hs'
# default values of short and long read thresholds
ShortReadDistThr=1000
LongReadDistThr=10000

CUT_ENZ=/data/shaojf/myReference/hg19_MboI_resfrag.bed
GENOME=/data/shaojf/myReference/bwa-index/hg19.fa
genomesize=/data/shaojf/myReference/star_index/hg19.chrom.sizes

workdir=/data/shaojf/LiLab.Seq/MCF7.SSBP1.HiChIP
fqDir=$workdir/FastqDir
fastqcDir=$workdir/FastQCDir
hicpro2higlass=$workdir/hicpro2higlass
# stats=$workdir/stats
mkdir $fastqcDir
mkdir $hicpro2higlass
# mkdir $stats
fastqc --threads 68 --outdir $fastqcDir $fastqDir/*.fq.gz &

inputfiles=$workdir/ngsrun.txt
while read line
do
	acc=`echo $line | awk '{print $1}'`
	name=`echo $line | awk '{print $2}'`

	FASTQ1=$fqDir/${acc}_R1_001.fastq.gz
	FASTQ2=$fqDir/${acc}_R2_001.fastq.gz
	PREFIX=$name
	# HiC-Pro pipeline
	HiCProInputDir=$workdir/HiCPro.input/$acc
	HiCProOutputDir=$workdir/HiCPro.output/$acc
	mkdir -p $HiCProInputDir/$acc
	#mkdir -p $HiCProOutputDir/$acc
	ln -s $FASTQ1 $HiCProInputDir/$acc/${acc}_1.fastq.gz
	ln -s $FASTQ2 $HiCProInputDir/$acc/${acc}_2.fastq.gz
	~/myTools/HiC-Pro_2.11.1/bin/HiC-Pro -i $HiCProInputDir -o $HiCProOutputDir -c $workdir/config-hicpro.txt
done < $inputfiles
wait

cd $fastqcDir
for f in *.zip
do
	unzip $f
done
echo "Total Sequences"
grep "Total Sequences" *_R1_001_fastqc/fastqc_data.txt | sed 's?_R1_001_fastqc?\t?' | awk '{print $1"\t"$4}'
echo "adapter%"
grep -A 2 "Overrepresented sequences" *_R1_001_fastqc/fastqc_data.txt | \
	sed -n '3~4p' | sed 's?_R1_001_fastqc?\t?' | awk '{print $1"\t"$4}'

# cd $stats
# ln -s $workdir/HiCPro.output/*/hic_results/stats/* .
echo "valid interaction"
grep -w "valid_interaction" $workdir/HiCPro.output/*/hic_results/stats/*/*_allValidPairs.mergestat 
grep -w "valid_interaction_rmdup" $workdir/HiCPro.output/*/hic_results/stats/*/*_allValidPairs.mergestat 


cd $hicpro2higlass
ln -s $workdir/HiCPro.output/*/hic_results/data/*/*.allValidPairs .
for f in *.allValidPairs
do
    hicpro2higlass.sh -i $f -r 20000 -c /data/shaojf/myReference/hg19.chromInfo.sim.txt &
done
