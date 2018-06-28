#!/bin/bash
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3319272/
# In order to characterize these targets globally, we developed two independent genome-wide R-loop profiling methods. The first one, which we termed DRIP (DNA:RNA ImmunoPrecipitation) relies on the intrinsic specificity of the S9.6 antibody for R-loop molecules. The second method, which we termed DRIVE (DNA:RNA in vitro Enrichment), makes use of a catalytically-deficient, but binding-competent, human RNASEH1 mutant protein in affinity pulldown assays.
# From 4.8 to 6 million mapped reads were obtained for each sequenced library. Alignment to the Hg19 build was carried out using BWA (Li and Durbin, 2009), and peak calling was done using MACS (Zhang et al., 2008). 
# For DRIP-seq, peaks were called using all mapped reads enforcing a greater than 10-fold enrichment above both input and RNase H-pre-treated control datasets. For DRIVE-seq, peaks were called by using uniquely mapped reads enforcing a greater than 5-fold enrichment above both input and RNase H-pre-treated samples, as well as an FDR of <0.1. 

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5030092/
# Raw reads from the DRIP experiments were aligned to the reference genome hg19/GRCh37 using bowtie2. An interval file of restriction fragments for the enzyme cocktail used in DRIP was obtained by verbose search of restriction site sequences with bowtie. Aligned reads with Q score over 10 were counted over the intervals separated by restriction enzyme cut sites using bedtools. Intervals with lower coverage (less than 10 counts per interval in any sample) were removed from the data set. Regions enriched over background were discovered with DESeq2. To perform differential analysis using triplicate data, we combined the set of regions with positive DRIP signal in at least one condition, together with the matched set of negative intervals and performed differential analysis with DESeq2. Read counts were normalized to the total number of mapped reads. To perform analysis of functional term enrichment, we ran GREAT analysis (McLean et al., 2010) on sets of differential DRIP regions.

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

name=GSE93368.Hela-DRIP-Seq
workDir=/home1/04935/shaojf/scratch/TIP60.project/GSE93368.Hela-DRIP-Seq
fastqDir=$workDir/fastqDir
mapDIR=$workDir/map.res
macsDIR=$workDir/macs2.res
visdir=$workDir/bigwigs

mkdir -p $mapDIR
mkdir -p $macsDIR
mkdir -p $visdir
cd $workDir

# prefetch SRR5164421 SRR5164422 SRR5683211 SRR5683212
# fastq-dump --split-files *.sra
# cut -f 6,9,28 RunInfo.txt > simRunInfo.txt

for f in $fastqDir/*_1.fastq.gz
do
	pre=`echo $f | awk -F"/" '{print $NF}' | sed 's/_1.fastq.gz//'`
	bowtie2 --local --very-sensitive -p 68 -x $HG19 -1 $fastqDir/${pre}_1.fastq.gz -2 $fastqDir/${pre}_2.fastq.gz -S $mapDIR/$pre.sam
	samtools view -1 -q 10 --threads 68 $mapDIR/$pre.sam | samtools sort --threads 68 > $mapDIR/$pre.sorted.bam
	makeTagDirectory $pre.mTD -genome hg19 -checkGC -tbp 1 $mapDIR/$pre.sorted.bam &
	samtools sort --threads 68 -n -o $mapDIR/$pre.namesorted.bam $mapDIR/$pre.sorted.bam
	samtools fixmate --threads 68 -m $mapDIR/$pre.namesorted.bam $mapDIR/$pre.fixmate.bam
	samtools sort --threads 68 -o $mapDIR/$pre.possorted.bam $mapDIR/$pre.fixmate.bam
	samtools markdup -r --threads 68 $mapDIR/$pre.possorted.bam $mapDIR/$pre.rmdup.bam &
done
wait
makeMultiWigHub.pl $name hg19 -url $visdir -webdir $visdir -d *.mTD &

for pair in SRR5164421-SRR5164422-S9.6DRIPSeq.1 SRR5683211-SRR5683212-S9.6DRIPSeq.2
do
	treat=`echo $pair | awk -F"-" '{print $1}'`
	input=`echo $pair | awk -F"-" '{print $2}'`
	pre=`echo $pair | awk -F"-" '{print $3}'`
	macs2 callpeak -t $mapDIR/$treat.rmdup.bam -c $mapDIR/$input.rmdup.bam -g hs -n $pre --outdir $macsDIR &
done
wait

samtools view -1 -q 10 $pre.sam | samtools sort > $pre.sorted.bam
samtools sort -n -o $pre.namesorted.bam $pre.sorted.bam
samtools fixmate -m $pre.namesorted.bam $pre.fixmate.bam
samtools sort -o $pre.possorted.bam $pre.fixmate.bam
samtools markdup -r $pre.possorted.bam $pre.rmdup.bam
