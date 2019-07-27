#!/bin/bash
#SBATCH -J myjob           # Job name
#SBATCH -o myjob.out       # Name of stdout output file
#SBATCH -e myjob.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 20:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

# Sequenced reads were trimmed for adaptor sequences using skewer 0.1.127 with parameters -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -t 16 -q 21 -l 21 -n -u -f sanger for bulk-population RNA-seq and parameters -x CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG -y CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT -t 16 -q 21 -l 21 -n -u -f sanger for single-cell RNA-seq. 
# Trimmed reads were then mapped to hg38 using STAR 2.4 with parameters  --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD  --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20 --alignIntronMax 1000000   --alignMatesGapMax 1000000  --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --twopassMode Basic --twopass1readsN -1 --outWigStrand Stranded
# Normalized expression values (TPM) and transcript counts were then calculated using RSEM 1.2.21 with parameters --bam --estimate-rspd --seed 12345 --paired-end --forward-prob 0
# Genome_build: hg38

# The --shift -100 --extsize 200 option centers a 200 bp window on the Tn5 binding site, which is more accurate for ATAC-seq data (the same option is also applicable for DNase-seq data). The 5' ends of reads represent the Tn5 (and DNase I) cut sites; so the 5' ends of reads represent the most accessible regions. If you have paired-end ATAC-seq reads, many of your read pairs will span at least one nucleosome. If you look at plots of distribution of fragment length, you should see many fragments that don't span nucleosomes (less 150-200 bp), but you will also see many fragments that do span nucleosomes (>200-300 bp). Therefore, you don't want to run MACS2 with the default model building parameters meant for ChIP-seq. These model building parameters assume that the binding site is in the middle of the fragment, which is accurate for ChIP-seq, but not for ATAC-seq. If you use the ChIP-seq options for ATAC-seq, many of your peaks will be centered on nucleosomes instead of centered on the accessible regions. The reason why a 200 bp window is often used (--shift 100 --extsize 200), is that it is reasonable to assume that ATAC-seq peaks are at least 200 bp long because this is about the size of a nucleosome-free region with a single nucleosome removed. Some people may use --shift -75 --extsize 150, with the assumption that the length of an accessible region with a single nucleosome removed is about 150 bp, which is also reasonable.

# ref to https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/
# and https://www.encodeproject.org/documents/c008d7bd-5d60-4a23-a833-67c5dfab006a/@@download/attachment/ATACSeqPipeline.pdf


#### additional analysis for ATAC-Seq
WorkDir=/home1/04935/shaojf/scratch/T21/NS26.ATAC1
FastqDir=$WorkDir/FastqDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
hubname="NS26.T21.ATAC"
runname="NS26.T21.ATAC"
Threads=272

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
blackListFileName=hg19.ENCFF001TDO.bed

mkdir -p $BamDir
mkdir -p $BwDir
Threads=272

for f in $FastqDir/UTA1-Xiaoyu-Genea21-NSC-P8-500k_S27_R1_001.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's?_R1_001.fastq.gz??'`
	bowtie2 --local --very-sensitive -p $Threads -x $HG19 \
		-1 $FastqDir/${i}_R1_001.fastq.gz \
		-2 $FastqDir/${i}_R2_001.fastq.gz -S $BamDir/$i.bt2.sam
	
	samtools view -1 -q 10 -bo $BamDir/$i.q10.bam --threads $Threads $BamDir/$i.bt2.sam
	samtools sort --threads 20 $BamDir/$i.q10.bam -o $BamDir/$i.sorted.bam
	samtools fixmate --threads $Threads -m $BamDir/$i.q10.bam $BamDir/$i.fixmate.bam
	samtools sort --threads 20 $BamDir/$i.fixmate.bam -o $BamDir/$i.fixmate.srt.bam
	samtools markdup -r -s -O BAM -@ $Threads $BamDir/$i.fixmate.srt.bam $BamDir/$i.rmdup.bam
	samtools sort -n --threads 20 $i.rmdup.bam -o $i.rmdup.srtbyn.bam

	samtools view -@ 50 -h $i.rmdup.bam | \
		grep -v "chrM" | samtools view -@ 50 -bo nochrM.$i.rmdup.bam -
	samtools index -@ 50 nochrM.$i.rmdup.bam

	bedtools bamtobed -bedpe -i $i.rmdup.srtbyn.bam > $i.rmdup.srtbyn.bedpe
	awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}if($9=="-"){print $1,$2-5,$6-5}}' \
		$i.rmdup.srtbyn.bedpe > $i.rmdup.srtbyn.shift.bedpe
	macs2 callpeak --verbose 3 -t $i.rmdup.srtbyn.bedpe -f BEDPE -g hs \
		-n $i.bedpe --keep-dup all --cutoff-analysis 1> $i.bedpe.log 2>&1 &
	macs2 callpeak --verbose 3 -t $i.rmdup.srtbyn.shift.bedpe -f BEDPE -g hs \
		-n $i.shift.bedpe --keep-dup all --cutoff-analysis 1> $i.shift.bedpe.log 2>&1 &
	macs2 callpeak --verbose 3 -t $i.rmdup.bam -f BAMPE -g hs \
		-n $i.bampe --keep-dup all --cutoff-analysis 1> $i.bampe.log 2>&1 &
	
	makeTagDirectory $i.mTD -fragLength 200 -genome hg19 -checkGC -tbp 1 nochrM.$i.rmdup.bam &
	# makeTagDirectory $i.mTD -fragLength pe -genome hg19 -checkGC -tbp 1 nochrM.$i.rmdup.bam &
	bamCoverage -b nochrM.$i.rmdup.bam -o $i.bigWig -of bigwig -bs 1 -p 50 \
        --extendReads --normalizeUsing CPM --effectiveGenomeSize 2864785220 \
        --minFragmentLength 150 --maxFragmentLength 1000 &
    bamCoverage -b nochrM.$i.rmdup.bam -o nolimit.$i.bigWig -of bigwig -bs 1 -p 50 \
        --extendReads --normalizeUsing CPM --effectiveGenomeSize 2864785220 &
    bamCoverage -b nochrM.$i.rmdup.bam -o RPGC.$i.bigWig -of bigwig -bs 1 -p 50 \
        --blackListFileName $blackListFileName --extendReads \
        --normalizeUsing RPGC --effectiveGenomeSize 2864785220 --ignoreForNormalization chr21 &
done
wait
makeMultiWigHub.pl $runname hg19 -url $BwDir -webdir $BwDir -d *.mTD


for f in *.bam
do
	i=`echo $f | sed 's?.bam??'`
	samtools sort -n --threads 20 $f -o $i.srtbyn.bam
	bedtools bamtobed -bedpe -i $i.srtbyn.bam > $i.srtbyn.bedpe
	# awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}if($9=="-"){print $1,$2-5,$6-5}}' \
	# 	$i.srtbyn.bedpe > $i.srtbyn.shift.bedpe
	# need to filter the unproper pairs
	macs2 callpeak --verbose 3 -t $i.srtbyn.bedpe -f BEDPE -g hs \
		-n $i.bedpe --keep-dup all --cutoff-analysis 1> $i.bedpe.log 2>&1 &

	alignmentSieve --numberOfProcessors 270 --ATACshift --bam $f -o $i.shift.bam
	macs2 callpeak --verbose 3 -t $i.shift.bam -f BAMPE -g hs \
		-n $i.shift.bampe --keep-dup all --cutoff-analysis 1> $i.shift.bampe.log 2>&1 &
done
