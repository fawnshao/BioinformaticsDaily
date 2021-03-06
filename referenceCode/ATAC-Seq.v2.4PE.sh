#!/bin/bash
#SBATCH -J myjob           # Job name
#SBATCH -o myjob.out       # Name of stdout output file
#SBATCH -e myjob.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 8:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

# Sequenced reads were trimmed for adaptor sequences using skewer 0.1.127 with parameters -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -t 16 -q 21 -l 21 -n -u -f sanger for bulk-population RNA-seq and parameters -x CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG -y CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT -t 16 -q 21 -l 21 -n -u -f sanger for single-cell RNA-seq.
# Trimmed reads were then mapped to hg38 using STAR 2.4 with parameters  --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD  --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000  --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --twopassMode Basic --twopass1readsN -1 --outWigStrand Stranded
# Normalized expression values (TPM) and transcript counts were then calculated using RSEM 1.2.21 with parameters --bam --estimate-rspd --seed 12345 --paired-end --forward-prob 0
# Genome_build: hg38

# The --shift -100 --extsize 200 option centers a 200 bp window on the Tn5 binding site, which is more accurate for ATAC-seq data (the same option is also applicable for DNase-seq data). The 5' ends of reads represent the Tn5 (and DNase I) cut sites; so the 5' ends of reads represent the most accessible regions. If you have paired-end ATAC-seq reads, many of your read pairs will span at least one nucleosome. If you look at plots of distribution of fragment length, you should see many fragments that don't span nucleosomes (less 150-200 bp), but you will also see many fragments that do span nucleosomes (>200-300 bp). Therefore, you don't want to run MACS2 with the default model building parameters meant for ChIP-seq. These model building parameters assume that the binding site is in the middle of the fragment, which is accurate for ChIP-seq, but not for ATAC-seq. If you use the ChIP-seq options for ATAC-seq, many of your peaks will be centered on nucleosomes instead of centered on the accessible regions. The reason why a 200 bp window is often used (--shift 100 --extsize 200), is that it is reasonable to assume that ATAC-seq peaks are at least 200 bp long because this is about the size of a nucleosome-free region with a single nucleosome removed. Some people may use --shift -75 --extsize 150, with the assumption that the length of an accessible region with a single nucleosome removed is about 150 bp, which is also reasonable.

#### additional analysis for ATAC-Seq
runname="ATAC_test"
workdir=/home1/04935/shaojf/scratch/ATAC_test
FastqDir=$workdir/FastqsDir
FastqcDir=$workdir/FastqcDir
TrimadapterDir=$workdir/TrimmedFastqsDir
VisDir=$workdir/bigwigs
Threads=200

HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

mkdir $FastqcDir
mkdir $TrimadapterDir
mkdir $VisDir

fastqc --threads $Threads --outdir $FastqcDir $FastqDir/*.fastq.gz &

module load python3/3.6.3
for f in $FastqDir/*_1.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/_1.fastq.gz//'`

	~/.local/bin/cutadapt -j $Threads --nextseq-trim=20 \
	-a tag1=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
	-a tag1rev=CTGTCTCTTATACACATCTGACGCTGCCGACGA \
	-a index1=GTCTCGTGGGCTCGG \
	-a index1rev=CCGAGCCCACGAGAC \
	-a tag2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
	-a tag2rev=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
	-a index2=TCGTCGGCAGCGTC \
	-a index2rev=GACGCTGCCGACGA \
	-A tag1=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
	-A tag1rev=CTGTCTCTTATACACATCTGACGCTGCCGACGA \
	-A index1=GTCTCGTGGGCTCGG \
	-A index1rev=CCGAGCCCACGAGAC \
	-A tag2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
	-A tag2rev=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
	-A index2=TCGTCGGCAGCGTC \
	-A index2rev=GACGCTGCCGACGA \
	-m 18 \
	-o $TrimadapterDir/$i.clean_1.fastq.gz -p $TrimadapterDir/$i.clean_2.fastq.gz \
	$FastqDir/${i}_1.fastq.gz $FastqDir/${i}_2.fastq.gz &
done
wait
###### my pipeline
for f in $FastqDir/*_1.fastq.gz
do
	echo $f
	i=`echo $f | awk -F"/" '{print $NF}' | sed 's/_1.fastq.gz//'`
	echo $i

	bowtie2 --local --very-sensitive -p $Threads -x $HG19 \
	-1 $TrimadapterDir/$i.clean_1.fastq.gz -2 $TrimadapterDir/$i.clean_2.fastq.gz -S $i.sam
	samtools view -1 -q 10 -bo $i.q10.bam --threads $Threads $i.sam
	samtools sort --threads 20 $i.q10.bam -o $i.sorted.bam
	makeTagDirectory $i.mTD -genome hg19 -checkGC -tbp 1 -format sam <(samtools view -@ 10 -h $i.sorted.bam | grep -v "chrM")
	findPeaks $i.mTD -style factor -o auto &
	samtools fixmate --threads $Threads -m $i.q10.bam $i.fixmate.bam
	samtools sort --threads 20 $i.fixmate.bam -o $i.fixmate.srt.bam
	samtools markdup -r -s -O BAM -@ $Threads $i.fixmate.srt.bam $i.rmdup.bam
	macs2 callpeak --verbose 3 -t $i.rmdup.bam -f BAM -g hs -n $i --keep-dup all --nomodel --shift -100 --extsize 200 --broad 1> $i.log 2>&1 &
	bedtools bamtobed -bedpe -i $i.rmdup.bam > $i.rmdup.bedpe
	macs2 callpeak --verbose 3 -t $i.rmdup.bedpe -f BEDPE -g hs -n $i.bedpe --keep-dup all --broad 1> $i.bedpe.log 2>&1 &
done
makeMultiWigHub.pl $runname hg19 -url $VisDir -webdir $VisDir -d *.mTD
wait
