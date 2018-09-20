#!/bin/bash
#SBATCH -J seqvis           # Job name
#SBATCH -o seqvis.out       # Name of stdout output file
#SBATCH -e seqvis.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 10:50:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

### for quick check and visualization
### for STAR
hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
hg19indexdir=/home1/04935/shaojf/scratch/star_index/hg19.star
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

WorkDir=/home1/04935/shaojf/scratch/RBMX_CLIP.GSE115429/
FastqDir=$WorkDir/FastqDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
hubname=RBMX_CLIP
# fastq-dump --split-files *.sra &
cd $WorkDir
mkdir -p $BamDir
mkdir -p $BwDir

# SRR7274971	RBMX_CLIP_rep1
# SRR7274972	RBMX_CLIP_rep2
# SRR7274973	RBMX_SMInput_rep1
# SRR7274974	RBMX_SMInput_rep2
# SRR7274975	RBMX_SEC_fraction6
# SRR7274976	RBMX_SEC_fraction20
date
while read line
do
	acc=`echo $line | awk '{print $1}'`
	pre=`echo $line | awk '{print $2}'`
	STAR --runMode alignReads --runThreadN 68 \
		--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
		--readFilesIn $FastqDir/${acc}_1.fastq $FastqDir/${acc}_2.fastq \
		--outSAMunmapped Within \
		--outFilterMultimapNmax 1 \
		--outFilterMultimapScoreRange 1 \
		--outFileNamePrefix $BamDir/${pre}. \
		--outSAMattributes All --outSAMtype BAM Unsorted \
		--outFilterType BySJout --outReadsUnmapped Fastx \
		--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
		--alignEndsType Local

	samtools view -@ 272 -q 10 -bo $BamDir/${pre}.q10.bam $BamDir/${pre}.Aligned.out.bam
	samtools sort -@ 12 -o $BamDir/${pre}.q10.srt.bam $BamDir/${pre}.q10.bam
	samtools index -@ 68 $BamDir/${pre}.q10.srt.bam
	# samtools markdup -r -@ 272 $BamDir/${pre}.q10.srt.bam $BamDir/${pre}.rmdup.bam
	makeTagDirectory $pre.mTD $BamDir/${pre}.q10.srt.bam -sspe -tbp 1 &
done < sim.info.txt
wait
makeMultiWigHub.pl $hubname hg19 -strand -url $BwDir -webdir $BwDir -d *.mTD
