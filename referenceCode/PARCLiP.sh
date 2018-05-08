#!/bin/sh

#SBATCH -J PARCLIP           # Job name
#SBATCH -o PARCLIP.out       # Name of stdout output file
#SBATCH -e PARCLIP.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 68              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 29:50:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

### for STAR
hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
hg19indexdir=/home1/04935/shaojf/scratch/star_index/hg19.star
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

# cut -f 1,9 SRP064870.runinfo.tsv | grep OTHER | cut -f 2 | xargs -n 1 prefetch
cd /home1/04935/shaojf/scratch/HEK293T_PARCLiP

date
while read line
do
	acc=`echo $line | awk '{print $1}'`
	pre=`echo $line | awk '{print $2}'`
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn ${acc}.fastq \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix ${pre}. \
	--outSAMattributes All --outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local

	samtools sort -@ 68 ${pre}.Aligned.out.srt.bam ${pre}.Aligned.out.bam
	samtools index -@ 68 ${pre}.Aligned.out.srt.bam
	make_bigwig_files.py --bam ${pre}.Aligned.out.srt.bam --genome $genomesize --bw_pos ${pre}.norm.pos.bw --bw_neg ${pre}.norm.neg.bw
	bedSort ${pre}.Aligned.out.srt.norm.pos.bg ${pre}.Aligned.out.srt.norm.pos.srt.bg &
	bedSort ${pre}.Aligned.out.srt.norm.neg.t.bg ${pre}.Aligned.out.srt.norm.neg.t.srt.bg &
	wait
	bedGraphToBigWig ${pre}.Aligned.out.srt.norm.pos.srt.bg $genomesize bedGraphToBigWig.${pre}.norm.pos.bw &
	bedGraphToBigWig ${pre}.Aligned.out.srt.norm.neg.t.srt.bg $genomesize bedGraphToBigWig.${pre}.norm.neg.bw &
done < sim.info.txt
wait
