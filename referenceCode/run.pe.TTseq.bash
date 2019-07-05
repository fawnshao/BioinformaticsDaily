#!/bin/bash
#----------------------------------------------------
# Sample SLURM job script
#   for TACC Stampede2 KNL nodes
#
#   *** Serial Job on Normal Queue ***
#
# Last revised: 27 Jun 2017
#
# Notes:
#
#   -- Copy/edit this script as desired.  Launch by executing
#      "sbatch knl.serial.slurm" on a Stampede2 login node.
#
#   -- Serial codes run on a single node (upper case N = 1).
#        A serial code ignores the value of lower case n,
#        but slurm needs a plausible value to schedule the job.
#
#   -- For a good way to run multiple serial executables at the
#        same time, execute "module load launcher" followed
#        by "module help launcher".

#----------------------------------------------------

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.out       # Name of stdout output file
#SBATCH -e myjob.err       # Name of stderr error file
#SBATCH -p normal         # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 1             # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 10:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all   # Send email at begin and end of job

### UMI
# HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
hg19indexdir=/home1/04935/shaojf/scratch/star_index/hg19.star

WorkDir=/home1/04935/shaojf/scratch/Xiaoyu.Seq/NS28.TTseq
FastqDir=$WorkDir/FastqDir
CleanDir=$WorkDir/CleanDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
hubname="NS28-TTSeq"
runname="NS28-TTSeq"
Threads=272

cd $WorkDir
mkdir -p $BamDir
mkdir -p $BwDir
# mkdir -p $CleanDir

# # fastqc --threads $Threads --outdir $FastqDir $FastqDir/*.fastq.gz &
# # module load python3/3.6.3
# ### --local 
# for f in $FastqDir/*_R1_001.fastq.gz
# do
# 	echo $f
# 	i=`echo $f | awk -F"/" '{print $NF}' | sed 's?_R1_001.fastq.gz??'`
	
# 	# ~/.local/bin/cutadapt --nextseq-trim=20 -m 18 -j $Threads \
# 	# 	-o $CleanDir/$i.R1.clean.fastq -p $CleanDir/$i.R2.clean.fastq \
# 	# 	$FastqDir/${i}_R1_001.fastq.gz $FastqDir/${i}_R2_001.fastq.gz
# 	# --extract-method=string 
# 	umi_tools extract -I $FastqDir/${i}_R2_001.fastq.gz \
# 		--bc-pattern=NNNNNNNNNN --read2-in=$FastqDir/${i}_R1_001.fastq.gz \
# 		--stdout=$CleanDir/$i.R2.umiprocessed.fastq.gz \
# 		--read2-out=$CleanDir/$i.R1.umiprocessed.fastq.gz \
# 		-L $CleanDir/$i.extract.log &
# done
# wait

# --readFilesIn $CleanDir/$i.R1.umiprocessed.fastq.gz \
# 		$CleanDir/$i.R2.umiprocessed.fastq.gz \
for f in $FastqDir/*_R1_001.fastq.gz
do
	echo $f
	pre=`echo $f | awk -F"/" '{print $NF}' | awk -F"_" '{print $1}'`
	i=`awk -v var=$pre '$1==var{print $2}' files.txt | sed 's?+?_?'`
	STAR --runMode alignReads --runThreadN 100 \
		--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
		--readFilesIn $FastqDir/${i}_R1_001.fastq.gz \
		$FastqDir/${i}_R2_001.fastq.gz \
		--outSAMunmapped Within \
		--outFilterMultimapNmax 1 \
		--outFilterMultimapScoreRange 1 \
		--outFilterScoreMinOverLread  0.3 \
		--outFilterMatchNminOverLread 0.3 \
		--readFilesCommand zcat \
		--outFileNamePrefix $BamDir/$i. \
		--outSAMattributes All --outSAMtype BAM Unsorted \
		--outFilterType BySJout --outReadsUnmapped Fastx \
		--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
		--alignEndsType Local

	# bowtie2 -5 3 -3 1 -p $Threads -x $HG19 \
	# 	-1 $CleanDir/$i.R1.umiprocessed.fastq.gz \
	# 	-2 $CleanDir/$i.R2.umiprocessed.fastq.gz -S $BamDir/$i.bt2.sam
	samtools view -q 10 -@ 272 -bo $BamDir/$i.q10.bam $BamDir/$i.Aligned.out.bam
	makeTagDirectory $i.mTD -fragLength pe -sspe -flip $BamDir/$i.q10.bam &
	# samtools sort -@ 68 -o $BamDir/$i.q10.srt.bam $BamDir/$i.q10.bam
	# samtools index $BamDir/$i.q10.srt.bam
	# umi_tools dedup -I $BamDir/$i.q10.srt.bam -S $BamDir/$i.q10.dedup.bam -L $BamDir/$i.umi_tools.log \
	# 	--paired --output-stats $i.umi_tools.stats
	# makeTagDirectory $i.mTD -sspe -flip $BamDir/$i.q10.dedup.bam &
done
wait
makeMultiWigHub.pl $hubname hg19 -url $BwDir -webdir $BwDir -d *.mTD -force -strand &

analyzeRepeats.pl rna hg19 -rpkm -strand + -condenseGenes -d *.mTD > $runname.gene.rpkm.txt &
# awk -vOFS="\t" '{print $4,$1,$2,$3,$6,$5}' MCF-7.p300.55905_countedputative.eRNAs.siCTL.bed > MCF-7.p300.55905_countedputative.eRNAs.txt
analyzeRepeats.pl HCT116.eRNA.sjf.notinblacklist.forHomer.merged.txt hg19 -rpkm -strand + -d *.mTD > $runname.eRNA.rpkm.txt &
wait