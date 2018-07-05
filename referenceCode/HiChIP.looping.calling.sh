#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Stampede2 SKX nodes
#
#   *** MPI Job on SKX Normal Queue ***
# 
# Last revised: 20 Oct 2017
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch skx.mpi.slurm" on Stampede2 login node.
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do not use mpirun or mpiexec.
#
#   -- Max recommended MPI ranks per SKX node: 48
#      (start small, increase gradually).
#
#   -- If you're running out of memory, try running
#      fewer tasks per node to give each task more memory.
#
#----------------------------------------------------

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.out       # Name of stdout output file
#SBATCH -e myjob.err       # Name of stderr error file
#SBATCH -p long            # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 2               # Total # of mpi tasks
#SBATCH -t 120:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-user=dreambetter@gmail.com
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A ncvar_function  # Allocation name (req'd if you have more than 1)
# Advanced Usage:
# /usr/local/bin/ibrun -n <number of processors> -o <offset into processor hostlist> executable <execuable_opions>

# In this case you can specify a subset of processors from the
# list of all hosts on which to run the executable, and an offset into the list
# of processors to start numbering.  This allows you to e.g. launch two different
# exectuables on two different subsets of the available hosts.

# For example, the following  batch environment will allocate 
#  two nodes with 28 tasks per node: 
#    #$ -n 56
#    #$ -N 2 

#  We can run two independent jobs on each node:
#  /usr/local/bin/ibrun -n 28 -o 0  ./mpihello &
#  /usr/local/bin/ibrun -n 28 -o 28 ./mpihello &
#  wait

# The first call launches a 28-task job on the first 28 cores in the hostfile,
# The second call launches a 28-task job on the second 28 cores in the hostfile.
# The shell 'wait' command waits for all processes to finish before the shell exits.

# hichipper
# http://aryee.mgh.harvard.edu/hichipper/
# http://hichipper.readthedocs.io/en/latest/
# conda install -c bioconda hichipper 

# FitHiChIP
# https://github.com/ay-lab/FitHiChIP

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101498
# https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP112520
# A primary HCASMC line derived from a normal human donor heart

# prefetch --max-size 200G SRR5831489

# /home1/04935/shaojf/stampede2/refs/TADs/TAD_lists/Mumbach_NG2017_HiChIP_GSE101498
# GSM2705041_GM_HiChIP_H3K27ac_B1_allValidPairs.txt.gz
# GSM2705042_GM_HiChIP_H3K27ac_B2_allValidPairs.txt.gz
# cut -f 8,14 run.info | sed 's/ /_/g' > run.info.sim
# /home1/04935/shaojf/myTools/HiChIP/HiC-Pro/bin/utils/digest_genome.py -r mboi -o hg19.MboI.bed ~/scratch/bwa-index/hg19.fa
# Analyzing /home1/04935/shaojf/scratch/bwa-index/hg19.fa
# Restriction site(s) GATC
# /home1/04935/shaojf/myTools/HiChIP/HiC-Pro/annotation/hg19.MboI.bed

# /home1/04935/shaojf/scratch/HiChIP.test/map.res/GM12878_cell_line.SRR5831489/FitHiChIP_Peak2ALL_b5000/Interaction_ALL


export PATH=/home1/04935/shaojf/myTools/HiChIP/FitHiChIP/Preprocess:/home1/04935/shaojf/myTools/HiChIP/FitHiChIP/bin:/home1/04935/shaojf/myTools/HiChIP/FitHiChIP/scripts/BAM:/home1/04935/shaojf/myTools/HiChIP/FitHiChIP/scripts/PAIRIX:/home1/04935/shaojf/myTools/HiChIP/FitHiChIP:$PATH
THREADS=200
MAPQ_Thr=30
GSIZE='hs'
# default values of short and long read thresholds
ShortReadDistThr=1000
LongReadDistThr=10000

CUT_ENZ=/home1/04935/shaojf/myTools/HiChIP/hg19_MboI_resfrag.bed
GENOME=/home1/04935/shaojf/scratch/bwa-index/hg19.fa
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
configtemp=/home1/04935/shaojf/scratch/HiChIP.test/configfile_BAM

workdir=/home1/04935/shaojf/scratch/HiChIP.test
rawDir=$workdir/raw.data
fqDir=$workdir/fastq.data
mapDir=$workdir/map.res
inputfiles=$workdir/raw.data/a
mkdir -p $mapDir
mkdir -p $fqDir
ln -s /home1/04935/shaojf/myTools/HiChIP/hichipper/RestrictionFragmentFiles
alias sort="/bin/sort -S 50% --parallel=$THREADS"

while read line
do
	sra=`echo $line | awk '{print $1}'`
	name=`echo $line | awk '{print $2}'`
	# fastq-dump --gzip --split-files $rawDir/$sra.sra
	# mv ${sra}_?.fastq.gz $fqDir
	FASTQ1=$fqDir/${sra}_1.fastq.gz
	FASTQ2=$fqDir/${sra}_2.fastq.gz
	PREFIX=$name"."$sra

	# HiC-Pro pipeline
	HiCProInputDir=$workdir/HiCPro.input
	HiCProOutputDir=$workdir/HiCPro.output/$sra
	mkdir -p $HiCProInputDir/$sra
	# mkdir -p $HiCProOutputDir/$sra
	ln -s $FASTQ1 $HiCProInputDir/$sra/${sra}_1.fastq.gz
	ln -s $FASTQ2 $HiCProInputDir/$sra/${sra}_2.fastq.gz
	HiC-Pro -i $HiCProInputDir -o $HiCProOutputDir -c $workdir/config-hicpro.txt
	FASTQFILE=$HiCProOutputDir/inputfiles_.txt; export FASTQFILE
	make --file /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/home1/04935/shaojf/scratch/HiChIP.test/config-hicpro.txt CONFIG_SYS=/home1/04935/shaojf/bin/HiC-Pro_2.10.0/config-system.txt all_sub 2>&1
	cd $HiCProOutputDir
	make --file /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/home1/04935/shaojf/scratch/HiChIP.test/config-hicpro.txt CONFIG_SYS=/home1/04935/shaojf/bin/HiC-Pro_2.10.0/config-system.txt all_persample 2>&1
	cd $workdir
	sleep 48h

	# FitHiChIP pipeline
	OutDir=$mapDir/$PREFIX
	BWADir=$OutDir'/Alignment_MAPQ'$MAPQ_Thr
	SegDir=$OutDir'/Segments_HiChIP'
	mkdir -p $BWADir
	mkdir -p $SegDir
	BWA_R1_Alignfile=$BWADir/bwa_out_R1_Chimeric.sam
	BWA_R2_Alignfile=$BWADir/bwa_out_R2_Chimeric.sam
	PairedFilteredReadFile=$BWADir/$PREFIX.paired.cis.RE.filtered.bam
	SortPrefix=$BWADir/$PREFIX.paired.cis.RE.filtered.sorted
	rmDupFile=$BWADir/$PREFIX.paired.cis.RE.filtered.sorted.nodup.bam
	ShortReadFile=$SegDir'/'$PREFIX'.cis.short.bam'
	LongReadFile=$SegDir'/'$PREFIX'.cis.long.bam'
	DumpedPairsSAMFile=$SegDir'/'$PREFIX'.cis.DumpedPairs.sam'
	DumpedPairsBAMFile=$SegDir'/'$PREFIX'.cis.DumpedPairs.bam'
	# bwa mem -t $THREADS $GENOME $FASTQ1 -o ${BWA_R1_Alignfile}.tmp.sam 2>$BWADir/$PREFIX.R1.bwa.log
	# bwa mem -t $THREADS $GENOME $FASTQ2 -o ${BWA_R2_Alignfile}.tmp.sam 2>$BWADir/$PREFIX.R2.bwa.log
	# chimeric.pl ${BWA_R1_Alignfile}.tmp.sam > $BWA_R1_Alignfile &
	# chimeric.pl ${BWA_R2_Alignfile}.tmp.sam > $BWA_R2_Alignfile &
	# wait
	bwa mem -t $THREADS $GENOME $FASTQ1 2>$BWADir/$PREFIX.R1.bwa.log | chimeric.pl - > $BWA_R1_Alignfile
	bwa mem -t $THREADS $GENOME $FASTQ2 2>$BWADir/$PREFIX.R2.bwa.log | chimeric.pl - > $BWA_R2_Alignfile
	pair_up_filter_HiChIP.py -f $BWA_R1_Alignfile -r $BWA_R2_Alignfile -q $MAPQ_Thr -W $DumpedPairsSAMFile 2>$BWADir/$PREFIX.PairFilterlog | samtools view -@ $THREADS -S -h -L $CUT_ENZ - | samtools view -@ $THREADS -bhS - > $PairedFilteredReadFile
	# pair_up_filter_HiChIP.py -f $BWA_R1_Alignfile -r $BWA_R2_Alignfile -q $MAPQ_Thr -W $DumpedPairsSAMFile 2>$BWADir/$PREFIX.PairFilterlog | samtools view -@ $THREADS -bhS - > $PairedFilteredReadFile
	samtools sort -@ $THREADS -o $SortPrefix'.bam' $PairedFilteredReadFile
	picard MarkDuplicates INPUT=$SortPrefix'.bam' OUTPUT=$rmDupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$BWADir/$PREFIX.picard_metrics.txt
	samtools index -@ $THREADS $rmDupFile
	bam2pairs $rmDupFile $BWADir'/'$PREFIX'.paired.cis.RE'
	samtools view -@ $THREADS -h $rmDupFile | awk -v dt="$ShortReadDistThr" 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) < dt) print }' | samtools view -@ $THREADS -bhS - > $ShortReadFile
	samtools index -@ $THREADS $ShortReadFile &
	samtools view -@ $THREADS -h $rmDupFile | awk -v var="$LongReadDistThr" 'function abs(v) {return v < 0 ? -v : v} { if(substr($1, 1, 1)=="@" || abs($9) >= var) print }' | samtools view -@ $THREADS -bhS - > $LongReadFile
	samtools index -@ $THREADS $LongReadFile &
	wait
	bam2pairs $LongReadFile $SegDir'/'$PREFIX'.cis.long' &
	bam2pairs $ShortReadFile $SegDir'/'$PREFIX'.cis.short' &
	wait
	samtools view -@ $THREADS -bhS $DumpedPairsSAMFile > $DumpedPairsBAMFile
	rm $DumpedPairsSAMFile
	macs2dir=$OutDir'/PeaksAnchors_ShortSegment_Fastq_Latest'
	mkdir -p $macs2dir
	macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $ShortReadFile -q 0.01
	grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -S 50% --parallel=200 -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed &
	macs2dir=$OutDir'/PeaksAnchors_CompleteAlignment_Fastq_Latest'
	mkdir -p $macs2dir
	macs2 callpeak -f AUTO -g $GSIZE --outdir $macs2dir -n $PREFIX --nomodel --extsize 147 --call-summits -t $rmDupFile -q 0.01
	grep -vwE "chrM" $macs2dir/$PREFIX'_peaks.narrowPeak' | sort -S 50% --parallel=200 -k5,5rn - | cut -f1-5 - > $macs2dir/$PREFIX.anchors.bed &
	wait
	sed "s?/path/to/bam/inp.bam?$rmDupFile?;s?/path/to/bam/inp.narrowPeak?$macs2dir/$PREFIX.anchors.bed?;s?/path/to/OutDir?$OutDir?;s?PREFIX=FitHiChIP?PREFIX=FitHiChIP.$PREFIX?" $configtemp > $PREFIX.config
	FitHiChIP.sh -C $PREFIX.config &

	# hichipper needs hicpro_output
	# cd $workdir
	################ relative directory
	echo "peaks:" > hichipper.$PREFIX.yaml
	# echo "  - $macs2dir/$PREFIX.anchors.bed" >> hichipper.$PREFIX.yaml
	# echo "  - $macs2dir/${PREFIX}_peaks.narrowPeak" >> hichipper.$PREFIX.yaml
	echo "  - map.res/$PREFIX/PeaksAnchors_CompleteAlignment_Fastq_Latest/${PREFIX}_peaks.narrowPeak" >> hichipper.$PREFIX.yaml
	echo "resfrags:" >> hichipper.$PREFIX.yaml
	echo "  - RestrictionFragmentFiles/hg19_MboI_resfrag.bed.gz" >> hichipper.$PREFIX.yaml
	echo "hicpro_output:" >> hichipper.$PREFIX.yaml
	echo "  - HiCPro.output/$sra" >> hichipper.$PREFIX.yaml
	hichipper --out hichipper.$PREFIX --read-length 76 hichipper.$PREFIX.yaml &
done < $inputfiles
wait
