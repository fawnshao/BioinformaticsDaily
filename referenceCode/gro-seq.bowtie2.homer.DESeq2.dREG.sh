#!/bin/bash
# use bowtie2 to mapping reads
# use dREG to calling peaks de novo 
# https://github.com/Danko-Lab/dREG
# use DESeq2 to find the DEG
#############################################
## set the input and output
HG19=/home1/04935/shaojf/scratch/bowtie2-index/hg19
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes

workdir=/home1/04935/shaojf/scratch/TIP60.project/test.gro.pipes
fastq_dir=$workdir/fastqs

fastqc_dir=$workdir/fastqc_res
clean_dir=$workdir/fastqs_clean
mapping_dir=$workdir/bowtie2_res
homer_dir=$workdir/homer_res
visdir=$workdir/bigwigs
dreg_dir=$workdir/dREG_res
url=http://sjf.dingding.biz/trackhubs
expname=FX.GRO-seq
#############################################

#############################################
## Preprocess data.  Remove adapters.  Trim.
## bowtie2 mapping
## homer
echo " "
echo "Preprocessing and mapping fastq files:"
mkdir $fastqc_dir
mkdir $clean_dir
mkdir $mapping_dir
mkdir $homer_dir
mkdir $visdir
mkdir $dreg_dir

for fastq in `ls $fastq_dir`
do
	name=`echo $fastq | awk -F"/" '{print $NF}' | cut -d "." -f 1`
	cutadapt --nextseq-trim=20 -a polyA=AAAAAAAAAAAAAAAAAAA -a truseq=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a truseqrev=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -m 18 -e 0.10 -o ${clean_dir}/${name}.clean.fastq.gz ${fastq_dir}/$fastq
	bowtie2 -5 3 -3 1 -p 68 -x $HG19 -U ${clean_dir}/${name}.clean.fastq.gz -S $mapping_dir/${name}.sam
	samtools view -q 10 -@ 68 $mapping_dir/${name}.sam | samtools sort -@ 68 -o $mapping_dir/${name}.sorted.bam
	makeTagDirectory $homer_dir/${name}.mTD -tbp 3 -fragLength 200 $mapping_dir/${name}.sorted.bam
done
makeMultiWigHub.pl expname hg19 -url $url -webdir $visdir -d $homer_dir/*.mTD -force -strand &
#############################################

#############################################
## Write out the bigWigs for dREG.
echo " "
echo "Writing bigWigs:"
for f in $mapping_dir/${name}.sorted.bam
do
	j=`echo $f | awk -F"/" '{print $NF}' | cut -d "." -f 1`
	echo $j

	bedtools bamtobed -i $f | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | gzip > ${dreg_dir}/$j.bed.gz
	echo 'Number of mappable reads:'
	echo `zcat ${dreg_dir}/$j.bed.gz | grep "" -c`

	## Convert to bedGraph ... Can't gzip these, unfortunately.
	bedtools genomecov -bg -i ${dreg_dir}/$j.bed.gz -g $genomesize -strand + > ${dreg_dir}/${j}_plus.bedGraph &
	bedtools genomecov -bg -i ${dreg_dir}/$j.bed.gz -g $genomesize -strand - > ${dreg_dir}/${j}_minus.noinv.bedGraph &
	wait

	## Invert minus strand.
	cat ${dreg_dir}/${j}_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${dreg_dir}/${j}_minus.bedGraph ## Invert read counts on the minus strand.

	## sort
	LC_COLLATE=C sort -k1,1 -k2,2n ${dreg_dir}/${j}_plus.bedGraph > ${dreg_dir}/${j}_plus.sorted.bedGraph &
	LC_COLLATE=C sort -k1,1 -k2,2n ${dreg_dir}/${j}_minus.bedGraph > ${dreg_dir}/${j}_minus.sorted.bedGraph &
	wait

	## Then to bigWig
	bedGraphToBigWig ${dreg_dir}/${j}_plus.sorted.bedGraph $genomesize ${dreg_dir}/${j}_plus.bw &
	bedGraphToBigWig ${dreg_dir}/${j}_minus.sorted.bedGraph $genomesize ${dreg_dir}/${j}_minus.bw &
	wait

	rm ${dreg_dir}/$j.bed.gz ${dreg_dir}/$j*.bedGraph
 done
#############################################


#############################################
# dREG
dregmodel=/home1/04935/shaojf/myTools/dREG/dREG-Model/asvm.gdm.6.6M.20170828.rdata
for gro in `ls ${dreg_dir}/${j}_plus.bw | awk -F"/" '{print $NF}' | cut -f 1 -d "."`
do
	negbw=${dreg_dir}/${gro}_minus.bw
	posbw=${dreg_dir}/${gro}_plus.bw
	bash $HOME/myTools/dREG/run_peakcalling.bsh $posbw $negbw ${dreg_dir}/${gro}.dREG.peak $dregmodel 68
	bash $HOME/myTools/dREG/run_dREG.bsh $posbw $negbw $gro.dREG ${dreg_dir}/${gro}.score $dregmodel 68
	bash $HOME/myTools/dREG/writeBed.bsh 0.25 ${dreg_dir}/${gro}.dREG.score.bedGraph.gz
done
#############################################

