#!/bin/sh
# add #!/usr/bin/perl, so we can call *.pl directly, instead of perl *.pl
# cp /home1/04935/shaojf/myTools/eCLIP/YeoLab.scripts/gscripts-1.1/perl_scripts/*  ~/.local/bin/
# change perl *.pl to *.pl in Peak_input_normalization_wrapper.pl
##### note: put every file in current working directory. don't move.
##### in the end, move the corresponding files into sub directories.
# prerequirments:
# fastqc, cutadapt, samtools, STAR, bedSort, bedGraphToBigWig

# rclone sync bigwigs/ mygoogle:CLIP/eCLIP_ENCODE_testpipeline/bigwigs
### input experiment name
rbp=RBFOX2
cell=K562
pre=K562.RBFOX2
inputfile=$pre.input
IPA=$pre.IPA
IPB=$pre.IPB

workdir=/home1/04935/shaojf/scratch/eCLIP.summary
barcodesfile=$workdir/all.barcodes.txt
# AAGCAAT	A01
# GGCTTGT	B06
# ACAAGTT	C01
# TGGTCCT	D08
# ATGACCNNNNT	A03
# TCCTGTNNNNT	G07
# CAGCTTNNNNT	A04
# GGATACNNNNT	F05
# NNNNNCCTATAT	X1A
# NNNNNTGCTATT	X1B
# NNNNNTATACTT	X2A
# NNNNNATCTTCT	X2B
### sub directories, use absolute path
fastq_dir=$workdir/fastqs
fastq_demux=$workdir/demux_paired_end_res
fastqc_dir=$workdir/fastqc_res
trimadapterdir=$workdir/fastqs_trimmed
rmrep=$workdir/STAR_rmRep
mappingtogenome=$workdir/STAR_genome
bw_dir=$workdir/bigwigs
clippeaks=$workdir/clipper_peaks
normalizedpeaks=$workdir/norm_peaks

# # for cutadapt, multiple processors
# module load python3
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
### for STAR
hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
hg19indexdir=/home1/04935/shaojf/scratch/star_index/hg19.star
repfile=/home1/04935/shaojf/scratch/star_index/RepBase23.01.human.fa
repindexdir=/home1/04935/shaojf/scratch/star_index/RepBase23.01.human.star
### STAR index, run once, used in future
# mkdir $hg19indexdir
# mkdir $repindexdir
# STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $hg19indexdir --genomeFastaFiles $hg19file
# STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $repindexdir --genomeFastaFiles $repfile

## step 1, FastQC
date
echo "step 1, FastQC ......"
mkdir $fastqc_dir
fastqc --threads 68 --outdir $fastqc_dir $fastq_dir/*.fastq* &
## Demultiplexing
# must do it; the random-mer was appended to the read name for later usage. 
# do remove PCR duplication?
# @NS500669:202:HHG7GBGX3:1:11101:5225:1082 1:N:0:TTAGGC
# to
# @NCCACCCGAA:NS500669:202:HHG7GBGX3:1:11101:5225:1082 1:N:0:TTAGGC

mkdir $fastq_demux
# gunzip -c fastqs/K562.RBFOX2.IPA_R1.fastq.gz | sed -n '2~4p' | cut -c 1-6 | head -1000 | sort | more
for pre in $inputfile $IPA $IPB
do
	demux_paired_end.py --fastq_1 $fastq_dir/${pre}_R1.fastq.gz --fastq_2 $fastq_dir/${pre}_R2.fastq.gz \
	-b $barcodesfile --out_file_1 ${pre}_R1.fastq.gz --out_file_2 ${pre}_R2.fastq.gz \
	--length 10 -m $fastq_demux/${pre}.metrics &
done
wait
mv ${inputfile}_R?.*.fastq.gz $fastq_demux/
mv ${IPA}_R?.*.fastq.gz $fastq_demux/
mv ${IPB}_R?.*.fastq.gz $fastq_demux/

####### merged
mkdir $trimadapterdir
for pre in $inputfile $IPA $IPB
do
	zcat $fastq_demux/${pre}_R1.*.fastq.gz | gzip >  $trimadapterdir/${pre}_R1.demux.fastq.gz &
	zcat $fastq_demux/${pre}_R2.*.fastq.gz | gzip >  $trimadapterdir/${pre}_R2.demux.fastq.gz &
done
wait

## step 2, cut adatper
date
echo "step 2, cut adatper ......"
# cutadapt searches and removes only the sequence as you provide it. The reverse-complement is not searched. This is because usually the two adapters used on both sides can be different from each other. If the adapters are the same ones in your case, you will need to reverse-complement the sequence yourself and provide it to cutadapt with a second -a option. Cutadapt will then remove the best-matching adapter from each read.
for pre in $inputfile $IPA $IPB
do
	cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 \
	--quality-cutoff 6 -m 18 --nextseq-trim=20 \
	-g r15p=CTACACGACGCTCTTCCGATCT \
	-G r25p=CTACACGACGCTCTTCCGATCT \
	-a rand103Tr3=NNNNNNNNNNAGATCGGAAGAGCACACGTCTG \
	-A comp_rand103Tr3=CAGACGTGTGCTCTTCCGATCTNNNNNNNNNN \
	-o $trimadapterdir/${pre}_R1.clean.fastq.gz \
	-p $trimadapterdir/${pre}_R2.clean.fastq.gz \
	$trimadapterdir/${pre}_R1.demux.fastq.gz $trimadapterdir/${pre}_R2.demux.fastq.gz \
	> $trimadapterdir/${pre}.metrics &
done
wait

## step 3, STAR rmRep
date
echo "step 3, STAR rmRep ......"
mkdir $rmrep
# STAR rmRep: Takes output from cutadapt round 2. Maps to human specific version of RepBase used to remove repetitive elements, helps control for spurious artifacts from rRNA (& other) repetitive reads.
# remove --outSAMunmapped Within --outStd BAM_Tnsorted > $rmrep/${pre}.round2.rep.bam
# for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
for pre in $inputfile $IPA $IPB
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $repindexdir --genomeLoad LoadAndRemove \
	--readFilesIn $trimadapterdir/${pre}_R1.clean.fastq.gz \
	$trimadapterdir/${pre}_R2.clean.fastq.gz \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 30 \
	--outFilterMultimapScoreRange 1 \
	--outFilterScoreMinOverLread  0.5 \
	--outFilterMatchNminOverLread 0.5 \
	--outFileNamePrefix $rmrep/${pre}.toRep. \
	--outSAMattributes All --readFilesCommand zcat \
	--outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local

	samtools view $rmrep/${pre}.toRep.Aligned.out.bam | count_aligned_from_sam.py > $rmrep/${pre}.toRep.metrics
	fastq-sort --id $rmrep/${pre}.toRep.Unmapped.out.mate1 > $rmrep/${pre}.toRep.Unmapped.out.sorted.mate1 &
	fastq-sort --id $rmrep/${pre}.toRep.Unmapped.out.mate2 > $rmrep/${pre}.toRep.Unmapped.out.sorted.mate2 &
done
wait

## step 4, STAR genome
date
echo "step 4, STAR genome ......"
mkdir $mappingtogenome
# "too short" means that the best alignments STAR found were too short to pass the filters.
# This is controlled by --outFilterScoreMinOverLread --outFilterMatchNminOverLread which by default are set to 0.66. which means that ~2/3 of the total read length (sum of mates) should be mapped.
# So what your output means is that 12% of the reads aligned unniquely, 7.7% aligned but multimapped and then 80% of your reads couldn't align with the above parameters. You can try to reduce these parameters to see how many more reads will be mapped. However, it looks like your data might just be contaminated with that alignment score :((
# --outSAMunmapped Within \
# for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
for pre in $inputfile $IPA $IPB
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn $rmrep/${pre}.toRep.Unmapped.out.sorted.mate1 \
	$rmrep/${pre}.toRep.Unmapped.out.sorted.mate2 \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFilterScoreMinOverLread  0.3 \
	--outFilterMatchNminOverLread 0.3 \
	--outFileNamePrefix $mappingtogenome/${pre}.rmRep.genome. \
	--outSAMattributes All --outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local
	# --alignEndsType EndToEnd
done

## step 5, process STAR results
# Barcode_collapse_pe: takes output from STAR genome mapping. Custom random-mer-aware script for PCR duplicate removal
date
echo "step 5, process STAR results ......"
mkdir $bw_dir
cd $mappingtogenome
# for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
for pre in $inputfile $IPA $IPB
do
	barcode_collapse_pe.py --bam ${pre}.rmRep.genome.Aligned.out.bam --out_file ${pre}.rmRep.genome.rmDup.bam --metrics_file ${pre}.rmRep.genome.rmDup.metrics
	picard SortSam INPUT=${pre}.rmRep.genome.rmDup.bam OUTPUT=${pre}.rmRep.genome.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
done
for pre in $inputfile $IPA $IPB
do
	samtools index -@ 68 ${pre}.rmRep.genome.sorted.bam
	samtools view -hb -f 128 -@ 68 ${pre}.rmRep.genome.sorted.bam > ${pre}.rmRep.genome.sorted.r2.bam
	samtools index -@ 68 $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam
	make_bigwig_files.py --bam ${pre}.rmRep.genome.sorted.r2.bam --genome $genomesize --bw_pos ${pre}.rmRep.genome.sorted.r2.norm.pos.bw --bw_neg ${pre}.rmRep.genome.sorted.r2.norm.neg.bw
	bedSort ${pre}.rmRep.genome.sorted.r2.norm.pos.bg ${pre}.rmRep.genome.sorted.r2.norm.pos.srt.bg &
	bedSort ${pre}.rmRep.genome.sorted.r2.norm.neg.t.bg ${pre}.rmRep.genome.sorted.r2.norm.neg.t.srt.bg &
	wait
	bedGraphToBigWig ${pre}.rmRep.genome.sorted.r2.norm.pos.srt.bg $genomesize bedGraphToBigWig.${pre}.rmRep.genome.sorted.r2.norm.pos.bw &
	bedGraphToBigWig ${pre}.rmRep.genome.sorted.r2.norm.neg.t.srt.bg $genomesize bedGraphToBigWig.${pre}.rmRep.genome.sorted.r2.norm.neg.bw &
done
# more ~/.local/lib/python2.7/site-packages/gscripts-0.1.6-py2.7.egg/gscripts/general/make_bigwig_files.py
# samtools view -h " $bamfile | awk 'BEGIN {OFS="\t"} {if(!!and($2,0x0080)) {if(!!and($2, 0x0004
# )) {$2 = $2 - 16} else {$2 = $2 + 16}}; print $0}' | samtools view -bS - | genomeCoverageBed -ibam stdin 

## step 6, calling peaks with clipper
# Clipper: Takes results from samtools view. Calls peaks on those files
# rm "--bonferroni --superlocal --threshold-method binomial "
date
echo "step 6, calling peaks with clipper ......"
mkdir $clippeaks
for pre in $inputfile $IPA $IPB
do
	clipper -b ${pre}.rmRep.genome.sorted.r2.bam -s hg19 -o ${pre}.rmRep.genome.sorted.r2.peaks.bed --save-pickle
	fix_scores.py --bed ${pre}.rmRep.genome.sorted.r2.peaks.bed --out_file ${pre}.rmRep.genome.sorted.r2.peaks.fixed.bed
	bedToBigBed ${pre}.rmRep.genome.sorted.r2.peaks.fixed.bed $genomesize ${pre}.rmRep.genome.sorted.r2.peaks.fixed.bb -type=bed6+4 &
done

## step 7, peak normalization
date
echo "step 7, peak normalization ......"
mkdir $normalizedpeaks
# uID	RBP	Cell line	CLIP_rep1	CLIP_rep2	INPUT
# 001	RBM39	MCF-7	${IPA}.rmRep.genome.sorted.r2.bam	${IPA}.rmRep.genome.sorted.r2.bam	$inputfile.rmRep.genome.sorted.bam
echo | awk -vOFS="\t" '{print "uID","RBP","Cell line","CLIP_rep1","CLIP_rep2","INPUT"}' > $IPA.$IPB.manifest.txt
echo "001 $rbp $cell ${IPA}.rmRep.genome.sorted.r2.bam ${IPB}.rmRep.genome.sorted.r2.bam $inputfile.rmRep.genome.sorted.bam" | awk -vOFS="\t" '{print $0}' >> $IPA.$IPB.manifest.txt
Peak_input_normalization_wrapper.pl $IPA.$IPB.manifest.txt $normalizedpeaks
wait
mv *.bw $bw_dir/
mv *r2.peaks* $clippeaks/
# peaks.l2inputnormnew.bed
# Chr \t start \t stop \t log10(p-value eCLIP vs SMInput) \t log2(fold-enrichmentin eCLIP vs SMInput) \t strand
# peaks.l2inputnormnew.bed.full
# Chr \t start \t stop \t peakID:CLIPer pvalue \t **** \t log10(p-value eCLIP vs SMInput) \t log2(fold-enrichmentin eCLIP vs SMInput)



### conda install -c bioconda PureCLIP
# hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
# peak: 61.2% memory
# Run time 11:31:24
for pre in $IPA $IPB
do
	pureclip -i ${pre}.rmRep.genome.sorted.r2.bam -bai ${pre}.rmRep.genome.sorted.r2.bam.bai -g $hg19file -o ${pre}.called_crosslinksites.bed -or ${pre}.called_bindingregions.bed -p ${pre}.learnedparameters.txt -ibam $inputfile.rmRep.genome.sorted.bam -ibai $inputfile.rmRep.genome.sorted.bam.bai -nt 68 -tmp ./${pre}.tmp/ &
done
wait

### conda install -c bioconda piranha
# peak: 2.2% memory
# Run time 00:12:50
for pre in $inputfile $IPA $IPB
do
	Piranha -o ${pre}.piranha.bed -b 10 -s ${pre}.rmRep.genome.sorted.r2.bam &
done
wait

for f in *.bed
do
	new=`echo $f | sed 's/bed$/bigbed/'`
	bedSort $f a.bed
	awk -vOFS="\t" '{print $1,$2,$3,$4"|"$5}' a.bed > b.bed
	bedToBigBed -type=bed4 b.bed $genomesize $new
	rm a.bed b.bed
done

