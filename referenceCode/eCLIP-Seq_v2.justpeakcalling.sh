#!/bin/sh
# add #!/usr/bin/perl, so we can call *.pl directly, instead of perl *.pl
# cp /home1/04935/shaojf/myTools/eCLIP/YeoLab.scripts/gscripts-1.1/perl_scripts/*  ~/.local/bin/
# change perl *.pl to *.pl in Peak_input_normalization_wrapper.pl
##### note: put every file in current working directory. don't move.
##### in the end, move the corresponding files into sub directories.
# prerequirments:
# fastqc, cutadapt, samtools, STAR, bedSort, bedGraphToBigWig

# rclone sync bigwigs/ mygoogle:CLIP/eCLIP_ENCODE_testpipeline/bigwigs

# only for mapping and visualization
### input experiment name
# AllInputs="TDP43.IPA.WT TDP43.IPA.Clone51 TDP43.IPA.Clone48 smINPUT1.WT smINPUT.Clone51 smINPUT.Clone48 TDP43.IPB.WT TDP43.IPB.Clone50 TDP43.IPB.Clone48 smINPUT2.WT"
AllInputs="smINPUT.WT TDP43.IPA.Clone48 smINPUT.Clone51 smINPUT.Clone48 TDP43.IPB.WT TDP43.IPB.Clone50 TDP43.IPB.Clone48"
# samtools merge smINPUT.WT.rmRep.genome.r2.bam smINPUT1.WT.rmRep.genome.sorted.r2.bam smINPUT2.WT.rmRep.genome.sorted.r2.bam
# picard SortSam INPUT=smINPUT.WT.rmRep.genome.r2.bam \
# 		OUTPUT=smINPUT.WT.rmRep.genome.sorted.r2.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
# NS19-Wenbo-eCLIP-IP1_S1	TDP43.IPA.WT
# NS19-Wenbo-eCLIP-IP4_S2	TDP43.IPA.Clone51
# NS19-Wenbo-eCLIP-IP7_S3	TDP43.IPA.Clone48
# NS19-Wenbo-eCLIP-INP1_S4	smINPUT1.WT
# NS19-Wenbo-eCLIP-INP4_S5	smINPUT.Clone51
# NS19-Wenbo-eCLIP-INP7_S6	smINPUT.Clone48
# NS19-Wenbo-eCLIP-IP2_S7	TDP43.IPB.WT
# NS19-Wenbo-eCLIP-IP5_S8	TDP43.IPB.Clone50
# NS19-Wenbo-eCLIP-IP8_S9	TDP43.IPB.Clone48
# NS19-Wenbo-eCLIP-INP2_S10	smINPUT2.WT

workdir=/home1/04935/shaojf/scratch/NS19.eCLIP
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
mappingtogenome=$workdir/STAR_genome
bw_dir=$workdir/bigwigs
clippeaks=$workdir/clipper_peaks
normalizedpeaks=$workdir/norm_peaks

# # for cutadapt, multiple processors
# module load python3
genomesize=/home1/04935/shaojf/scratch/star_index/mm10/mm10.chrom.sizes
### for STAR
hg19file=/home1/04935/shaojf/scratch/star_index/mm10/mm10.fa
hg19indexdir=/home1/04935/shaojf/scratch/star_index/mm10/mm10.star
### STAR index, run once, used in future
# mkdir $hg19indexdir
# mkdir $repindexdir
# STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $hg19indexdir --genomeFastaFiles $hg19file
# STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $repindexdir --genomeFastaFiles $repfile

## step 6, calling peaks with clipper
# Clipper: Takes results from samtools view. Calls peaks on those files
# rm "--bonferroni --superlocal --threshold-method binomial "
# --gtfFile=/home1/04935/shaojf/stampede2/refs/GENCODE/gencode.vM19.annotation.gtf \
# pre=TDP43.IPB.WT
# clipper -b ${pre}.rmRep.genome.sorted.r2.bam \
# 	-s mm10 -g ENSMUSG00000092341 \
# 	-o ${pre}.test.bed --bonferroni \
# 	--superlocal --threshold-method binomial --binomial=0.05 --save-pickle

date
echo "step 6, calling peaks with clipper ......"
mkdir $clippeaks
for pre in $AllInputs
do
	clipper -b ${pre}.rmRep.genome.sorted.r2.bam \
		--gtfFile=/home1/04935/shaojf/stampede2/refs/GENCODE/gencode.vM19.annotation.gtf \
		-o ${pre}.rmRep.genome.sorted.r2.peaks.bed --bonferroni \
		--superlocal --threshold-method binomial --binomial=0.05 --save-pickle
	fix_scores.py --bed ${pre}.rmRep.genome.sorted.r2.peaks.bed --out_file ${pre}.rmRep.genome.sorted.r2.peaks.fixed.bed
	bedToBigBed ${pre}.rmRep.genome.sorted.r2.peaks.fixed.bed $genomesize ${pre}.rmRep.genome.sorted.r2.peaks.fixed.bb -type=bed6+4 &
done

## step 7, peak normalization
#### please be careful, list the full path for bam file
date
echo "step 7, peak normalization ......"
# mkdir $normalizedpeaks
# uID	RBP	Cell line	CLIP_rep1	CLIP_rep2	INPUT
# 001	RBM39	MCF-7	${IPA}.rmRep.genome.sorted.r2.bam	${IPA}.rmRep.genome.sorted.r2.bam	$inputfile.rmRep.genome.sorted.bam
# echo | awk -vOFS="\t" '{print "uID","RBP","Cell line","CLIP_rep1","CLIP_rep2","INPUT"}' > $IPA.$IPB.manifest.txt
# echo "001 $rbp $cell ${IPA}.rmRep.genome.sorted.r2.bam ${IPB}.rmRep.genome.sorted.r2.bam $inputfile.rmRep.genome.sorted.r2.bam" | awk -vOFS="\t" '{print $0}' >> $IPA.$IPB.manifest.txt
# Peak_input_normalization_wrapper.pl $IPA.$IPB.manifest.txt $normalizedpeaks
# awk '$4>3 && $5>=1' 001_01.basedon_001_01.peaks.l2inputnormnew.bed > ${IPA}.significant.peaks.bed
# awk '$4>3 && $5>=1' 001_02.basedon_001_02.peaks.l2inputnormnew.bed > ${IPB}.significant.peaks.bed
# mv ${IPA}.significant.peaks.bed $clippeaks/
# mv ${IPB}.significant.peaks.bed $clippeaks/
# pairs="smINPUT.WT:TDP43.IPB.WT smINPUT.Clone48:TDP43.IPA.Clone48:TDP43.IPB.Clone48 smINPUT.Clone51:TDP43.IPB.Clone50"
IPA=TDP43.IPA.Clone48
IPB=TDP43.IPB.Clone48
inputfile=smINPUT.Clone48
echo | awk -vOFS="\t" '{print "uID","RBP","Cell line","CLIP_rep1","CLIP_rep2","INPUT"}' > $IPA.$IPB.manifest.txt
echo "Clone48 TDP43 mouse ${IPA}.rmRep.genome.sorted.r2.bam ${IPB}.rmRep.genome.sorted.r2.bam $inputfile.rmRep.genome.sorted.r2.bam" | \
		awk -vOFS="\t" '{print $0}' >> $IPA.$IPB.manifest.txt
Peak_input_normalization_wrapper.pl $IPA.$IPB.manifest.txt .
awk '$4>3 && $5>=1' Clone48_01.basedon_Clone48_01.peaks.l2inputnormnew.bed > ${IPA}.significant.peaks.bed
awk '$4>3 && $5>=1' Clone48_02.basedon_Clone48_02.peaks.l2inputnormnew.bed > ${IPB}.significant.peaks.bed

echo | awk -vOFS="\t" '{print "uID","RBP","Cell line","CLIP","INPUT"}' > onerep.manifest.txt
IPA=TDP43.IPB.WT
inputfile=smINPUT.WT
echo "WT TDP43 mouse ${IPA}.rmRep.genome.sorted.r2.bam $inputfile.rmRep.genome.sorted.r2.bam" | \
		awk -vOFS="\t" '{print $0}' >> onerep.manifest.txt
IPA=TDP43.IPB.Clone50
inputfile=smINPUT.Clone51
echo "Clone50 TDP43 mouse ${IPA}.rmRep.genome.sorted.r2.bam $inputfile.rmRep.genome.sorted.r2.bam" | \
		awk -vOFS="\t" '{print $0}' >> onerep.manifest.txt
Peak_input_normalization_wrapper.pl onerep.manifest.txt .
awk '$4>3 && $5>=1' WT_01.basedon_WT_01.peaks.l2inputnormnew.bed > TDP43.WT.significant.peaks.bed
awk '$4>3 && $5>=1' Clone50_01.basedon_Clone50_01.peaks.l2inputnormnew.bed > TDP43.Clone50.significant.peaks.bed



# peaks.l2inputnormnew.bed
# Chr \t start \t stop \t log10(p-value eCLIP vs SMInput) \t log2(fold-enrichmentin eCLIP vs SMInput) \t strand
# peaks.l2inputnormnew.bed.full
# Chr \t start \t stop \t peakID:CLIPer pvalue \t **** \t log10(p-value eCLIP vs SMInput) \t log2(fold-enrichmentin eCLIP vs SMInput)
# wait
# mv $mappingtogenome/*.rmRep.genome.sorted.r2.norm.pos.bw $bw_dir/
# ########################################################################

# ### conda install -c bioconda PureCLIP
# # hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
# # peak: 83.3% memory
# # Run time 11:31:24
# for pre in $IPA $IPB
# do
# 	pureclip -i ${pre}.rmRep.genome.sorted.r2.bam \
# 		-bai ${pre}.rmRep.genome.sorted.r2.bam.bai \
# 		-g $hg19file -o ${pre}.called_crosslinksites.bed \
# 		-or ${pre}.called_bindingregions.bed \
# 		-p ${pre}.learnedparameters.txt \
# 		-ibam $inputfile.rmRep.genome.sorted.bam \
# 		-ibai $inputfile.rmRep.genome.sorted.bam.bai -nt 272 -tmp ./${pre}.tmp/ &
# done
# wait

# ### conda install -c bioconda piranha
# # peak: 2.2% memory
# # Run time 00:12:50
# for pre in $inputfile $IPA $IPB
# do
# 	Piranha -o ${pre}.piranha.bed -b 10 -s ${pre}.rmRep.genome.sorted.r2.bam &
# done
# wait

# for f in *.bed
# do
# 	new=`echo $f | sed 's/bed$/bigbed/'`
# 	bedSort $f a.bed
# 	awk -vOFS="\t" '{print $1,$2,$3,$4"|"$5}' a.bed > b.bed
# 	bedToBigBed -type=bed4 b.bed $genomesize $new
# 	rm a.bed b.bed
# done

