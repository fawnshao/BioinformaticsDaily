#!/bin/sh
# add #!/usr/bin/perl, so we can call *.pl directly, instead of perl *.pl
# cp /home1/04935/shaojf/myTools/eCLIP/YeoLab.scripts/gscripts-1.1/perl_scripts/*  ~/.local/bin/
# change perl *.pl to *.pl in Peak_input_normalization_wrapper.pl
##### note: put every file in current working directory. don't move.
##### in the end, move the corresponding files into sub directories.
# prerequirments:
# fastqc, cutadapt, samtools, STAR, bedSort, bedGraphToBigWig

### input experiment name
inputfile=MCF7-CLIP-input
IPA=MCF7-CLIP-RBM39-IPA
IPB=MCF7-CLIP-RBM39-IPB

barcodesfile=/home1/04935/shaojf/scratch/eCLIP_MCF-7/all.barcodes.txt
# adapterfile=/home1/04935/shaojf/scratch/eCLIP_MCF-7/adapter.to.cut.txt

### sub directories, use absolute path
fastq_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqs
fastq_demux=/home1/04935/shaojf/scratch/eCLIP_MCF-7/demux_paired_end_res
fastqc_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqc_res
# trimadapter=/home1/04935/shaojf/myTools/Trimmomatic-0.36/trimmomatic-0.36.jar
trimadapterdir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqs_trimmed
rmrep=/home1/04935/shaojf/scratch/eCLIP_MCF-7/STAR_rmRep
mappingtogenome=/home1/04935/shaojf/scratch/eCLIP_MCF-7/STAR_genome
bw_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/bigwigs
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
clippeaks=/home1/04935/shaojf/scratch/eCLIP_MCF-7/clipper_peaks
normalizedpeaks=/home1/04935/shaojf/scratch/eCLIP_MCF-7/norm_peaks

# # for cutadapt, multiple processors
# module load python3

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
mkdir $fastq_demux
# NNNNNCCTATAT	X1A
# NNNNNTGCTATT	X1B
# NNNNNTATACTT	X2A
# NNNNNATCTTCT	X2B
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

## step 2, cut adatper
date
echo "step 2, cut adatper ......"
mkdir $trimadapterdir
pre=$inputfile
ln -s $fastq_demux/${pre}_R1.unassigned.fastq.gz $trimadapterdir/${pre}_R1.raw.fastq.gz
ln -s $fastq_demux/${pre}_R1.unassigned.fastq.gz $trimadapterdir/${pre}_R2.raw.fastq.gz
pre=$IPA
ln -s $fastq_demux/${pre}_R1.X1A.fastq.gz $trimadapterdir/${pre}_X1A_R1.raw.fastq.gz
ln -s $fastq_demux/${pre}_R2.X1A.fastq.gz $trimadapterdir/${pre}_X1A_R2.raw.fastq.gz
ln -s $fastq_demux/${pre}_R1.X1B.fastq.gz $trimadapterdir/${pre}_X1B_R1.raw.fastq.gz
ln -s $fastq_demux/${pre}_R2.X1B.fastq.gz $trimadapterdir/${pre}_X1B_R2.raw.fastq.gz
pre=$IPB
ln -s $fastq_demux/${pre}_R1.X2A.fastq.gz $trimadapterdir/${pre}_X2A_R1.raw.fastq.gz
ln -s $fastq_demux/${pre}_R2.X2A.fastq.gz $trimadapterdir/${pre}_X2A_R2.raw.fastq.gz
ln -s $fastq_demux/${pre}_R1.X2B.fastq.gz $trimadapterdir/${pre}_X2B_R1.raw.fastq.gz
ln -s $fastq_demux/${pre}_R2.X2B.fastq.gz $trimadapterdir/${pre}_X2B_R2.raw.fastq.gz
	# -G AATGATACGGCGACCACCGAGATCTCTCTTTCCCTACACGACGCTCTTCCGATCT \
	# -A CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT \
	# -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	# -g CTACACGACGCTCTTCCGATCT \
	# -G CTACACGACGCTCTTCCGATCT \
	# -A rand103Tr3=NNNNNNNNNNAGATCGGAAGAGCACACGTCTG \
	# -A RiL19=AGATCGGAAGAGCGTCGTGTAG \
	# -A X1A=^NNNNNCCTATAT \
	# -A X1B=^NNNNNTGCTATT \
	# -A X2A=^NNNNNTATACTT \
	# -A X2B=^NNNNNATCTTCT \
	# -A AGATCGGAAGAGCGTCGTGTAG \
# cutadapt searches and removes only the sequence as you provide it. The reverse-complement is not searched. This is because usually the two adapters used on both sides can be different from each other. If the adapters are the same ones in your case, you will need to reverse-complement the sequence yourself and provide it to cutadapt with a second -a option. Cutadapt will then remove the best-matching adapter from each read.
for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
do
	cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 \
	--quality-cutoff 6 -m 18 --nextseq-trim=20 \
	-g CTACACGACGCTCTTCCGATCT \
	-G CTACACGACGCTCTTCCGATCT \
	-A NNNNNNNNNNAGATCGGAAGAGCACACGTCTG \
	-A NNNNNCCTATAT \
	-A NNNNNTGCTATT \
	-A NNNNNTATACTT \
	-A NNNNNATCTTCT \
	-A CAGACGTGTGCTCTTCCGATCTNNNNNNNNNN \
	-A ATATAGGNNNNN \
	-A AATAGCANNNNN \
	-A AAGTATANNNNN \
	-A AGAAGATNNNNN \
	-o $trimadapterdir/${pre}_R1.clean.fastq.gz \
	-p $trimadapterdir/${pre}_R2.clean.fastq.gz \
	$trimadapterdir/${pre}_R1.raw.fastq.gz $trimadapterdir/${pre}_R2.raw.fastq.gz \
	> $trimadapterdir/${pre}.metrics &
done
wait

## step 3, STAR rmRep
date
echo "step 3, STAR rmRep ......"
mkdir $rmrep
# STAR rmRep: Takes output from cutadapt round 2. Maps to human specific version of RepBase used to remove repetitive elements, helps control for spurious artifacts from rRNA (& other) repetitive reads.
# remove --outSAMunmapped Within --outStd BAM_Tnsorted > $rmrep/${pre}.round2.rep.bam
for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $repindexdir --genomeLoad LoadAndRemove \
	--readFilesIn $trimadapterdir/${pre}_R1.clean.fastq.gz \
	$trimadapterdir/${pre}_R2.clean.fastq.gz \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 30 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix $rmrep/${pre}.toRep. \
	--outSAMattributes All --readFilesCommand zcat \
	--outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local
	# --alignEndsType EndToEnd

	samtools view $rmrep/${pre}.toRep.Aligned.out.bam | count_aligned_from_sam.py > $rmrep/${pre}.toRep.metrics
	# fastqc --outdir $fastqc_dir $rmrep/${pre}.toRep.Unmapped.out.mate1 > $rmrep/${pre}.toRep.Unmapped.out.mate1.dummy_fastqc &
	# fastqc --outdir $fastqc_dir $rmrep/${pre}.toRep.Unmapped.out.mate2 > $rmrep/${pre}.toRep.Unmapped.out.mate2.dummy_fastqc &
	fastq-sort --id $rmrep/${pre}.toRep.Unmapped.out.mate1 > $rmrep/${pre}.toRep.Unmapped.out.sorted.mate1 &
	fastq-sort --id $rmrep/${pre}.toRep.Unmapped.out.mate2 > $rmrep/${pre}.toRep.Unmapped.out.sorted.mate2 &
done
wait

## step 4, STAR genome
date
echo "step 4, STAR genome ......"
mkdir $mappingtogenome
# --outSAMunmapped Within \
for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn $rmrep/${pre}.toRep.Unmapped.out.sorted.mate1 \
	$rmrep/${pre}.toRep.Unmapped.out.sorted.mate2 \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
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
for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
do
	barcode_collapse_pe.py --bam ${pre}.rmRep.genome.Aligned.out.bam --out_file ${pre}.rmRep.genome.rmDup.bam --metrics_file ${pre}.rmRep.genome.rmDup.metrics
	picard SortSam INPUT=${pre}.rmRep.genome.rmDup.bam OUTPUT=${pre}.rmRep.genome.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
done
samtools merge -@ 68 ${IPA}.rmRep.genome.sorted.bam ${IPA}_X1A.rmRep.genome.sorted.bam ${IPA}_X1B.rmRep.genome.sorted.bam
samtools merge -@ 68 ${IPB}.rmRep.genome.sorted.bam ${IPB}_X2A.rmRep.genome.sorted.bam ${IPB}_X2B.rmRep.genome.sorted.bam
for pre in $inputfile $IPA $IPB
do
	samtools index -@ 68 ${pre}.rmRep.genome.sorted.bam
	samtools view -hb -f 128 -@ 68 ${pre}.rmRep.genome.sorted.bam > ${pre}.rmRep.genome.sorted.r2.bam
	# samtools sort -@ 68 -o $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.bam $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam
	# samtools index -@ 68 $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.bam
	make_bigwig_files.py --bam ${pre}.rmRep.genome.sorted.r2.bam --genome $genomesize --bw_pos ${pre}.rmRep.genome.sorted.r2.norm.pos.bw --bw_neg ${pre}.rmRep.genome.sorted.r2.norm.neg.bw
	bedSort ${pre}.rmRep.genome.sorted.r2.norm.pos.bg ${pre}.rmRep.genome.sorted.r2.norm.pos.srt.bg &
	bedSort ${pre}.rmRep.genome.sorted.r2.norm.neg.t.bg ${pre}.rmRep.genome.sorted.r2.norm.neg.t.srt.bg &
	wait
	bedGraphToBigWig ${pre}.rmRep.genome.sorted.r2.norm.pos.srt.bg $genomesize bedGraphToBigWig.${pre}.rmRep.genome.sorted.r2.norm.pos.bw &
	bedGraphToBigWig ${pre}.rmRep.genome.sorted.r2.norm.neg.t.srt.bg $genomesize bedGraphToBigWig.${pre}.rmRep.genome.sorted.r2.norm.neg.bw &
done

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
echo | awk -vOFS="\t" '{print "001","RBM39","MCF-7","ABCD1","ABCD2","ABCD3"}' >> $IPA.$IPB.manifest.txt
sed -i "s?ABCD1?${IPA}.rmRep.genome.sorted.r2.bam?;s?ABCD2?${IPB}.rmRep.genome.sorted.r2.bam?;s?ABCD3?$inputfile.rmRep.genome.sorted.bam?" $IPA.$IPB.manifest.txt
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

genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
for f in *.bed
do
	new=`echo $f | sed 's/bed$/bigbed/'`
	bedSort $f a.bed
	awk -vOFS="\t" '{print $1,$2,$3,$4"|"$5}' a.bed > b.bed
	bedToBigBed -type=bed4 b.bed $genomesize $new
	rm a.bed b.bed
done

