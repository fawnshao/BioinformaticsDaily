#!/bin/sh
# adapters=/home1/04935/shaojf/scratch/eCLIP_MCF-7/adapters.txt
# barcodes=/home1/04935/shaojf/scratch/eCLIP_MCF-7/primers.barcodes.txt
######## eCLIP adapters ########
# >RNA_A01
# AUUGCUUAGAUCGGAAGAGCGUCGUGUAG
# >RNA_B06
# ACAAGCCAGAUCGGAAGAGCGUCGUGUAG
# >RNA_C01
# AACUUGUAGAUCGGAAGAGCGUCGUGUAG
# >RNA_D08
# AGGACCAAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X1A
# AUAUAGGNNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X1B
# AAUAGCANNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X2A
# AAGUAUANNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X2B
# AGAAGAUNNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RiL19
# AGAUCGGAAGAGCGUCGUG
# >AR17
# ACACGACGCTCTTCCGA
# >rand103T3
# NNNNNNNNNNAGATCGGAAGAGCACACGTCTG
# >PCR_F_D501
# AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
# >PCR_F_D502
# AATGATACGGCGACCACCGAGATCTACACATAGAGGCACACTCTTTCCCTACACGACGCTCTTCCGATCT
# >PCR_R_D701
# CAAGCAGAAGACGGCATACGAGATCGAGTAATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
# >PCR_R_D702
# CAAGCAGAAGACGGCATACGAGATTCTCCGGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC

# >Index_1_(i7)_Adapters
# CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
# >Index_2_(i5)_Adapters
# AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
# >Adapter_Trimming
# CTGTCTCTTATACACATCT


### sub directories
fastq_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqs
adapterTrimmed=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqs_adapterTrim
fastqc_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqc_res
rmrep=/home1/04935/shaojf/scratch/eCLIP_MCF-7/STAR_rmRep
mappingtogenome=/home1/04935/shaojf/scratch/eCLIP_MCF-7/STAR_genome
bw_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/bigwigs
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
clippeaks=/home1/04935/shaojf/scratch/eCLIP_MCF-7/clipper_peaks

# for cutadapt, multiple processors
module load python3

### for STAR
hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
hg19indexdir=/home1/04935/shaojf/scratch/star_index/hg19.star
repfile=/home1/04935/shaojf/scratch/star_index/RepBase23.01.human.fa
repindexdir=/home1/04935/shaojf/scratch/star_index/RepBase23.01.human.star
### STAR index, run once, used in future
mkdir $hg19indexdir
mkdir $repindexdir
STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $hg19indexdir --genomeFastaFiles $hg19file
STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $repindexdir --genomeFastaFiles $repfile

## step 1, FastQC
date
echo "step 1, FastQC ......"
mkdir $fastqc_dir
# fastqc --contaminants $barcodes --adapters $adapters --outdir fastqc_res $fastq_dir/*.fastq
fastqc --outdir $fastqc_dir $fastq_dir/*.fastq
### Demultiplexing
## don't need to do this
# for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
# do
# 	mkdir demux_paired_end_res
# 	demux_paired_end.py --fastq_1=$fastq_dir/${pre}_R1.fastq --fastq_2=$fastq_dir/${pre}_R1.fastq \
# 	--out_file_1=demux_paired_end_res/${pre}_R1.fastq --out_file_2=demux_paired_end_res/${pre}_R1.fastq
# done

## step 2, Cutadapt
date
echo "step 2, Cutadapt ......"
mkdir $adapterTrimmed
# Cutadapt round 1: Takes output from demultiplexed files. Run to trim off both 5’ and 3’ adapters on both reads
# with cutadapt v1.15, not use "-f fastq --match-read-wildcards"
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	/home1/04935/shaojf/.local/bin/cutadapt --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 --nextseq-trim=20 --cores=21 -A a0=NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g g1=CTTCCGATCTACAAGTT -g g2=CTTCCGATCTTGGTCCT -A a1=AACTTGTAGATCGGA -A a2=AGGACCAAGATCGGA -A a3=ACTTGTAGATCGGAA -A a4=GGACCAAGATCGGAA -A a5=CTTGTAGATCGGAAG -A a6=GACCAAGATCGGAAG -A a7=TTGTAGATCGGAAGA -A a8=ACCAAGATCGGAAGA -A a9=TGTAGATCGGAAGAG -A a10=CCAAGATCGGAAGAG -A a11=GTAGATCGGAAGAGC -A a12=CAAGATCGGAAGAGC -A a13=TAGATCGGAAGAGCG -A a14=AAGATCGGAAGAGCG -A a15=AGATCGGAAGAGCGT -A a16=GATCGGAAGAGCGTC -A a17=ATCGGAAGAGCGTCG -A a18=TCGGAAGAGCGTCGT -A a19=CGGAAGAGCGTCGTG -A a20=GGAAGAGCGTCGTGT -o $adapterTrimmed/${pre}_R1.adapterTrim.fastq -p $adapterTrimmed/${pre}_R2.adapterTrim.fastq $fastq_dir/${pre}_R1.fastq $fastq_dir/${pre}_R2.fastq > $adapterTrimmed/${pre}.fastq.adapterTrim.metrics &
done
wait
# Takes output from cutadapt round 1. Run to trim off the 3’ adapters on read 2, to control for double ligation events
# add "--nextseq-trim=20" for our experiment
# --cores=21 for fastq. 
# Make also sure that you have pigz (parallel gzip) installed if you use multiple cores and write to a .gz output file. Otherwise, compression of the output will be done in a single thread and therefore be the main bottleneck.
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	/home1/04935/shaojf/.local/bin/cutadapt --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 --cores=21 -A a1=AACTTGTAGATCGGA -A a2=AGGACCAAGATCGGA -A a3=ACTTGTAGATCGGAA -A a4=GGACCAAGATCGGAA -A a5=CTTGTAGATCGGAAG -A a6=GACCAAGATCGGAAG -A a7=TTGTAGATCGGAAGA -A a8=ACCAAGATCGGAAGA -A a9=TGTAGATCGGAAGAG -A a10=CCAAGATCGGAAGAG -A a11=GTAGATCGGAAGAGC -A a12=CAAGATCGGAAGAGC -A a13=TAGATCGGAAGAGCG -A a14=AAGATCGGAAGAGCG -A a15=AGATCGGAAGAGCGT -A a16=GATCGGAAGAGCGTC -A a17=ATCGGAAGAGCGTCG -A a18=TCGGAAGAGCGTCGT -A a19=CGGAAGAGCGTCGTG -A a20=GGAAGAGCGTCGTGT -o $adapterTrimmed/${pre}_R1.adapterTrim.round2.fastq -p $adapterTrimmed/${pre}_R2.adapterTrim.round2.fastq $adapterTrimmed/${pre}_R1.adapterTrim.fastq $adapterTrimmed/${pre}_R2.adapterTrim.fastq > $adapterTrimmed/${pre}.fastq.adapterTrim.round2.metrics &
done
wait

## step 3, STAR rmRep
date
echo "step 3, STAR rmRep ......"
mkdir $rmrep
# STAR rmRep: Takes output from cutadapt round 2. Maps to human specific version of RepBase used to remove repetitive elements, helps control for spurious artifacts from rRNA (& other) repetitive reads.
# remove --outSAMunmapped Within --outStd BAM_Unsorted > $rmrep/${pre}.round2.rep.bam
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $repindexdir --genomeLoad LoadAndRemove \
	--readFilesIn $adapterTrimmed/${pre}_R1.adapterTrim.round2.fastq \
	$adapterTrimmed/${pre}_R2.adapterTrim.round2.fastq \
	--outFilterMultimapNmax 30 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix $rmrep/${pre}.round2.rep. \
	--outSAMattributes All --outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local
	# --alignEndsType EndToEnd

	samtools view $rmrep/${pre}.round2.rep.Aligned.out.bam | count_aligned_from_sam.py > $rmrep/${pre}.round2.rep.metrics
	fastqc --outdir $fastqc_dir $rmrep/${pre}.round2.rep.Unmapped.out.mate1 > $rmrep/${pre}.round2.rep.Unmapped.out.mate1.dummy_fastqc &
	fastqc --outdir $fastqc_dir $rmrep/${pre}.round2.rep.Unmapped.out.mate2 > $rmrep/${pre}.round2.rep.Unmapped.out.mate2.dummy_fastqc &
	fastq-sort --id $rmrep/${pre}.round2.rep.Unmapped.out.mate1 > $rmrep/${pre}.round2.rep.Unmapped.out.sorted.mate1 &
	fastq-sort --id $rmrep/${pre}.round2.rep.Unmapped.out.mate2 > $rmrep/${pre}.round2.rep.Unmapped.out.sorted.mate2 &
done
wait

## step 4, STAR genome
date
echo "step 4, STAR genome ......"
mkdir $mappingtogenome
# --outSAMunmapped Within \
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn $rmrep/${pre}.round2.rep.Unmapped.out.sorted.mate1 \
	$rmrep/${pre}.round2.rep.Unmapped.out.sorted.mate2 \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix $mappingtogenome/${pre}.round2.rmRep.genome. \
	--outSAMattributes All --outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:foo \
	--alignEndsType Local
	# --alignEndsType EndToEnd
done

## step 5, process STAR results
# Barcode_collapse_pe: takes output from STAR genome mapping. Custom random-mer-aware script for PCR duplicate removal
# WARNING: Genome (-g) files are ignored when BAM input is provided. 
# rm "--genome $genomesize "
date
echo "step 5, process STAR results ......"
mkdir $bw_dir
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	barcode_collapse_pe.py --bam $mappingtogenome/${pre}.round2.rmRep.genome.Aligned.out.bam --out_file $mappingtogenome/${pre}.round2.rmRep.genome.rmDup.bam --metrics_file $mappingtogenome/${pre}.round2.rmRep.genome.rmDup.metrics
	picard SortSam INPUT=$mappingtogenome/${pre}.round2.rmRep.genome.rmDup.bam OUTPUT=$mappingtogenome/${pre}.round2.rmRep.genome.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
	samtools index -@ 68 $mappingtogenome/${pre}.round2.rmRep.genome.sorted.bam
	samtools view -hb -f 128 -@ 68 $mappingtogenome/${pre}.round2.rmRep.genome.sorted.bam > $mappingtogenome/${pre}.round2.rmRep.genome.sorted.r2.bam
	####
	make_bigwig_files.py --bam $mappingtogenome/${pre}.round2.rmRep.genome.sorted.r2.bam --bw_pos $bw_dir/${pre}.round2.rmRep.genome.sorted.r2.norm.pos.bw --bw_neg $bw_dir/${pre}.round2.rmRep.genome.sorted.r2.norm.neg.bw
	####
done

## step 6, calling peaks with clipper
# Clipper: Takes results from samtools view. Calls peaks on those files
# rm "--bonferroni --superlocal --threshold-method binomial "
date
echo "step 6, calling peaks with clipper ......"
mkdir $clippeaks
for pre in MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	clipper -b $mappingtogenome/${pre}.round2.rmRep.genome.sorted.r2.bam -s hg19 -o $clippeaks/${pre}.r2.peaks.bed --save-pickle --processors 64 -p
	fix_scores.py --bed $clippeaks/${pre}.r2.peaks.bed --out_file $clippeaks/${pre}.r2.peaks.fixed.bed
	bedToBigBed $clippeaks/${pre}.r2.peaks.fixed.bed $genomesize $clippeaks/${pre}.r2.peaks.fixed.bb -type=bed6+4
done


###### test
hisat2index=/home1/04935/shaojf/scratch/hisat2_index/hg19/genome
trimadapter=/home1/04935/shaojf/myTools/Trimmomatic-0.36/trimmomatic-0.36.jar
trimadapterdir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/trimadapter_dir
barcodesfile=/home1/04935/shaojf/scratch/eCLIP_MCF-7/test.adapter.txt
mkdir $trimadapterdir
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	fwdinput=$fastq_dir/${pre}_R1.fastq
	revinput=$fastq_dir/${pre}_R2.fastq
	i=1
	while read line
	do
		echo $line | awk '{print ">"$1"\n"$2}' > $pre.test.fa
		java -jar $trimadapter PE -threads 68 -phred33 \
		$fwdinput $revinput\
		$trimadapterdir/${pre}_R1.t${i}.fastq $trimadapterdir/${pre}_R1.t${i}.up.fastq \
		$trimadapterdir/${pre}_R2.t${i}.fastq $trimadapterdir/${pre}_R2.t${i}.up.fastq \
		ILLUMINACLIP:$pre.test.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		fastqc -t 16 -o $fastqc_dir $trimadapterdir/${pre}_R1.t${i}.fastq &
		fastqc -t 16 -o $fastqc_dir $trimadapterdir/${pre}_R2.t${i}.fastq &
		fwdinput=$trimadapterdir/${pre}_R1.t${i}.fastq
		revinput=$trimadapterdir/${pre}_R2.t${i}.fastq
		i=$((i+1))
	done < $barcodesfile
	rm $pre.test.fa
done
wait

for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $repindexdir --genomeLoad LoadAndRemove \
	--readFilesIn $trimadapterdir/${pre}_R1.t4.fastq \
	$trimadapterdir/${pre}_R2.t4.fastq \
	--outFilterMultimapNmax 30 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix $rmrep/${pre}.t4.rep. \
	--outSAMattributes All --outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:$pre \
	--alignEndsType Local
	# --alignEndsType EndToEnd

	samtools view -@ 64 $rmrep/${pre}.t4.rep.Aligned.out.bam | count_aligned_from_sam.py > $rmrep/${pre}.t4.rep.metrics
	fastqc -t 16 --outdir $fastqc_dir $rmrep/${pre}.t4.rep.Unmapped.out.mate1 > $rmrep/${pre}.t4.rep.Unmapped.out.mate1.dummy_fastqc &
	fastqc -t 16 --outdir $fastqc_dir $rmrep/${pre}.t4.rep.Unmapped.out.mate2 > $rmrep/${pre}.t4.rep.Unmapped.out.mate2.dummy_fastqc &
	fastq-sort --id $rmrep/${pre}.t4.rep.Unmapped.out.mate1 > $rmrep/${pre}.t4.rep.Unmapped.out.sorted.mate1 &
	fastq-sort --id $rmrep/${pre}.t4.rep.Unmapped.out.mate2 > $rmrep/${pre}.t4.rep.Unmapped.out.sorted.mate2 &
done
wait

## step 4, STAR genome
date
echo "step 4, STAR genome ......"
mkdir $mappingtogenome
# --outSAMunmapped Within \
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	STAR --runMode alignReads --runThreadN 68 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn $rmrep/${pre}.t4.rep.Unmapped.out.sorted.mate1 \
	$rmrep/${pre}.t4.rep.Unmapped.out.sorted.mate2 \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix $mappingtogenome/${pre}.t4.rmRep.genome. \
	--outSAMattributes All --outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:$pre \
	--alignEndsType Local
	# --alignEndsType EndToEnd
done

# cat STAR_genome/*t4.rmRep.genome.Log.final.out | grep "Uniquely mapped reads %"
#                         Uniquely mapped reads % |	33.84%
#                         Uniquely mapped reads % |	13.09%
#                         Uniquely mapped reads % |	18.96%

for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	barcode_collapse_pe.py --bam $mappingtogenome/${pre}.t4.rmRep.genome.Aligned.out.bam --out_file $mappingtogenome/${pre}.t4.rmRep.genome.rmDup.bam --metrics_file $mappingtogenome/${pre}.t4.rmRep.genome.rmDup.metrics
	picard SortSam INPUT=$mappingtogenome/${pre}.t4.rmRep.genome.rmDup.bam OUTPUT=$mappingtogenome/${pre}.t4.rmRep.genome.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
	samtools index -@ 68 $mappingtogenome/${pre}.t4.rmRep.genome.sorted.bam
	samtools view -hb -f 128 -@ 68 $mappingtogenome/${pre}.t4.rmRep.genome.sorted.bam > $mappingtogenome/${pre}.t4.rmRep.genome.sorted.r2.bam
done

for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	####
	make_bigwig_files.py --bam $mappingtogenome/${pre}.t4.rmRep.genome.sorted.r2.bam --genome $genomesize --bw_pos $bw_dir/${pre}.t4.rmRep.genome.sorted.r2.norm.pos.bw --bw_neg $bw_dir/${pre}.t4.rmRep.genome.sorted.r2.norm.neg.bw &
	####
	clipper -b $mappingtogenome/${pre}.t4.rmRep.genome.sorted.r2.bam -s hg19 -o $clippeaks/${pre}.r2.peaks.bed --save-pickle # --processors=64 -p
	fix_scores.py --bed $clippeaks/${pre}.r2.peaks.bed --out_file $clippeaks/${pre}.r2.peaks.fixed.bed
	bedToBigBed $clippeaks/${pre}.r2.peaks.fixed.bed $genomesize $clippeaks/${pre}.r2.peaks.fixed.bb -type=bed6+4
done
wait

# for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
# do
# 	hisat2 --rg ID:$pre -p 68 -x $hisat2index \
# 	-1 $rmrep/${pre}.round2.rep.Unmapped.out.sorted.mate1 \
# 	-2 $rmrep/${pre}.round2.rep.Unmapped.out.sorted.mate2 \
# 	-S $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sam
# done

# for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
# do
# 	hisat2 --rg ID:$pre -p 68 -x $hisat2index \
# 	-1 $trimadapterdir/${pre}_R1.t4.fastq \
# 	-2 $trimadapterdir/${pre}_R2.t4.fastq \
# 	-S $mappingtogenome/${pre}.hisat2.sam
# done


for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	samtools view -b -o $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.bam -@ 68 $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sam
	barcode_collapse_pe.py --bam $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.bam --out_file $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.rmDup.bam --metrics_file $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.rmDup.metrics
	picard SortSam INPUT=$mappingtogenome/${pre}.round2.rmRep.genome.hisat2.rmDup.bam OUTPUT=$mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
	samtools index -@ 68 $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.bam
	samtools view -hb -f 128 -@ 68 $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.bam > $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.r2.bam
	####
	make_bigwig_files.py --bam $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.r2.bam --bw_pos $bw_dir/${pre}.round2.rmRep.genome.hisat2.sorted.r2.norm.pos.bw --bw_neg $bw_dir/${pre}.round2.rmRep.genome.hisat2.sorted.r2.norm.neg.bw
	####
done