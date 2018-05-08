#!/bin/sh
# add #!/usr/bin/perl
# cp /home1/04935/shaojf/myTools/eCLIP/YeoLab.scripts/gscripts-1.1/perl_scripts  ~/.local/bin/
### input experiment name
inputfile=MCF7-CLIP-input
IPA=MCF7-CLIP-RBM39-IPA
IPB=MCF7-CLIP-RBM39-IPB

barcodesfile=/home1/04935/shaojf/scratch/eCLIP_MCF-7/all.barcodes.txt
adapterfile=/home1/04935/shaojf/scratch/eCLIP_MCF-7/adapter.to.cut.txt

### sub directories
fastq_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqs
fastq_demux=/home1/04935/shaojf/scratch/eCLIP_MCF-7/demux_paired_end_res
fastqc_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqc_res
trimadapter=/home1/04935/shaojf/myTools/Trimmomatic-0.36/trimmomatic-0.36.jar
trimadapterdir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqs_trimmed
rmrep=/home1/04935/shaojf/scratch/eCLIP_MCF-7/STAR_rmRep
mappingtogenome=/home1/04935/shaojf/scratch/eCLIP_MCF-7/STAR_genome
bw_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/bigwigs
genomesize=/home1/04935/shaojf/scratch/star_index/hg19.chrom.sizes
chrfile=/home1/04935/shaojf/scratch/star_index/hg19.onlychr.bed
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
fastqc --threads 64 --outdir $fastqc_dir $fastq_dir/*.fastq* &
## Demultiplexing
# must do it; the random-mer was appended to the read name for later usage. 
mkdir $fastq_demux
# pre=$inputfile
# ln -s $fastq_dir/${pre}_R1.fastq.gz $fastq_demux/${pre}_R1.fastq.gz
# ln -s $fastq_dir/${pre}_R1.fastq.gz $fastq_demux/${pre}_R2.fastq.gz
# NNNNNCCTATAT	X1A
# NNNNNTGCTATT	X1B
# NNNNNTATACTT	X2A
# NNNNNATCTTCT	X2B
# pre=$inputfile
# demux_paired_end.py --fastq_1 $fastq_dir/${pre}_R1.fastq.gz --fastq_2 $fastq_dir/${pre}_R2.fastq.gz \
# 	--out_file_1 ${pre}_R1.fastq.gz --out_file_2 ${pre}_R2.fastq.gz \
# 	--length 10 -m $fastq_demux/${pre}.metrics &
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
# ln -s $fastq_dir/${pre}_R1.fastq.gz $trimadapterdir/${pre}_R1.raw.fastq.gz
# ln -s $fastq_dir/${pre}_R1.fastq.gz $trimadapterdir/${pre}_R2.raw.fastq.gz
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
# WARNING: Genome (-g) files are ignored when BAM input is provided. 
# rm "--genome $genomesize "
date
echo "step 5, process STAR results ......"
mkdir $bw_dir
for pre in $inputfile ${IPA}_X1A ${IPA}_X1B ${IPB}_X2A ${IPB}_X2B
do
	barcode_collapse_pe.py --bam $mappingtogenome/${pre}.rmRep.genome.Aligned.out.bam --out_file $mappingtogenome/${pre}.rmRep.genome.rmDup.bam --metrics_file $mappingtogenome/${pre}.rmRep.genome.rmDup.metrics
	picard SortSam INPUT=$mappingtogenome/${pre}.rmRep.genome.rmDup.bam OUTPUT=$mappingtogenome/${pre}.rmRep.genome.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
done
samtools merge -@ 68 $mappingtogenome/${IPA}.rmRep.genome.sorted.bam $mappingtogenome/${IPA}_X1A.rmRep.genome.sorted.bam $mappingtogenome/${IPA}_X1B.rmRep.genome.sorted.bam
samtools merge -@ 68 $mappingtogenome/${IPB}.rmRep.genome.sorted.bam $mappingtogenome/${IPB}_X2A.rmRep.genome.sorted.bam $mappingtogenome/${IPB}_X2B.rmRep.genome.sorted.bam
for pre in $inputfile $IPA $IPB
do
	samtools index -@ 68 $mappingtogenome/${pre}.rmRep.genome.sorted.bam
	samtools view -hb -f 128 -@ 68 $mappingtogenome/${pre}.rmRep.genome.sorted.bam > $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam
	####
	# picard SortSam INPUT=$mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam OUTPUT=$mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
	samtools sort -@ 68 -o $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.bam $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam
	samtools index -@ 68 $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.bam
	make_bigwig_files.py --bam $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.bam --genome $genomesize --bw_pos $bw_dir/${pre}.rmRep.genome.sorted.r2.norm.pos.bw --bw_neg $bw_dir/${pre}.rmRep.genome.sorted.r2.norm.neg.bw
	bedSort $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.norm.pos.bg $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.norm.pos.srt.bg &
	bedSort $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.norm.neg.t.bg $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.norm.net.t.srt.bg &
	wait
	bedGraphToBigWig $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.norm.pos.srt.bg $genomesize $bw_dir/bedGraphToBigWig.${pre}.rmRep.genome.sorted.r2.norm.pos.bw &
	bedGraphToBigWig $mappingtogenome/${pre}.rmRep.genome.sorted.r2.srt.norm.neg.t.srt.bg $genomesize $bw_dir/bedGraphToBigWig.${pre}.rmRep.genome.sorted.r2.norm.neg.bw &
	### MCF7-CLIP-RBM39-IPB.rmRep.genome.sorted.r2.srt.norm.neg.t.bg is not case-sensitive sorted at line 1128908.  Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again.
	### /home1/04935/shaojf/anaconda2/bin/bedtools genomecov -ibam stdin -bg -strand - -split -g ../../star_index/hg19.chrom.sizes
	####
done
# for pre in $inputfile $IPA $IPB
# do
# 	samtools view -b -f 128 -F 16 $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam -@ 68 > $mappingtogenome/${pre}.fwd1.bam
# 	samtools index -@ 68 $mappingtogenome/${pre}.fwd1.bam
# 	bamCoverage -b $mappingtogenome/${pre}.fwd1.bam -o $bw_dir/${pre}.fwd.bigWig --normalizeUsingRPKM &
# 	samtools view -b -f 144 $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam -@ 68 > $mappingtogenome/${pre}.rev1.bam
# 	samtools index -@ 68 $mappingtogenome/${pre}.rev1.bam
# 	bamCoverage -b $mappingtogenome/${pre}.rev1.bam -o $bw_dir/${pre}.rev.bigWig --normalizeUsingRPKM &
	
# 	samtools view -@ 68 -L $chrfile -o $mappingtogenome/${pre}.rmRep.genome.sorted.r2.onlychr.bam $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam
# 	# samtools sort -@ 68 -o $mappingtogenome/${pre}.rmRep.genome.sorted.r2.onlychr.srt.bam $mappingtogenome/${pre}.rmRep.genome.sorted.r2.onlychr.bam
# 	# samtools index -@ 68 $mappingtogenome/${pre}.rmRep.genome.sorted.r2.onlychr.srt.bam
# 	picard SortSam INPUT=$mappingtogenome/${pre}.rmRep.genome.sorted.r2.onlychr.bam OUTPUT=$mappingtogenome/${pre}.rmRep.genome.sorted.r2.onlychr.srt.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
# 	make_bigwig_files.py --bam $mappingtogenome/${pre}.rmRep.genome.sorted.r2.onlychr.srt.bam --genome $genomesize --bw_pos $bw_dir/${pre}.rmRep.genome.sorted.r2.norm.pos.bw --bw_neg $bw_dir/${pre}.rmRep.genome.sorted.r2.norm.neg.bw &
# done
# wait

# pre=MCF7-CLIP-RBM39-IPB
# # include reads that are 2nd in a pair (128);
# # exclude reads that are mapped to the reverse strand (16)
# samtools view -b -f 128 -F 16 ${pre}.rmRep.genome.sorted.r2.srt.bam -@ 68 > a.fwd1.bam

# # exclude reads that are mapped to the reverse strand (16) and
# # first in a pair (64): 64 + 16 = 80
# samtools view -b -f 80 ${pre}.rmRep.genome.sorted.r2.srt.bam -@ 68 > a.fwd2.bam

# # combine the temporary files
# samtools merge -@ 68 -f a.fwd.bam a.fwd1.bam a.fwd2.bam

# # index the filtered BAM file
# samtools index -@ 68 a.fwd.bam

# # run bamCoverage
# bamCoverage -b a.fwd.bam -o a.fwd.bigWig

# # remove the temporary files
# rm a.fwd*.bam

# # include reads that map to the reverse strand (128)
# # and are second in a pair (16): 128 + 16 = 144
# $ samtools view -b -f 144 ${pre}.rmRep.genome.sorted.r2.srt.bam > a.rev1.bam

# # include reads that are first in a pair (64), but
# # exclude those ones that map to the reverse strand (16)
# $ samtools view -b -f 64 -F 16 a.bam > a.rev2.bam

# # merge the temporary files
# $ samtools merge -f rev.bam rev1.bam rev2.bam

# # index the merged, filtered BAM file
# $ samtools index rev.bam

# # run bamCoverage
# $ bamCoverage -b rev.bam -o a.rev.bw

# # remove temporary files
# $ rm a.rev*.bam

###### use only chr as genome?
# login1.stampede2(1040)$ sed -n '1128907,1128909p' MCF7-CLIP-RBM39-IPB.rmRep.genome.sorted.r2.srt.norm.neg.t.bg
# chr17_ctg5_hap1	1493545	1493594	-0.118322933339
# chr17	30373	30389	-0.118322933339
# chr17	34755	34774	-0.118322933339
# login1.stampede2(1041)$ sed -n '1128906,1128909p' MCF7-CLIP-RBM39-IPB.rmRep.genome.sorted.r2.srt.norm.neg.t.bg
# chr17_ctg5_hap1	1482764	1482817	-0.118322933339
# chr17_ctg5_hap1	1493545	1493594	-0.118322933339
# chr17	30373	30389	-0.118322933339
# chr17	34755	34774	-0.118322933339


## step 6, calling peaks with clipper
# Clipper: Takes results from samtools view. Calls peaks on those files
# rm "--bonferroni --superlocal --threshold-method binomial "
date
echo "step 6, calling peaks with clipper ......"
mkdir $clippeaks
for pre in $inputfile $IPA $IPB
do
	clipper -b $mappingtogenome/${pre}.rmRep.genome.sorted.r2.bam -s hg19 -o $clippeaks/${pre}.r2.peaks.bed --save-pickle
	fix_scores.py --bed $clippeaks/${pre}.r2.peaks.bed --out_file $clippeaks/${pre}.r2.peaks.fixed.bed
	bedToBigBed $clippeaks/${pre}.r2.peaks.fixed.bed $genomesize $clippeaks/${pre}.r2.peaks.fixed.bb -type=bed6+4
done

## step 7, peak normalization
date
echo "step 7, peak normalization ......"
mkdir $normalizedpeaks
# uID	RBP	Cell line	CLIP_rep1	CLIP_rep2	INPUT
# 001	RBM39	MCF-7	$mappingtogenome/${IPA}.rmRep.genome.sorted.r2.bam	$mappingtogenome/${IPA}.rmRep.genome.sorted.r2.bam	$mappingtogenome/$inputfile.rmRep.genome.sorted.bam
echo | awk -vOFS="\t" '{print "uID","RBP","Cell line","CLIP_rep1","CLIP_rep2","INPUT"}' > $IPA.$IPB.manifest.txt
echo | awk -vOFS="\t" '{print "001","RBM39","MCF-7","ABCD1","ABCD2","ABCD3"}' >> $IPA.$IPB.manifest.txt
sed -i "s?ABCD1?$mappingtogenome/${IPA}.rmRep.genome.sorted.r2.bam?;s?ABCD2?$mappingtogenome/${IPB}.rmRep.genome.sorted.r2.bam?;s?ABCD3?$mappingtogenome/$inputfile.rmRep.genome.sorted.bam?" $IPA.$IPB.manifest.txt
Peak_input_normalization_wrapper.pl $IPA.$IPB.manifest.txt $normalizedpeaks

# AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
# demux_paired_end.py --fastq_1 MCF7-CLIP-RBM39-IPB_R1.fastq --fastq_2 MCF7-CLIP-RBM39-IPB_R2.fastq -b all.barcodes.txt --out_file_1 MCF7-CLIP-RBM39-IPB_R1.fastq.gz --out_file_2 MCF7-CLIP-RBM39-IPB_R2.fastq.gz --length 5 -m MCF7-CLIP-RBM39-IPB.metrics &
# demux_paired_end.py --fastq_1 MCF7-CLIP-RBM39-IPA_R1.fastq --fastq_2 MCF7-CLIP-RBM39-IPA_R2.fastq -b all.barcodes.txt --out_file_1 MCF7-CLIP-RBM39-IPA_R1.fastq.gz --out_file_2 MCF7-CLIP-RBM39-IPA_R2.fastq.gz --length 5 -m MCF7-CLIP-RBM39-IPA.metrics &

## step 2, trimmomatic
date
echo "step 2, trimmomatic ......"
mkdir $trimadapterdir
for pre in $inputfile $IPA $IPB
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
		ILLTMINACLIP:$pre.test.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		fastqc -t 16 -o $fastqc_dir $trimadapterdir/${pre}_R1.t${i}.fastq &
		fastqc -t 16 -o $fastqc_dir $trimadapterdir/${pre}_R2.t${i}.fastq &
		fwdinput=$trimadapterdir/${pre}_R1.t${i}.fastq
		revinput=$trimadapterdir/${pre}_R2.t${i}.fastq
		i=$((i+1))
	done < $adapterfile
	rm $pre.test.fa
done
wait



## step 2, Cutadapt
date
echo "step 2, Cutadapt ......"
mkdir $adapterTrimmed
# Cutadapt round 1: Takes output from demultiplexed files. Run to trim off both 5’ and 3’ adapters on both reads
# with cutadapt v1.15, not use "-f fastq --match-read-wildcards"
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 --nextseq-trim=20 --cores=21 -A a0=NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g g1=CTTCCGATCTACAAGTT -g g2=CTTCCGATCTTGGTCCT -A a1=AACTTGTAGATCGGA -A a2=AGGACCAAGATCGGA -A a3=ACTTGTAGATCGGAA -A a4=GGACCAAGATCGGAA -A a5=CTTGTAGATCGGAAG -A a6=GACCAAGATCGGAAG -A a7=TTGTAGATCGGAAGA -A a8=ACCAAGATCGGAAGA -A a9=TGTAGATCGGAAGAG -A a10=CCAAGATCGGAAGAG -A a11=GTAGATCGGAAGAGC -A a12=CAAGATCGGAAGAGC -A a13=TAGATCGGAAGAGCG -A a14=AAGATCGGAAGAGCG -A a15=AGATCGGAAGAGCGT -A a16=GATCGGAAGAGCGTC -A a17=ATCGGAAGAGCGTCG -A a18=TCGGAAGAGCGTCGT -A a19=CGGAAGAGCGTCGTG -A a20=GGAAGAGCGTCGTGT -o $adapterTrimmed/${pre}_R1.adapterTrim.fastq -p $adapterTrimmed/${pre}_R2.adapterTrim.fastq $fastq_dir/${pre}_R1.fastq $fastq_dir/${pre}_R2.fastq > $adapterTrimmed/${pre}.fastq.adapterTrim.metrics &
done
wait
# Takes output from cutadapt round 1. Run to trim off the 3’ adapters on read 2, to control for double ligation events
# add "--nextseq-trim=20" for our experiment
# --cores=21 for fastq. 
# Make also sure that you have pigz (parallel gzip) installed if you use multiple cores and write to a .gz output file. Otherwise, compression of the output will be done in a single thread and therefore be the main bottleneck.
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 --cores=21 -A a1=AACTTGTAGATCGGA -A a2=AGGACCAAGATCGGA -A a3=ACTTGTAGATCGGAA -A a4=GGACCAAGATCGGAA -A a5=CTTGTAGATCGGAAG -A a6=GACCAAGATCGGAAG -A a7=TTGTAGATCGGAAGA -A a8=ACCAAGATCGGAAGA -A a9=TGTAGATCGGAAGAG -A a10=CCAAGATCGGAAGAG -A a11=GTAGATCGGAAGAGC -A a12=CAAGATCGGAAGAGC -A a13=TAGATCGGAAGAGCG -A a14=AAGATCGGAAGAGCG -A a15=AGATCGGAAGAGCGT -A a16=GATCGGAAGAGCGTC -A a17=ATCGGAAGAGCGTCG -A a18=TCGGAAGAGCGTCGT -A a19=CGGAAGAGCGTCGTG -A a20=GGAAGAGCGTCGTGT -o $adapterTrimmed/${pre}_R1.adapterTrim.round2.fastq -p $adapterTrimmed/${pre}_R2.adapterTrim.round2.fastq $adapterTrimmed/${pre}_R1.adapterTrim.fastq $adapterTrimmed/${pre}_R2.adapterTrim.fastq > $adapterTrimmed/${pre}.fastq.adapterTrim.round2.metrics &
done
wait



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
		ILLTMINACLIP:$pre.test.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
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
	--outSAMattributes All --outSAMtype BAM Tnsorted \
	--outFilterType BySJout --outReadsTnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:$pre \
	--alignEndsType Local
	# --alignEndsType EndToEnd

	samtools view -@ 64 $rmrep/${pre}.t4.rep.Aligned.out.bam | count_aligned_from_sam.py > $rmrep/${pre}.t4.rep.metrics
	fastqc -t 16 --outdir $fastqc_dir $rmrep/${pre}.t4.rep.Tnmapped.out.mate1 > $rmrep/${pre}.t4.rep.Tnmapped.out.mate1.dummy_fastqc &
	fastqc -t 16 --outdir $fastqc_dir $rmrep/${pre}.t4.rep.Tnmapped.out.mate2 > $rmrep/${pre}.t4.rep.Tnmapped.out.mate2.dummy_fastqc &
	fastq-sort --id $rmrep/${pre}.t4.rep.Tnmapped.out.mate1 > $rmrep/${pre}.t4.rep.Tnmapped.out.sorted.mate1 &
	fastq-sort --id $rmrep/${pre}.t4.rep.Tnmapped.out.mate2 > $rmrep/${pre}.t4.rep.Tnmapped.out.sorted.mate2 &
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
	--readFilesIn $rmrep/${pre}.t4.rep.Tnmapped.out.sorted.mate1 \
	$rmrep/${pre}.t4.rep.Tnmapped.out.sorted.mate2 \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFileNamePrefix $mappingtogenome/${pre}.t4.rmRep.genome. \
	--outSAMattributes All --outSAMtype BAM Tnsorted \
	--outFilterType BySJout --outReadsTnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:$pre \
	--alignEndsType Local
	# --alignEndsType EndToEnd
done

# cat STAR_genome/*t4.rmRep.genome.Log.final.out | grep "Tniquely mapped reads %"
#                         Tniquely mapped reads % |	33.84%
#                         Tniquely mapped reads % |	13.09%
#                         Tniquely mapped reads % |	18.96%

for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	barcode_collapse_pe.py --bam $mappingtogenome/${pre}.t4.rmRep.genome.Aligned.out.bam --out_file $mappingtogenome/${pre}.t4.rmRep.genome.rmDup.bam --metrics_file $mappingtogenome/${pre}.t4.rmRep.genome.rmDup.metrics
	picard SortSam INPTT=$mappingtogenome/${pre}.t4.rmRep.genome.rmDup.bam OTTPTT=$mappingtogenome/${pre}.t4.rmRep.genome.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
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
# 	-1 $rmrep/${pre}.round2.rep.Tnmapped.out.sorted.mate1 \
# 	-2 $rmrep/${pre}.round2.rep.Tnmapped.out.sorted.mate2 \
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
	picard SortSam INPUT=$mappingtogenome/${pre}.round2.rmRep.genome.hisat2.rmDup.bam OTTPTT=$mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.bam VALIDATION_STRINGENCY=SILENT SO=coordinate CREATE_INDEX=true
	samtools index -@ 68 $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.bam
	samtools view -hb -f 128 -@ 68 $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.bam > $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.r2.bam
	####
	make_bigwig_files.py --bam $mappingtogenome/${pre}.round2.rmRep.genome.hisat2.sorted.r2.bam --bw_pos $bw_dir/${pre}.round2.rmRep.genome.hisat2.sorted.r2.norm.pos.bw --bw_neg $bw_dir/${pre}.round2.rmRep.genome.hisat2.sorted.r2.norm.neg.bw
	####
done