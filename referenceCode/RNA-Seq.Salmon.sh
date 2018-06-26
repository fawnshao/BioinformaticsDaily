#!/bin/bash
# wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz  
# wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh37.75.ncrna.fa.gz
# gunzip -c Homo_sapiens.GRCh37.75.cdna.all.fa.gz Homo_sapiens.GRCh37.75.ncrna.fa.gz > Homo_sapiens.GRCh37.75.cdna.ncrna.fa
# cat Homo_sapiens.GRCh37.75.cdna.ncrna.fa ERCC92.fa > Homo_sapiens.GRCh37.75.cdna.ncrna.ERCC92.fa
# salmon index -p 68 -t Homo_sapiens.GRCh37.75.cdna.ncrna.ERCC92.fa -i Homo_sapiens.GRCh37.75_quasi_index 

# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz
# zcat gencode.v19.pc_transcripts.fa.gz gencode.v19.lncRNA_transcripts.fa.gz | cat - ERCC92.fa > gencode.v19.ERCC.fa
# salmon index -p 68 --gencode -t gencode.v19.ERCC.fa -i gencode.v19.ERCC_quasi_index 
# awk '$3=="transcript"{print $12"\t"$10}' gencode.v19.annotation.gtf | sed 's/"//g;s/;//g' | awk 'BEGIN{print "TXNAME\tGENEID";}{print $0}' > tx2gene.gencode.v19.tsv

salmonIndex=/home1/04935/shaojf/scratch/Salmon.index/gencode.v19.ERCC_quasi_index
gtf=/home1/04935/shaojf/scratch/Salmon.index/gencode.v19.annotation.gtf
workDir=/home1/04935/shaojf/scratch/YY1.mut.RNAseq
fastqDir=$workDir/fastqDir
salmonDIR=$workDir/salmon.res
mkdir -p $salmonDIR
cd $workDir

# salmon quant --help-reads
# --gcBias  
for f in $fastqDir/*_1.fastq
do
	pre=`echo $f | awk -F"/" '{print $NF}' | sed 's/_1.fastq//'`
	salmon quant -i $salmonIndex -l A \
	-1 $fastqDir/${pre}_1.fastq -2 $fastqDir/${pre}_2.fastq \
	-p 68 -o $salmonDIR/$pre.salmon # -g $gtf
done

echo "samples files" > $workDir/tximport.samples.txt
for f in $fastqDir/*_1.fastq
do
	pre=`echo $f | awk -F"/" '{print $NF}' | sed 's/_1.fastq//'`
	echo "$pre $salmonDIR/$pre.salmon/quant.sf" >> $workDir/tximport.samples.txt
done

Rscript 
