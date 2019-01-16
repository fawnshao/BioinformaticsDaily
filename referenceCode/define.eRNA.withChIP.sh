#!/bin/bash
###################################################
# use ENCODE Hela-S3 H3K27ac and EP300 peaks to define eRNA
# https://www.nature.com/articles/nature14906
# Peak analysis of RNAPII ChIP-seq data after EGF stimulation was performed using HOMER 4.6 (run in ‘factor’ mode). Next, we used the BEDtools suite to discard any peak overlapping to: (i) all exons from Hg19 UCSC Known Genes (with additional 2 kb surrounding every exon); (ii) RNA Genes (from the Hg18 genome annotation table, plus additional 1 kb); (iii) tRNA Genes (Hg19, plus additional 1 kb). We further selected peaks overlapping (±400 bp) with H3K27ac peaks from the ENCODE ChIP-seq obtained in HeLa-S3 (GEO GSE31477). The analysis resulted in 2,029 regions that were further examined for their transcriptional response to EGF. Briefly, we centred a 6-kb window at the midst of the RNAPII peak and we used HOMER 4.6 to calculate RPKM across the entire eRNA locus using chromRNA-seq data before and after EGF induction. We selected a group of 225 EGF-inducible eRNAs displaying a fold change greater than 2 (ctrl versus EGF) and identified the nearest EGF regulated gene (fold change RPKM >1.6). 91 EGF-induced enhancer RNAs located within 500 kb from the nearest EGF-responsive protein-coding genes were selected for further analysis.
# gunzip -c hg19.RefSeq.all.gtf.gz | awk '$3=="exon"' | awk -vOFS="\t" '{print $1,$4,$5,$10,$8,$7}' | sed 's/;//g;s/"//g;s/\.[0-9]*//' > hg19.RefSeq.exon.bed
# gunzip -c hg19.RefSeq.bed.gz | sed 's/\.[0-9]*//' | cut -f 1-6 > hg19.RefSeq.gene.bed
# gunzip -c hg19.lincRNA.gtf.gz | awk '$3=="exon"' | awk -vOFS="\t" '{print $1,$4,$5,$10,$8,$7}' | sed 's/;//g;s/"//g' > hg19.lincRNA.exon.bed
# gunzip -c hg19.lincRNA.bed.gz | cut -f 1-6 > hg19.lincRNA.gene.bed
# gunzip -c hg19.tRNA.gtf.gz | awk -vOFS="\t" '{print $1,$4,$5,$10,$8,$7}' | sed 's/;//g;s/"//g;' > hg19.tRNA.exon.bed

genebed=/home1/04935/shaojf/stampede2/refs/RefSeq/GRCh37/hg19.RefSeq.gene.bed
exonbed=/home1/04935/shaojf/stampede2/refs/RefSeq/GRCh37/hg19.RefSeq.exon.bed
lincbed=/home1/04935/shaojf/stampede2/refs/RefSeq/GRCh37/hg19.lincRNA.gene.bed
trnabed=/home1/04935/shaojf/stampede2/refs/RefSeq/GRCh37/hg19.tRNA.exon.bed

p300peaks=hg19.MichaelSnyder.p300.rep12opt.ENCFF784YVX.bed
k27acpeaks=GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak
pol2peaks=hg19.RichardMyers.pol2.rep12opt.ENCFF417LLT.bed

bedtools intersect -wa -u -a $p300peaks -b <(awk -vOFS="\t" '{a=$2-400;if(a<0){a=0}print $1,a,$3+400,$4,$5,$6}' $k27acpeaks) | bedtools intersect -v -a - -b <(awk -vOFS="\t" '{a=$2-2000;if(a<0){a=0}print $1,a,$3+2000,$4,$5,$6}' $exonbed) | bedtools intersect -v -a - -b <(awk -vOFS="\t" '{a=$2-3000;if(a<0){a=0}print $1,a,$2,$4,$5,$6"\n"$1,$3,$3+3000,$4,$5,$6}' $genebed) | bedtools intersect -v -a - -b <(awk -vOFS="\t" '{a=$2-1000;if(a<0){a=0}print $1,a,$3+1000,$4,$5,$6}' $lincbed) | bedtools intersect -v -a - -b <(awk -vOFS="\t" '{a=$2-1000;if(a<0){a=0}print $1,a,$3+1000,$4,$5,$6}' $trnabed) > enhancer.p300.peaks

bedtools intersect -v -a enhancer.p300.peaks -b $genebed > intergenic.enhancer.p300.peaks
bedtools intersect -wo -a enhancer.p300.peaks -b $genebed > intragenic.enhancer.p300.peaks

awk -vOFS="\t" '{mid=int(($2+$3)/2);print $1,mid-2000,mid,"-\n"$1,mid,mid+2000,"+"}' intergenic.enhancer.p300.peaks | bedtools sort -i - | uniq | awk -vOFS="\t" '{print $1,$2,$3,"intergenic_"NR,".",$4}' > intergenic.eRNA.bed

bedtools intersect -wo -a <(cut -f 1-3,16 intragenic.enhancer.p300.peaks | sort | uniq | cut -f 1-3 | uniq -c | awk -vOFS="\t" '$1==1{print $2,$3,$4}') -b $genebed | cut -f 1-3,9 | sort | uniq | awk -vOFS="\t" '{mid=int(($2+$3)/2);if($4=="-"){print $1,mid,mid+2000,"intragenic_"NR,".","+"}if($4=="+"){print $1,mid-2000,mid,"intragenic_"NR,".","-"}}' | bedtools sort -i - > intragenic.eRNA.bed
