#!/bin/sh
bedtools intersect -wa -u -a <(gunzip -c hg19.Snyder.K562.EP300.opt.ENCFF549TYR.bed.gz) -b <(gunzip -c hg19.K562.DNase.1.ENCFF121EHW.bed.gz) > hg19.Snyder.K562.EP300.opt.ENCFF549TYR.DNase1.txt
bedtools intersect -wa -u -a <(gunzip -c hg19.Snyder.K562.EP300.opt.ENCFF549TYR.bed.gz) -b <(gunzip -c hg19.K562.DNase.2.ENCFF355WGC.bed.gz) > hg19.Snyder.K562.EP300.opt.ENCFF549TYR.DNase2.txt
### 28288 -> 22329,17325

### use the merged p300 peaks from the too replicates
cat hg19.Snyder.K562.EP300.opt.ENCFF549TYR.DNase* | cut -f 1-3 | sort | uniq -c | awk -vOFS="\t" '{print $2,$3,$4,"p300_"NR,$1,"."}' > hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.bed
### 28288 -> 22708

### only keep peaks from Intergenic and intron
annotatePeaks.pl hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.bed hg19 > hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.homer.txt
awk -F"\t" -vOFS="\t" '$8=="Intergenic" && $10 < -3000 {print $2,$3,$4+2000,$1"_F",$6,"+\n"$2,$3-2000,$4,$1"_R",$6,"-"}' hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.homer.txt > candidate.Intergenic.bed
awk -F"\t" -vOFS="\t" '$8~/intron/ && $10 > 3000 {print $8,$2,$3,$4+2000,$1"_F",$6,"+\n"$8,$2,$3-2000,$4,$1"_R",$6,"-"}' hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.homer.txt | perl /home1/04935/shaojf/myScripts/add_any_2files_together.pl hg19.basic.annotation /dev/stdin 0 0 | awk -F"\t" '$7!=$12' | cut -f 2-7 > candidate.intron.bed

### count the gro-seq signal
for pre in candidate.Intergenic.bed candidate.intron.bed
do
	annotatePeaks.pl $pre hg19 -fpkm -strand + -d K562.GRO-seq/ K562.GRO-cap.ctl/ K562.GRO-cap/ 1> $pre.gro 2> $pre.gro.log &
done

for pre in candidate.Intergenic.bed candidate.intron.bed
do
	cut -f 1-6,20- $pre.gro > $pre.gro.sim
done

for pre in candidate.Intergenic.bed candidate.intron.bed
do
	cut -f 1-6,20- $pre.gro > $pre.gro.sim
done

for pre in candidate.Intergenic.bed candidate.intron.bed
do
	awk -F"\t" '$7>0.1{print $1}' $pre.gro.sim | sed 's/_F//;s/_R//' | sort | uniq -c | sort -k1,1nr > $pre.gro.sim.count
done

cat *gro.sim | grep -v "PeakID" | awk -vOFS="\t" '$7>0.1{print $2,$3,$4,$1,$7,$5}' | bedtools sort -i - > candidate.eRNA.withGRO.bed

### count the gro-seq signal for homer annotation, to do fisher test
cat coding.ann.txt cpgIsland.ann.txt exons.ann.txt intergenic.ann.txt introns.ann.txt miRNA.ann.txt ncRNA.ann.txt promoters.ann.txt protein-coding.ann.txt tts.ann.txt utr3.ann.txt utr5.ann.txt | awk -F"\t" -vOFS="\t" '{if($4-$3==0){print $1,$2,$3-1,$4,$5,$6}else{print $0}}' > hg19.homer.ann.txt
# pre=hg19.homer.ann.txt
# annotatePeaks.pl $pre hg19 -fpkm -strand + -d K562.GRO-seq/ K562.GRO-cap.ctl/ K562.GRO-cap/ 1> $pre.gro 2> $pre.gro.log &
# cut -f 1-5,8,20- hg19.homer.ann.txt.gro > hg19.homer.ann.txt.gro.sim
# awk -F"\t" '$7>0.1{print $1}'  hg19.homer.ann.txt.gro.sim
for pre in coding.ann.txt cpgIsland.ann.txt exons.ann.txt intergenic.ann.txt introns.ann.txt miRNA.ann.txt ncRNA.ann.txt promoters.ann.txt protein-coding.ann.txt tts.ann.txt utr3.ann.txt utr5.ann.txt
do
	annotatePeaks.pl <(awk -F"\t" -vOFS="\t" '{if($4-$3==0){print $1,$2,$3-1,$4,$5,$6}else{print $0}}' /work/04935/shaojf/stampede2/myTools/HOMER/data/genomes/hg19/annotations/basic/$pre) hg19 -fpkm -strand + -d K562.GRO-seq/ K562.GRO-cap.ctl/ K562.GRO-cap/ 1> $pre.gro 2> $pre.gro.log &
done
for pre in coding.ann.txt cpgIsland.ann.txt exons.ann.txt intergenic.ann.txt introns.ann.txt miRNA.ann.txt ncRNA.ann.txt promoters.ann.txt protein-coding.ann.txt tts.ann.txt utr3.ann.txt utr5.ann.txt
do
	awk -F"\t" '$20>0.1{print $1}' $pre.gro > $pre.gro.sim
done


#### just rna
pre=refseq
analyzeRepeats.pl rna hg19 -fpkm -strand + -d K562.GRO-seq/ K562.GRO-cap.ctl/ K562.GRO-cap/ 1> $pre.gro 2> $pre.gro.log &
cut -f 1-6,9- $pre.gro | awk -F"\t" '$7>0.1' > $pre.gro.sim
grep -v "Transcript" refseq.gro.sim | awk -vOFS="\t" '{a=$3;b=$3+2000;if($5=="-"){a=$4-2000;b=$4}{print $2,a,b,$1,$7,$5}}' > refseq.utr5.bed
grep -v "Transcript" refseq.gro.sim | awk -vOFS="\t" '{a=$3;b=$3+2000;if($5=="+"){a=$4-2000;b=$4}{print $2,a,b,$1,$7,$5}}' > refseq.utr3.bed
# #### for reference
# for pre in candidate.Intergenic.bed candidate.intron.bed
# do
# 	awk -vOFS="\t" '{print $4,$1,$2,$3,$6,$5}' $pre > a.$pre.peaks
# 	analyzeRepeats.pl a.$pre.peaks hg19 -count pausing -fpkm -strand + -d K562.GRO-seq/ K562.GRO-cap.ctl/ K562.GRO-cap/ 1> repeat.$pre.gro 2> repeat.$pre.gro.log &
# done

# ### add strand
# ###awk -vOFS="\t" '{print $1,$2,$3,"p300_F_"NR,$4,"+\n"$1,$2,$3,"p300_R_"NR,$4,"-"}' hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.bed > hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.withstrand.bed
# awk -vOFS="\t" '{print $1,$2,$3+2000,"p300_F_"NR,$4,"+\n"$1,$2-2000,$3,"p300_R_"NR,$4,"-"}' hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.bed > hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.withstrand.bed
# ### remove promoter and those in the same direction
# annotatePeaks.pl hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.withstrand.bed hg19 -strand + > hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.withstrand.homer.txt

# awk -F"\t" -vOFS="\t" '$8=="Intergenic" && $10 < -3000 {print $2,$3,$4,$1,$6,$5}' hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.withstrand.homer.txt > candidate.Intergenic.bed
# awk -F"\t" -vOFS="\t" '$8~/intron/ && $10 > 3000 {print $8,$2,$3,$4,$1,$6,$5}' hg19.Snyder.K562.EP300.opt.ENCFF549TYR.withDHSin1rep.withstrand.homer.txt | perl /home1/04935/shaojf/myScripts/add_any_2files_together.pl hg19.basic.annotation /dev/stdin 0 0 | awk -F"\t" '$7!=$12' | cut -f 2-7 > candidate.intron.bed

# bigWigToWig bigwigs.from.GEO/GSM1480325_K562_GROseq_minus.bigWig GSM1480325_K562_GROseq_minus.wig
# bigWigToWig bigwigs.from.GEO/GSM1480327_K562_PROseq_plus.bigWig GSM1480327_K562_PROseq_plus.wig
# bigWigToWig bigwigs.from.GEO/GSM1480327_K562_PROseq_plus.bw GSM1480327_K562_PROseq_plus.wig
# bigWigToWig bigwigs.from.GEO/GSM1480327_K562_PROseq_minus.bw GSM1480327_K562_PROseq_minus.wig
# for pre in candidate.Intergenic.bed candidate.intron.bed
# do
# 	for post in *.wig
# 	do
# 		annotatePeaks.pl $pre hg19 -wig $post -strand + 1> $pre.$post 2> $pre.$post.log &
# 	done
# done

############ for eRNA <-> eCLIP
awk -F"\t" '$2=="K562" && $7=="hg19"{print $1".bed.gz\t"$7"."$2"."$3"."$4}' sim.desc.tsv | while read line
do
	pre=`echo $line | awk '{print $1}'`
	post=`echo $line | awk '{print $2}'`
	mv $pre $post.$pre
done
for f in hg19.K562.*; do gunzip -c $f | head -1; done
zcat hg19.K562.* > hg19.K562.eCLIP.normpeaks.bed

####
annotatePeaks.pl hg19.K562.eCLIP.normpeaks.bed hg19 1> hg19.K562.eCLIP.normpeaks.annotation 2> hg19.K562.eCLIP.normpeaks.annotation.log
#### should force same strand
# bedtools intersect -wao -a hg19.K562.eCLIP.normpeaks.bed -b candidate.eRNA.withGRO.bed > eCLIP.eRNA.txt
# bedtools intersect -wao -a candidate.eRNA.withGRO.bed -b hg19.K562.eCLIP.normpeaks.bed > eRNA.eCLIP.txt

# awk '$7!="."' eRNA.eCLIP.txt | cut -f 1-4,10 | sed 's/_K562_rep02//;s/_K562_rep01//' | sort | uniq > eRNA.eCLIP.interaction
# cut -f 4 eRNA.eCLIP.interaction | sort | uniq -c | sort -k1,1nr > eRNA.eCLIP.interaction.count 

####### force same strand
for f in candidate.eRNA.withGRO.bed refseq.utr?.bed
do
	bedtools intersect -wo -s -a $f -b hg19.K562.eCLIP.normpeaks.bed > $f.eCLIP.txt &
done
##### for eRNA: 200783 -> 31996. decreased a lot
##### eRNA.eCLIP.txt doesn't consider about strand.
for f in candidate.eRNA.withGRO.bed refseq.utr?.bed
do
	awk '$7!="."'  $f.eCLIP.txt | cut -f 1-4,10 | sed 's/_K562_rep02//;s/_K562_rep01//' | sort | uniq > $f.eCLIP.interaction
	cut -f 4 $f.eCLIP.interaction | sort | uniq -c | sort -k1,1nr > $f.eCLIP.interaction.count 
done
for f in candidate.eRNA.withGRO.bed refseq.utr?.bed
do
	cut -f 5 $f.eCLIP.interaction | sort | uniq -c > eCLIP.$f.interaction.count 
done
paste eCLIP.candidate.eRNA.withGRO.bed.interaction.count eCLIP.refseq.utr3.bed.interaction.count eCLIP.refseq.utr5.bed.interaction.count | awk -vOFS="\t" '{print $2,$1,$3,$5}' > eCLIP.eRNA.utr3.utr5.count

  #   629 candidate.eRNA.withGRO.bed.eCLIP.interaction.count
  # 21457 refseq.utr3.bed.eCLIP.interaction.count
  # 33862 refseq.utr5.bed.eCLIP.interaction.count
  #     7531 candidate.eRNA.withGRO.bed
  #    35919 refseq.utr3.bed
  #    35919 refseq.utr5.bed



















### chr6:16,192,628-16,341,494

for f in hg19.K562.*
# for f in hg19.K562.YWHAG-human*.bed.gz hg19.K562.Z*.bed.gz
do
	gunzip -c $f | annotatePeaks.pl /dev/stdin hg19 1> $f.annotation 2> $f.log &
done

grep -A 15 "Finding Closest TSS" *.log | grep Intergenic | cut -f 4 > stats.Intergenic
grep -A 15 "Finding Closest TSS" *.log | grep Intron | cut -f 4 > stats.Intron
grep -A 15 "Finding Closest TSS" *.log | grep 3UTR | cut -f 4 > stats.3UTR
grep -A 15 "Finding Closest TSS" *.log | grep 5UTR | cut -f 4 > stats.5UTR
grep -A 15 "Finding Closest TSS" *.log | grep TTS | cut -f 4 > stats.TTS
grep -A 15 "Finding Closest TSS" *.log | grep Promoter | cut -f 4 > stats.Promoter
grep -A 15 "Finding Closest TSS" *.log | grep ncRNA | cut -f 4 > stats.ncRNA
grep -A 15 "Finding Closest TSS" *.log | grep Exon | cut -f 4 > stats.exon
ll hg19.K562.*.gz | awk '{print $9}' > stats.files
awk '$11!="."{print $4}' eCLIP.eRNA.txt | sort | uniq -c | awk '{print $2"\t"$1}' > eCLIP.eRNA.count
echo "RBP	TotalPeaks	Intergenic	Intron	3UTR	5UTR	TTS	Promoter	ncRNA	exon	eRNA" | sed "\s" "\t" > stats.all.txt
paste stats.files stats.Intergenic stats.Intron stats.3UTR stats.5UTR stats.TTS stats.Promoter stats.ncRNA stats.exon <(cut -f 2eCLIP.eRNA.count) >> stats.all.txt
# 18433,528700,50356,84547,62512,62512,23368
# wc -l ~/myTools/HOMER/data/genomes/hg19/annotations/basic/*.txt
   #     24 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/centromeres.ann.txt
   # 550757 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/coding.ann.txt
   #  28691 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/cpgIsland.ann.txt
   # 685660 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/exons.ann.txt
   #    457 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/gaps.ann.txt
   #  18433 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/intergenic.ann.txt
   # 528700 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/introns.ann.txt
   #   2918 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/miRNA.ann.txt
   #  23368 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/ncRNA.ann.txt
   #  62512 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/promoters.ann.txt
   # 651012 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/protein-coding.ann.txt
   #   7700 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/pseudo.ann.txt
   #      3 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/rRNA.ann.txt
   #      1 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/scRNA.ann.txt
   #    500 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/snoRNA.ann.txt
   #    156 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/snRNA.ann.txt
   #  62512 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/tts.ann.txt
   #  50356 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/utr3.ann.txt
   #  84547 /home1/04935/shaojf/myTools/HOMER/data/genomes/hg19/annotations/basic/utr5.ann.txt
