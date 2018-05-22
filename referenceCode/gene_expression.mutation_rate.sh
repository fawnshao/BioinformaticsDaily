#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/TIP60.KAT5/ mygoogle:TIP60.KAT5/
rboxplot=/home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/myboxplot.R
rquantile=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/data.quantile.R
myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl
gene=KAT5
grep KAT5 ~/stampede2/refs/hg19.refGene.tss.uniq.srt.bed > KAT5.hg19.refGene.tss
grep KAT5 <(gunzip -c ~/stampede2/refs/GENCODE/gencode.v19.annotation.gtf.gz) > KAT5.gencode.v19.annotation.gtf
# gunzip -c ../../mutations/ICGC.release.27/specimen.all_projects.tsv.gz | cut -f 1-2,7 > specimen.info

cut -f 4 KAT5.hg19.refGene.tss | tr ";" "\n" | tr "|" "\n" | sort | uniq > KAT5.ids
echo "ENSG00000172977" >> KAT5.ids
awk '{print $12}' KAT5.gencode.v19.annotation.gtf | sed 's/"//g;s/;//' | sort | uniq | grep ENST | cut -f 1 -d"." >> KAT5.ids

###########
# cd exprs
for f in exp_*.tsv.gz
do
	gunzip -c $f | cut -f 2-3,8-9 | grep -wf KAT5.ids > KAT5.$f.txt
done

# for f in `ll KAT5.exp_* | awk '$5==0{print $9}' | sed 's/KAT5.//;s/.txt//'`
# do
# 	gunzip -c $f | cut -f 1-2,8-9 | grep -wf KAT5.ids > KAT5.$f.txt
# done

for f in `ll KAT5.exp_*.tsv.gz.txt | awk '$5 > 0 {print $9}'`
do
	Rscript $rboxplot $f
done

for f in `ll KAT5.exp_*.tsv.gz.txt | awk '$5 > 0 {print $9}'`
do
	perl $myperl specimen.info $f 0 1 | cut -f 1-4,7 > $f.specimen
done

for f in KAT5.exp_*.tsv.gz.txt.specimen
do
	Rscript $rquantile $f
done

###########
# cd muts
for f in simple_somatic_mutation.open.*.tsv.gz.srt.bed
do
	pre=`echo $f | sed 's/simple_somatic_mutation.open.//;s/.tsv.gz.srt.bed//'`
	annotatePeaks.pl <(awk '{print "chr"$0}' $f) hg19 1> $pre.homer 2>$pre.homer.log &
done

# grep "Done annotating peaks file" *.log
# for pre in `grep "terminate" *.log | cut -f 1 -d "."`
# grep "Done annotating peaks file" *.log | cut -d"." -f 1 > a
# for pre in `ls *.log | cut -f 1 -d"." | grep -vwf a`
# for pre in `grep "Out of memory" *.log | cut -d"." -f1`
# do
# 	f=simple_somatic_mutation.open.$pre.tsv.gz.srt.bed
# 	annotatePeaks.pl <(awk '{print "chr"$0}' $f) hg19 1> $pre.homer 2>$pre.homer.log &
# done

# grep "Out of memory" *.log
# BRCA-FR.homer.log:Out of memory!
# ESAD-UK.homer.log:Out of memory!
# LINC-JP.homer.log:Out of memory!
# MELA-AU.homer.log:Out of memory!
# PAEN-AU.homer.log:Out of memory!
# SKCA-BR.homer.log:Out of memory!

echo "cancers regions counts" | tr "\s" "\t" > ICGC.r27.mut.stats
for f in *.homer.log
do
	pre=`echo $f | sed 's/.homer.log//'`
	grep -A 11 "3UTR" $f | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%d\n", var,$1,$2}' >> ICGC.r27.mut.stats
done

echo "cancers regions counts class" | tr "\s" "\t" > ICGC.r27.mut.byexpr.stats
for f in *.homer
do
	pre=`echo $f | sed 's/.homer//'`
	cut -f 8 $pre.homer | cut -d "(" -f 1 | tail -n +2 | sort | uniq -c | sed 's/ UTR/UTR/' | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%d\n", var,$2,$1,"total"}' >> ICGC.r27.mut.byexpr.stats
	# cut -f 8 $pre.homer | cut -d "(" -f 1 | tail -n +2 | sort | uniq -c | sed 's/ UTR/UTR/' | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%d\n", var,$2,$1,"total"}' >> ICGC.r27.mut.byexpr.stats
done

###########
# cd /home1/04935/shaojf/stampede2/TIP60.KAT5/association_expr_mut
# non-coding is exon for non-coding RNA
for f in *.homer
do
	pre=`echo $f | sed 's/.homer//'`
	# if [ ! -f KAT5.exp_array.$pre.tsv.gz.txt.specimen -a ! -f KAT5.exp_seq.$pre.tsv.gz.txt.specimen ]
	if [ ! -f KAT5.exp_array.$pre.tsv.gz.txt.specimen ] && [ ! -f KAT5.exp_seq.$pre.tsv.gz.txt.specimen ]
	then
		rm $pre
	fi
done

# in homer annotation, NA is from chrMT
echo "cancers specimens regions counts" | tr " " "\t" > WGS.mut.byspeciman.stats
# for f in `ls *.homer | grep -v "MALY-DE" | grep -v "PBCA-US"`
for f in *.homer
do
	pre=`echo $f | sed 's/.homer//'`
	awk -F"\t" '$1~/\|WGS/{print $1"\t"$8}' $f | sed 's/ UTR/UTR/;s/|/\t/' | awk -F"\t" '{print $1,$3}' | awk '{print $1"\t"$2}' | sort | uniq -c | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%s\t%d\n", var,$3,$2,$1}' >> WGS.mut.byspeciman.stats
done

for f in KAT5.exp_*.tsv.gz.txt.specimen
do
	pre=`echo $f | sed 's/.tsv.gz.txt.specimen//'`
	perl $myperl WGS.mut.byspeciman.stats $f 2 1 | grep -v "/" | cut -f 1-5,7,9 > WGS.mut.$pre.tsv
done


