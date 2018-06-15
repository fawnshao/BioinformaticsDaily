#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/TIP60.KAT5/ mygoogle:TIP60.KAT5/
# how about CNV
rboxplot=/home1/04935/shaojf/myTools/BioinformaticsDaily/RVisualization/myboxplot.R
rquantile=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/data.quantile.R
rscatterplot=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/exp.mut.scatter.plot.R
rquantileboxplot=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/exp.mut.boxplot.R
rquantileviolin=/home1/04935/shaojf/myTools/BioinformaticsDaily/gene_expression.mutation_rate.R/exp.mut.violinplot.R
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

# /bin/sbatch -A Functional-and-struc run.sh 
# split -l 200000 simple_somatic_mutation.open.MELA-AU.tsv.gz.srt.bed
# for f in x*
# do
# 	annotatePeaks.pl <(awk '{print "chr"$0}' $f) hg19 1> $f.homer 2>$f.homer.log &
# done


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
	awk -F"\t" '$1~/\|WGS/{print $1"\t"$8}' $f | sed 's/ UTR/UTR/;s/|/\t/' | awk -F"\t" '{print $1,$3}' | awk '{print $1"\t"$2}' | sort | uniq -c | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%s\t%d\n", var,$2,$3,$1}' >> WGS.mut.byspeciman.stats
done

for f in KAT5.exp_*.tsv.gz.txt.specimen
do
	pre=`echo $f | sed 's/.tsv.gz.txt.specimen//'`
	perl $myperl WGS.mut.byspeciman.stats $f 1 1 | grep -v "/" | cut -f 1-5,8,9 > WGS.mut.$pre.tsv
done

for f in `ll WGS.mut.KAT5.exp_* | awk '$5 > 100{print $9}'`
do
	Rscript $rscatterplot $f
done

# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen.tumors.*.tsv .
# rm quantile.exp.mut.homer.KAT5.exp_*
# rm *tsv.gz.txt.specimen.tumors.*er.tsv
# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv .
for f in KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv
do
	pre=`echo $f | cut -f 1-3 -d "."`
	class=`echo $f | awk -F"." '{print $(NF-1)}'`
	perl $myperl WGS.mut.byspeciman.stats <(sed 's/"//g' $f) 1 1 | grep -v "/" | cut -f 1-5,8,9 | awk -v var=$class '{print $0"\t"var}' >> quantile.exp.mut.homer.$pre.tsv
done

for f in `wc -l quantile.exp.mut.homer.KAT5.exp_*.tsv | awk '$1 > 100 && $2!="total"{print $2}'`
do
	Rscript $rquantileboxplot $f
done

for f in `wc -l quantile.exp.mut.homer.KAT5.exp_*.tsv | awk '$1 > 100 && $2!="total"{print $2}'`
do
	Rscript $rquantileviolin $f
done

###########
# cd /home1/04935/shaojf/stampede2/TIP60.KAT5/sv
# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen .
echo "cancers specimens svtypes seqtypes counts" | tr " " "\t" > sv.byspeciman.stats
for f in structural_somatic_mutation.*.tsv.gz
do
	pre=`echo $f | sed 's/structural_somatic_mutation.//;s/.tsv.gz//'`
	if [ ! -f KAT5.exp_array.$pre.tsv.gz.txt.specimen ] && [ ! -f KAT5.exp_seq.$pre.tsv.gz.txt.specimen ]
	then
		rm $f
	else
		gunzip -c $f | cut -f 2-3,7,23 | tail -n +2 | sed 's/ /./g' | sort | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$5,$1}' >> sv.byspeciman.stats
	fi
done

for f in KAT5.exp_*.tsv.gz.txt.specimen
do
	pre=`echo $f | sed 's/.tsv.gz.txt.specimen//'`
	perl $myperl sv.byspeciman.stats $f 1 1 | grep -v "/" | cut -f 1-5,8,10 > exp.sv.$pre.tsv
done

for f in `ll exp.sv.KAT5.exp_* | awk '$5 > 100{print $9}'`
do
	Rscript $rscatterplot $f
done

for pre in `ll exp.sv.KAT5.exp_*.tsv | awk '$5 > 100{print $9}' | cut -f 5 -d"."`
do
	f=structural_somatic_mutation.$pre.tsv.gz
	annotatePeaks.pl <(gunzip -c $f | awk -F"\t" -vOFS="\t" '{print "chr"$12,$13-1,$13,$3"|"$2"|"$7"|"$22"\nchr"$17,$18-1,$18,$3"|"$2"|"$7"|"$22}' | tail -n +3 | sed 's/ /./g') hg19 1> $pre.homer 2>$pre.homer.log &
done

echo "cancers specimens regions counts" | tr " " "\t" > sv.homer.byspeciman.stats
for f in *.homer
do
	pre=`echo $f | sed 's/.homer//'`
	awk -F"\t" '{print $1"\t"$8}' $f | tail -n +2 | sed 's/ UTR/UTR/;s/|/\t/' | awk -F"\t" '{print $1,$3}' | awk '{print $1"\t"$2}' | sort | uniq -c | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%s\t%d\n", var,$2,$3,$1}' >> sv.homer.byspeciman.stats
done

# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen.tumors.*.tsv .
# rm quantile.*.homer.KAT5.exp_*
# rm *tsv.gz.txt.specimen.tumors.*er.tsv
# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv .
for f in KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv
do
	pre=`echo $f | cut -f 1-3 -d "."`
	# class=`echo $f | cut -f 9 -d "."`
	class=`echo $f | awk -F"." '{print $(NF-1)}'`
	perl $myperl sv.homer.byspeciman.stats <(sed 's/"//g' $f) 1 1 | grep -v "/" | cut -f 1-5,8,9 | awk -v var=$class '{print $0"\t"var}' >> quantile.sv.homer.$pre.tsv
done

for f in `wc -l quantile.sv.homer.KAT5.exp_*.tsv | awk '$1 > 100 && $2!="total"{print $2}'`
do
	Rscript $rquantileboxplot $f
done

for f in `wc -l quantile.sv.homer.KAT5.exp_*.tsv | awk '$1 > 100 && $2!="total"{print $2}'`
do
	Rscript $rquantileviolin $f
done


###########
# cd /home1/04935/shaojf/stampede2/TIP60.KAT5/cnv
# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen .
echo "cancers specimens cnvtypes seqtypes counts" | tr " " "\t" > cnv.byspeciman.stats
for f in copy_number_somatic_mutation.*.tsv.gz
do
	pre=`echo $f | sed 's/copy_number_somatic_mutation.//;s/.tsv.gz//'`
	if [ ! -f KAT5.exp_array.$pre.tsv.gz.txt.specimen ] && [ ! -f KAT5.exp_seq.$pre.tsv.gz.txt.specimen ]
	then
		rm $f
	else
		gunzip -c $f | cut -f 2-3,8,20 | tail -n +2 | sed 's/ /./g' | sort | uniq -c | awk -v OFS="\t" '{print $2,$3,$4,$5,$1}' >> cnv.byspeciman.stats
	fi
done

awk '$4=="WGS"' cnv.byspeciman.stats > WGS.cnv.byspeciman.stats
awk '$4=="non-NGS"' cnv.byspeciman.stats > non-NGS.cnv.byspeciman.stats
for f in KAT5.exp_*.tsv.gz.txt.specimen
do
	pre=`echo $f | sed 's/.tsv.gz.txt.specimen//'`
	perl $myperl WGS.cnv.byspeciman.stats $f 1 1 | grep -v "/" | cut -f 1-5,8,10 > exp.WGS.cnv.$pre.tsv
	perl $myperl non-NGS.cnv.byspeciman.stats $f 1 1 | grep -v "/" | cut -f 1-5,8,10 > exp.non-NGS.cnv.$pre.tsv
done

for f in `ll exp.*.cnv.KAT5.exp_*.tsv | awk '$5 > 100{print $9}'`
do
	Rscript $rscatterplot $f
done

for pre in `ll exp.*.cnv.KAT5.exp_*.tsv | awk '$5 > 100{print $9}' | cut -f 6 -d"." | sort | uniq`
do
	f=copy_number_somatic_mutation.$pre.tsv.gz
	annotatePeaks.pl <(gunzip -c $f | awk -F"\t" -vOFS="\t" '{print "chr"$12,$13,$14,$3"|"$2"|"$15"|"$8"|"$20}' | tail -n +2 | sed 's/ /./g') hg19 1> $pre.homer 2>$pre.homer.log &
done

echo "cancers specimens regions counts" | tr " " "\t" > WGS.cnv.homer.byspeciman.stats
echo "cancers specimens regions counts" | tr " " "\t" > non-NGS.cnv.homer.byspeciman.stats
for f in *.homer
do
	pre=`echo $f | sed 's/.homer//'`
	awk -F"\t" '$1~/\|WGS/{print $1"\t"$8}' $f | sed 's/ UTR/UTR/;s/|/\t/' | awk -F"\t" '{print $1,$3}' | awk '{print $1"\t"$2}' | sort | uniq -c | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%s\t%d\n", var,$2,$3,$1}' >> WGS.cnv.homer.byspeciman.stats
	awk -F"\t" '$1~/\|non-NGS/{print $1"\t"$8}' $f | sed 's/ UTR/UTR/;s/|/\t/' | awk -F"\t" '{print $1,$3}' | awk '{print $1"\t"$2}' | sort | uniq -c | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%s\t%d\n", var,$2,$3,$1}' >> non-NGS.cnv.homer.byspeciman.stats
done

for f in KAT5.exp_*.tsv.gz.txt.specimen
do
	pre=`echo $f | sed 's/.tsv.gz.txt.specimen//'`
	perl $myperl WGS.cnv.homer.byspeciman.stats $f 1 1 | grep -v "/" | cut -f 1-5,8,9 > WGS.cnv.homer.$pre.tsv
	perl $myperl non-NGS.cnv.homer.byspeciman.stats $f 1 1 | grep -v "/" | cut -f 1-5,8,9 > non-NGS.cnv.homer.$pre.tsv
done

for f in `ll *cnv.homer.KAT5.exp_*.tsv | awk '$5 > 100{print $9}'`
do
	Rscript $rscatterplot $f
done


# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen.tumors.*.tsv .
# rm quantile.*.homer.KAT5.exp_*
# rm *tsv.gz.txt.specimen.tumors.*er.tsv
# ln -s ../exprs/KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv .
for f in KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv
do
	pre=`echo $f | cut -f 1-3 -d "."`
	# class=`echo $f | cut -f 9 -d "."`
	class=`echo $f | awk -F"." '{print $(NF-1)}'`
	perl $myperl WGS.cnv.homer.byspeciman.stats <(sed 's/"//g' $f) 1 1 | grep -v "/" | cut -f 1-5,8,9 | awk -v var=$class '{print $0"\t"var}' >> quantile.WGS.cnv.homer.$pre.tsv
	perl $myperl non-NGS.cnv.homer.byspeciman.stats <(sed 's/"//g' $f) 1 1 | grep -v "/" | cut -f 1-5,8,9 | awk -v var=$class '{print $0"\t"var}' >> quantile.non-NGS.cnv.homer.$pre.tsv
done

for f in `wc -l quantile.*cnv.homer.KAT5.exp_*.tsv | awk '$1 > 100 && $2!="total"{print $2}'`
do
	Rscript $rquantileboxplot $f
done

for f in `wc -l quantile.*cnv.homer.KAT5.exp_*.tsv | awk '$1 > 100 && $2!="total"{print $2}'`
do
	Rscript $rquantileviolin $f
done

##########################################################
# cd /home1/04935/shaojf/stampede2/TIP60.KAT5/lianghan.enhancers
ln -s ../muts/simple_somatic_mutation.open.*.bed .
ln -s ../exprs/KAT5.exp_*.tsv.gz.txt.specimen.*.tumors.*er.tsv .

for pre in `ls KAT5.exp_* | cut -f 3 -d"." | sort | uniq`
do
	f=simple_somatic_mutation.open.$pre.tsv.gz.srt.bed
	bedtools intersect -wo -a <(awk '{print "chr"$0}' $f) -b lianghan.2018cell.enhancer.bed > lianghan.2018cell.enhancer.$pre &
done

echo "cancers specimens regions counts" | tr " " "\t" > lianghan.2018cell.enhancer.byspeciman.stats
for f in lianghan.2018cell.enhancer.*-*
do
	pre=`echo $f | sed 's/lianghan.2018cell.enhancer.//'`
	cut -f 4 $f | cut -f1 -d"|" | sort | uniq -c | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%s\t%d\n", var,$2,"enhancer",$1}' >> lianghan.2018cell.enhancer.byspeciman.stats
done


for f in KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv
do
	pre=`echo $f | cut -f 1-3 -d "."`
	class=`echo $f | awk -F"." '{print $(NF-1)}'`
	perl $myperl lianghan.2018cell.enhancer.byspeciman.stats <(sed 's/"//g' $f) 1 1 | grep -v "/" | cut -f 1-5,8,9 | awk -v var=$class '{print $0"\t"var}' >> lianghan.2018cell.enhancer.counts.$pre.tsv
done

for f in `wc -l lianghan.2018cell.enhancer.counts.*.tsv | awk '$1 > 10 && $2!="total"{print $2}'`
do
	Rscript $rquantileboxplot $f
done

##########################################################
# cd /home1/04935/shaojf/stampede2/TIP60.KAT5/HEDD.enhancers
ln -s ../muts/simple_somatic_mutation.open.*.bed .
ln -s ../exprs/KAT5.exp_*.tsv.gz.txt.specimen.*.tumors.*er.tsv .

for pre in `ls KAT5.exp_* | cut -f 3 -d"." | sort | uniq`
do
	f=simple_somatic_mutation.open.$pre.tsv.gz.srt.bed
	bedtools intersect -wo -a <(awk '{print "chr"$0}' $f) -b Enhancer.bed > HEDD.enhancer.$pre &
done

echo "cancers specimens regions counts" | tr " " "\t" > HEDD.enhancer.byspeciman.stats
for f in HEDD.enhancer.*-*
do
	pre=`echo $f | sed 's/HEDD.enhancer.//'`
	cut -f 4 $f | cut -f1 -d"|" | sort | uniq -c | awk -v var=$pre -vOFS="\t" '{printf "%s\t%s\t%s\t%d\n", var,$2,"enhancer",$1}' >> HEDD.enhancer.byspeciman.stats
done

for f in KAT5.exp_*tsv.gz.txt.specimen.*.tumors.*er.tsv
do
	pre=`echo $f | cut -f 1-3 -d "."`
	class=`echo $f | awk -F"." '{print $(NF-1)}'`
	perl $myperl HEDD.enhancer.byspeciman.stats <(sed 's/"//g' $f) 1 1 | grep -v "/" | cut -f 1-5,8,9 | awk -v var=$class '{print $0"\t"var}' >> HEDD.enhancer.counts.$pre.tsv
done

for f in `wc -l HEDD.enhancer.counts.*.tsv | awk '$1 > 10 && $2!="total"{print $2}'`
do
	Rscript $rquantileboxplot $f
done
