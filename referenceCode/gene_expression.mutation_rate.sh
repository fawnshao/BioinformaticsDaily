#!/bin/sh
# rclone sync /home1/04935/shaojf/stampede2/TIP60.KAT5/ mygoogle:TIP60.KAT5/
gene=KAT5
grep KAT5 ~/stampede2/refs/hg19.refGene.tss.uniq.srt.bed > KAT5.hg19.refGene.tss
grep KAT5 <(gunzip -c ~/stampede2/refs/GENCODE/gencode.v19.annotation.gtf.gz) > KAT5.gencode.v19.annotation.gtf

cut -f 4 KAT5.hg19.refGene.tss | tr ";" "\n" | tr "|" "\n" | sort | uniq > KAT5.ids
echo "ENSG00000172977" >> KAT5.ids
awk '{print $12}' KAT5.gencode.v19.annotation.gtf | sed 's/"//g;s/;//' | sort | uniq | grep ENST | cut -f 1 -d"." >> KAT5.ids

for f in exp_*.tsv.gz
do
	gunzip -c $f | cut -f 1-2,8-9 | grep -wf KAT5.ids > KAT5.$f.txt &
done

# for f in `ll KAT5.exp_* | awk '$5==0{print $9}' | sed 's/KAT5.//;s/.txt//'`
# do
# 	gunzip -c $f | cut -f 1-2,8-9 | grep -wf KAT5.ids > KAT5.$f.txt
# done

