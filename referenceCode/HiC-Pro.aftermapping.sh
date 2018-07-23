#!/bin/bash
workdir=/home1/04935/shaojf/scratch/HiChIP.test
sra=SRR5831512
HiCProInputDir=$workdir/HiCPro.input/$sra
HiCProOutputDir=$workdir/HiCPro.output/$sra
cd $HiCProOutputDir
FASTQFILE=$HiCProOutputDir/inputfiles_.txt; export FASTQFILE
# make --file /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/home1/04935/shaojf/scratch/HiChIP.test/config-hicpro.txt CONFIG_SYS=/home1/04935/shaojf/bin/HiC-Pro_2.10.0/config-system.txt all_sub 2>&1
make --file /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/home1/04935/shaojf/scratch/HiChIP.test/config-hicpro.txt CONFIG_SYS=/home1/04935/shaojf/bin/HiC-Pro_2.10.0/config-system.txt proc_hic 2>&1
make --file /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/home1/04935/shaojf/scratch/HiChIP.test/config-hicpro.txt CONFIG_SYS=/home1/04935/shaojf/bin/HiC-Pro_2.10.0/config-system.txt all_persample 2>&1
# /home1/04935/shaojf/anaconda2/bin/python /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/mergeSAM.py -q 0 -t -v -f bowtie_results/bwt2/$sra/${sra}_1_hg19.bwt2merged.bam -r bowtie_results/bwt2/$sra/${sra}_2_hg19.bwt2merged.bam -o bowtie_results/bwt2/$sra/${sra}_hg19.bwt2pairs.bam > logs/$sra/mergeSAM.log
# # fastq Quality checks - Mapping results ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' bwtDir='bowtie_results/bwt2/$sra' sampleName='$sra' r1tag='_1' r2tag='_2'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_mapping_portion.R logs/$sra/plot_mapping_portion.Rout
# # Quality Cheks - Pairing results ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' bwtDir='bowtie_results/bwt2/$sra' sampleName='$sra' rmMulti='1' rmSingle='1'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_pairing_portion.R logs/$sra/plot_pairing_portion.Rout
# # Quality checks - Hi-C processing ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' hicDir='hic_results/data/$sra' sampleName='$sra'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_hic_fragment.R logs/$sra/plot_hic_fragment.Rout
# # Quality checks - Hi-C contact maps ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' hicDir='hic_results/data/$sra' sampleName='$sra'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_hic_contacts.R logs/$sra/plot_hic_contacts.Rout
# # fastq Quality checks - Mapping results ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' bwtDir='bowtie_results/bwt2/$sra' sampleName='$sra' r1tag='_1' r2tag='_2'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_mapping_portion.R logs/$sra/plot_mapping_portion.Rout
# # Quality Cheks - Pairing results ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' bwtDir='bowtie_results/bwt2/$sra' sampleName='$sra' rmMulti='1' rmSingle='1'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_pairing_portion.R logs/$sra/plot_pairing_portion.Rout
# # Quality checks - Hi-C processing ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' hicDir='hic_results/data/$sra' sampleName='$sra'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_hic_fragment.R logs/$sra/plot_hic_fragment.Rout
# # Quality checks - Hi-C contact maps ...
# /home1/04935/shaojf/anaconda2/bin/R CMD BATCH --no-save --no-restore "--args picDir='hic_results/pic/$sra' hicDir='hic_results/data/$sra' sampleName='$sra'" /home1/04935/shaojf/bin/HiC-Pro_2.10.0/scripts/plot_hic_contacts.R logs/$sra/plot_hic_contacts.Rout
