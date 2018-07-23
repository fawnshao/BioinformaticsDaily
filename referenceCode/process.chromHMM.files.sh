#!/bin/bash
# /home1/04935/shaojf/stampede2/refs/ChromHMM/egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final
#  DATA SOURCE
# 12-mark/127-reference epigenome/25-state Imputation Based Chromatin State Model

# The states are as follows
# STATE NO.	MNEMONIC	DESCRIPTION	COLOR NAME	COLOR CODE
# 1	TssA	Active TSS	Red	255,0,0
# 2	PromU	Promoter Upstream TSS	Orange Red	255,69,0
# 3	PromD1	Promoter Downstream TSS 1	Orange Red	255,69,0
# 4	PromD2	Promoter Downstream TSS 2	Orange Red	255,69,0
# 5	Tx5	Transcribed - 5' preferential	Green	0,128,0
# 6	Tx	Strong transcription	Green	0,128,0
# 7	Tx3	Transcribed - 3' preferential	Green	0,128,0
# 8	TxWk	Weak transcription	Lighter Green	0,150,0
# 9	TxReg	Transcribed & regulatory (Prom/Enh)	Electric Lime	194,225,5
# 10	TxEnh5	Transcribed 5' preferential and Enh	Electric Lime	194,225,5
# 11	TxEnh3	Transcribed 3' preferential and Enh	Electric Lime	194,225,5
# 12	TxEnhW	Transcribed and Weak Enhancer	Electric Lime	194,225,5
# 13	EnhA1	Active Enhancer 1	Orange	255,195,77
# 14	EnhA2	Active Enhancer 2	Orange	255,195,77
# 15	EnhAF	Active Enhancer Flank	Orange	255,195,77
# 16	EnhW1	Weak Enhancer 1	Yellow	255,255,0
# 17	EnhW2	Weak Enhancer 2	Yellow	255,255,0
# 18	EnhAc	Primary H3K27ac possible Enhancer	Yellow	255,255,0
# 19	DNase	Primary DNase	Lemon	255,255,102
# 20	ZNF/Rpts	ZNF genes & repeats	Aquamarine	102,205,170
# 21	Het	Heterochromatin	Light Purple	138,145,208
# 22	PromP	Poised Promoter	Pink	230,184,183
# 23	PromBiv	Bivalent Promoter	Dark Purple	112,48,160
# 24	ReprPC	Repressed Polycomb	Gray	128,128,128
# 25	Quies	Quiescent/Low	White	255,255,255
# A chromatin state model based on the imputed data for 12 marks, H3K4me1, H3K4me2, H3K4me3, H3K9ac, H3K27ac, H4K20me1, H3K79me2, H3K36me3, H3K9me3, H3K27me3, H2A.Z, and DNase, across all 127 reference epigenomes with 25-states was learned. Within the directory for this model one can find both plain bed files containing the segmentation information and versions that can be viewed in the browser. An additional chromatin state model based on 29-marks across just the seven class 1 reference epigenomes was also learned.

myperl=/home1/04935/shaojf/myTools/BioinformaticsDaily/textProcess/add_any_2files_together.pl

for f in *_25_imputed12marks_dense.bed.gz
do
	id=`echo $f | awk -F"_" '{print $1}'`
	name=`awk -v var=$id -F"\t" '$1==var{print $2}' EIDlegend.txt`
	experiment=$id".$name"
	gunzip -c $f | cut -f 1-4 | tail -n +2 | awk -v var="$experiment" '{print $0"\t"var}' >> hg19.roadmap.25_imputed12marks.bed
done

# /home1/04935/shaojf/stampede2/refs/ChromHMM/encodeDCC.wgEncodeBroadHmm
for f in wgEncodeBroadHmm*HMM.bed.gz
do
	id=`echo $f | sed 's/wgEncodeBroadHmm//;s/HMM.bed.gz//'`
	experiment="wgEncodeBroadHmm."$id
	gunzip -c $f | cut -f 1-4 | awk -v var="$experiment" '{print $0"\t"var}' >> hg19.wgEncodeBroadHmm.bed
done

