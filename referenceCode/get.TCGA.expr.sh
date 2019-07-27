#!/bin/bash
genelist=$1
infocol=$2
expr=/home1/04935/shaojf/stampede2/refs/xenahubs/cohort.TCGA.PanCanAtlas/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena
annotation=/home1/04935/shaojf/stampede2/refs/xenahubs/cohort.TCGA.PanCanAtlas/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info
# head_line EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info 
# 1 sample
# 2 _PATIENT
# 3 cancer type abbreviation
# 4 age_at_initial_pathologic_diagnosis
# 5 gender
# 6 race
# 7 ajcc_pathologic_tumor_stage
# 8 clinical_stage
# 9 histological_type
# 10 histological_grade
# 11 initial_pathologic_dx_year
# 12 menopause_status
# 13 birth_days_to
# 14 vital_status
# 15 tumor_status
# 16 last_contact_days_to
# 17 death_days_to
# 18 cause_of_death
# 19 new_tumor_event_type
# 20 new_tumor_event_site
# 21 new_tumor_event_site_other
# 22 new_tumor_event_dx_days_to
# 23 treatment_outcome_first_course
# 24 margin_status
# 25 residual_tumor
# 26 _EVENT
# 27 _TIME_TO_EVENT
# 28 OS
# 29 OS.time
# 30 DSS
# 31 DSS.time
# 32 DFI
# 33 DFI.time
# 34 PFI
# 35 PFI.time
# 36 Redaction
# awk -F"\t" -vOFS="\t" '{a=$7;if($7==""){a=$8;}print $1,$3,a}' EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info > xena.stage.info.txt

head -1 $expr > xena.mygene.txt
grep -wf $genelist $expr >> xena.mygene.txt
# cut -f 1,3,7-8 EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.info | grep PAAD
cut -f 1,3,$infocol
