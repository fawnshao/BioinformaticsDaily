setwd("~/Data/myworkData/pancancer.data/")
library(data.table)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
mygenes <- fread(input = "NUPs.gene.txt")
myexprs <- fread(input = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", header = T)
# log2(TPM+1)
# myexprs$sample[1:30]
# [1] "100130426" "100133144" "100134869" "10357"     "10431"     "136542"    "155060"    "26823"     "280660"    "317712"   
# [11] "340602"    "388795"    "390284"    "391343"    "391714"    "404770"    "441362"    "442388"    "553137"    "57714"    
# [21] "645851"    "652919"    "653553"    "728045"    "728603"    "728788"    "729884"    "8225"      "90288"     "A1BG"    
# genes like MARCH* SEPT* is coverted into number!
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# Code	Definition	Short Letter Code
# 01	Primary Solid Tumor	TP
# 02	Recurrent Solid Tumor	TR
# 03	Primary Blood Derived Cancer - Peripheral Blood	TB
# 04	Recurrent Blood Derived Cancer - Bone Marrow	TRBM
# 05	Additional - New Primary	TAP
# 06	Metastatic	TM
# 07	Additional Metastatic	TAM
# 08	Human Tumor Original Cells	THOC
# 09	Primary Blood Derived Cancer - Bone Marrow	TBM
# 10	Blood Derived Normal	NB
# 11	Solid Tissue Normal	NT
# 12	Buccal Cell Normal	NBC
# 13	EBV Immortalized Normal	NEBV
# 14	Bone Marrow Normal	NBM
# 15	sample type 15	15SH
# 16	sample type 16	16SH
# 20	Control Analyte	CELLC
# 40	Recurrent Blood Derived Cancer - Peripheral Blood	TRB
# 50	Cell Lines	CELL
# 60	Primary Xenograft Tissue	XP
# 61	Cell Line Derived Xenograft Tissue	XCL
# 99	sample type 99	99SH
myinfo <- fread(input = "Survival_SupplementalTable_S1_20171025_xena_sp", header = T)
mygenes.exprs <- myexprs[sample %in% mygenes$targets]
mygenes.exprs.mt <- melt.data.table(mygenes.exprs, id.vars = "sample")
mygenes.exprs.mt[, Cancer := myinfo[match(mygenes.exprs.mt$variable, Tumor_Sample_Barcode)]$`cancer type abbreviation`]
mygenes.exprs.mt[, Stage := myinfo[match(mygenes.exprs.mt$variable, Tumor_Sample_Barcode)]$ajcc_pathologic_tumor_stage]
mygenes.exprs.mt[, STcode := substr(variable, 14,15)]
# unique(mygenes.exprs.mt$STcode)
# [1] "01" "11" "02" "06" "07" "05" "03"
mygenes.exprs.mt[STcode == "01", Type := "Primary Solid Tumor"]
mygenes.exprs.mt[STcode == "02", Type := "Recurrent Solid Tumor"]
mygenes.exprs.mt[STcode == "03", Type := "Primary Blood Derived Cancer"]
mygenes.exprs.mt[STcode == "05", Type := "Additional - New Primary"]
mygenes.exprs.mt[STcode == "06", Type := "Metastatic"]
mygenes.exprs.mt[STcode == "07", Type := "Additional Metastatic"]
mygenes.exprs.mt[STcode == "11", Type := "Solid Tissue Normal"]


mygenes.exprs.mt$Type <- factor(mygenes.exprs.mt$Type, 
                                levels = c("Solid Tissue Normal", "Primary Solid Tumor", "Recurrent Solid Tumor", "Primary Blood Derived Cancer", 
                                           "Additional - New Primary", "Metastatic", "Additional Metastatic"))
cancerwithnormal <- unique(na.omit(mygenes.exprs.mt[STcode=="11"]$Cancer))
mygenes.exprs.mt$Cancer <- factor(mygenes.exprs.mt$Cancer, levels = sort(unique(mygenes.exprs.mt$Cancer)))
p1 <- ggboxplot(data = mygenes.exprs.mt[sample == "NUP93" & Cancer %in% cancerwithnormal & STcode %in% c("01","11")], 
                x = "Cancer", y = "value", 
                color = "Type", fill = "Type", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
                title = "NUP93 in TCGA (with matched normal tissues)")
p11 <- ggpar(p1, x.text.angle = 45)
p2 <- ggboxplot(data = mygenes.exprs.mt[sample == "CTCF" & Cancer %in% cancerwithnormal & STcode %in% c("01","11")], 
                x = "Cancer", y = "value", 
                color = "Type", fill = "Type", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
                title = "CTCF in TCGA (with matched normal tissues)")
p22 <- ggpar(p2, x.text.angle = 45)
ggsave(filename = "ggboxplot.NUP93.CTCF.expression.pdf", width = 12, height = 10, 
       ggarrange(p11, p22, nrow = 2))


################
library(survival)
library(survminer)
for(i in 1:length(cancers)){
    for(j in 1:length(genes.withCNV)){
        mygenes.cnv.cancer <- mygenes.cnvs.mt[Sample == genes.withCNV[j] & Cancer == cancers[i]]
        mygenes.cnv.cancer[value < 0, CNV  := "Del"]
        mygenes.cnv.cancer[value > 0, CNV  := "Amp"]
        mygenes.cnv.cancer[value == 0, CNV  := "NoChange"]
        mygenes.cnv.cancer[, OS := myinfo[match(mygenes.cnv.cancer$variable, Tumor_Sample_Barcode)]$OS]
        mygenes.cnv.cancer[, OS.time := myinfo[match(mygenes.cnv.cancer$variable, Tumor_Sample_Barcode)]$OS.time]
        # mygenes.cnv.cancer[, DFI := myinfo[match(mygenes.cnv.cancer$variable, Tumor_Sample_Barcode)]$DFI]
        # mygenes.cnv.cancer[, DFI.time := myinfo[match(mygenes.cnv.cancer$variable, Tumor_Sample_Barcode)]$DFI.time]
        # mygenes.cnv.cancer[, PFI := myinfo[match(mygenes.cnv.cancer$variable, Tumor_Sample_Barcode)]$PFI]
        # mygenes.cnv.cancer[, PFI.time := myinfo[match(mygenes.cnv.cancer$variable, Tumor_Sample_Barcode)]$PFI.time]
        fit <- survfit(Surv(OS.time, OS) ~ CNV, data = mygenes.cnv.cancer)
        p1 <- ggsurvplot(fit,
                         pval = TRUE, conf.int = FALSE,
                         risk.table = TRUE, 
                         risk.table.col = "strata", 
                         linetype = "strata", 
                         # surv.median.line = "hv",
                         # ggtheme = theme_bw(),
                         palette = "npg", title = paste(genes.withCNV[j], cancers[i], sep = " in "))
        pdf(file = paste("batch.ggsurvplot.CNV", cancers[i], genes.withCNV[j], "pdf", sep = "."), width = 10, height = 8)
        print(p1)
        dev.off()
    }
}