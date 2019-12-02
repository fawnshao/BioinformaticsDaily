setwd("~/Data/myworkData/pancancer.data/")
library(data.table)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
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
# We then define a hotspot as an amino acid position in a protein-coding gene mutated more frequently than would be expected in the absence of selection.
# https://www.nature.com/articles/nbt.3391
mygenes <- fread(input = "NUPs.gene.txt")
myhotspots <- fread("nbt.3391.taylorlab.hotspots/publication_hotspots.tsv", header = T)
myexprs <- fread(input = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", header = T)
# log2(TPM+1)
# myexprs$sample[1:30]
# [1] "100130426" "100133144" "100134869" "10357"     "10431"     "136542"    "155060"    "26823"     "280660"    "317712"   
# [11] "340602"    "388795"    "390284"    "391343"    "391714"    "404770"    "441362"    "442388"    "553137"    "57714"    
# [21] "645851"    "652919"    "653553"    "728045"    "728603"    "728788"    "729884"    "8225"      "90288"     "A1BG"    
# genes like MARCH* SEPT* is converted into number!
mycnvs <- fread(input = "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", header = T)
mymuts <- fread(input = "formatted.mc3.v0.2.8.PUBLIC.xena.maf", header = T)
myinfo <- fread(input = "Survival_SupplementalTable_S1_20171025_xena_sp", header = T)
# length(unique(mygenes.exprs.mt$variable))
# length(unique(colnames(myexprs)[-1]))
# 11060
# length(unique(mymuts$Tumor_Sample_Barcode))
# 9104
# length(unique(colnames(mycnvs)[-1]))
# 10845
# length(unique(myinfo$Tumor_Sample_Barcode))
# 12591
mutationload <- as.data.table(table(unique(mymuts[, c(2:4,9)])$Tumor_Sample_Barcode))
# mutationload.genes <- as.data.table(table(unique(mymuts[, c(1:4,9)])[,c(1,5)]))
# mutationload.genes <- mutationload.genes[N > 0]

# from cBioportal
# TCGA-D1-A17Q-01 5945
# TCGA-E6-A1LX-01 12696
# > mutationload[V1=="TCGA-D1-A17Q-01"]
# V1    N
# 1: TCGA-D1-A17Q-01 8578
# > mutationload[V1=="TCGA-E6-A1LX-01"]
# V1     N
# 1: TCGA-E6-A1LX-01 22062
# https://dcc.icgc.org/donors/DO41903/mutations
# 8628 in ICGC
mygenes.cnvs <- mycnvs[Sample %in% mygenes$targets]
mygenes.cnvs.mt <- melt.data.table(mygenes.cnvs, id.vars = "Sample")
mygenes.exprs <- myexprs[sample %in% mygenes$targets]
mygenes.exprs.mt <- melt.data.table(mygenes.exprs, id.vars = "sample")
mygenes.mut.mt <- as.data.table(table(unique(mymuts[Hugo_Symbol %in% mygenes$targets, c(1:4,9)])[,c(1,5)]))

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

colnames(mygenes.exprs.mt)[1] <- "Gene"
colnames(mygenes.cnvs.mt)[1] <- "Gene"
mygenes.exprs.cnv <- merge(x = mygenes.exprs.mt, y = mygenes.cnvs.mt, all.x = TRUE, by = c("Gene", "variable"))
# mygenes.exprs.cnv <- merge(mygenes.exprs.mt, mygenes.cnvs.mt, by.x = c("sample", "variable"), by.y = c("Sample", "variable"))
mygenes.exprs.cnv.mut <- merge(x = mygenes.exprs.cnv, y = mygenes.mut.mt, 
                               all.x = TRUE, by.x = c("Gene", "variable"), by.y = c("Hugo_Symbol", "Tumor_Sample_Barcode"))
colnames(mygenes.exprs.cnv.mut) <- c("Gene", "Tumor_Sample_Barcode", "RNA", "Cancer", "Stage", "STcode", "Type", "CNV", "Mutation")
# mygenes.exprs.cnv.mut[Mutation > 1]
mygenes.exprs.cnv.mut[, MutationLoad := mutationload[match(mygenes.exprs.cnv.mut$Tumor_Sample_Barcode, V1)]$N]

mygenes.exprs.cnv.mut$Type <- factor(mygenes.exprs.cnv.mut$Type, 
                                     levels = c("Solid Tissue Normal", "Primary Solid Tumor", 
                                                "Recurrent Solid Tumor", "Primary Blood Derived Cancer", 
                                                "Additional - New Primary", "Metastatic", "Additional Metastatic"))
cancerwithnormal <- unique(na.omit(mygenes.exprs.cnv.mut[STcode=="11"]$Cancer))
mygenes.exprs.cnv.mut$Cancer <- factor(mygenes.exprs.cnv.mut$Cancer, levels = sort(unique(mygenes.exprs.cnv.mut$Cancer)))
mygenes.exprs.cnv.mut[is.na(CNV), CNVStatus := "Unknown"]
mygenes.exprs.cnv.mut[CNV < 0, CNVStatus := "Del"]
mygenes.exprs.cnv.mut[CNV == 0, CNVStatus := "WT"]
mygenes.exprs.cnv.mut[CNV > 0, CNVStatus := "Amp"]
mygenes.exprs.cnv.mut[is.na(MutationLoad), MutationStatus := "Unknown"]
mygenes.exprs.cnv.mut[is.na(Mutation) & MutationLoad > 0, MutationStatus := "WT"]
mygenes.exprs.cnv.mut[Mutation == 0, MutationStatus := "WT"]
mygenes.exprs.cnv.mut[Mutation == 1, MutationStatus := "singleMut"]
mygenes.exprs.cnv.mut[Mutation > 1, MutationStatus := "multipleMut"]
mygenes.exprs.cnv.mut$CNVStatus <- factor(mygenes.exprs.cnv.mut$CNVStatus, levels = c("Del", "WT", "Amp", "Unknown"))
mygenes.exprs.cnv.mut$MutationStatus <- factor(mygenes.exprs.cnv.mut$MutationStatus, levels = c("WT", "singleMut", "multipleMut", "Unknown"))

mygenes.exprs.cnv.mut[, PATIENT := myinfo[match(mygenes.exprs.cnv.mut$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$`_PATIENT`]
mygenes.exprs.cnv.mut[, gene.pid := paste(Gene, PATIENT)]
mygenes.exprs.cnv.mut.t <- mygenes.exprs.cnv.mut[STcode != "11"]
mygenes.exprs.cnv.mut[STcode == "11", CNVStatus := mygenes.exprs.cnv.mut.t[match(mygenes.exprs.cnv.mut[STcode == "11"]$gene.pid, gene.pid)]$CNVStatus]
normcounts <- table(mygenes.exprs.cnv.mut[STcode == "11" & Gene == "NUP93"]$Cancer)
p1 <- ggboxplot(data = mygenes.exprs.cnv.mut[Gene == "NUP93" & PATIENT %in% unique(mygenes.exprs.cnv.mut[STcode=="11"]$PATIENT) & 
                                                 STcode %in% c("01","11") & Cancer %in% names(normcounts[normcounts > 20]) & 
                                                 !is.na(CNVStatus) & CNVStatus != "Unknown"], 
                x = "Type", y = "RNA", 
                color = "CNVStatus", fill = "CNVStatus", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
                add = "jitter", width = 0.8, 
                title = "NUP93 in TCGA (with matched normal tissues & normal samples > 20)") + border()
# p11 <- ggpar(facet(p1, facet.by = "Cancer", nrow = 1), x.text.angle = 45)
p11 <- ggpar(facet(p1, facet.by = "Cancer", nrow = 1), font.xtickslab = NA)
p2 <- ggboxplot(data = mygenes.exprs.cnv.mut[Gene == "CTCF" & PATIENT %in% unique(mygenes.exprs.cnv.mut[STcode=="11"]$PATIENT) & 
                                                 STcode %in% c("01","11") & Cancer %in% names(normcounts[normcounts > 20]) & 
                                                 !is.na(CNVStatus) & CNVStatus != "Unknown"], 
                x = "Type", y = "RNA", 
                color = "CNVStatus", fill = "CNVStatus", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
                add = "jitter", width = 0.8, 
                title = "CTCF in TCGA (with matched normal tissues & normal samples > 20)") + border()
p22 <- ggpar(facet(p2, facet.by = "Cancer", nrow = 1), font.xtickslab = NA)
ggsave(filename = "ggboxplot.matchedNT.NUP93.CTCF.expression.pdf", width = 25, height = 10, 
       ggarrange(p11, p22, nrow = 2))

# THCA no del * decreased???
paired <- mygenes.exprs.cnv.mut[PATIENT %in% unique(mygenes.exprs.cnv.mut[STcode=="11"]$PATIENT) & 
                                    PATIENT %in% unique(mygenes.exprs.cnv.mut[STcode=="01"]$PATIENT) & 
                                    Cancer %in% names(normcounts[normcounts > 20]) & STcode %in% c("01","11")]
paired.t <- paired[STcode == "01"]
paired.n <- paired[STcode == "11"]
paired.nt <- merge(x = paired.n, y = paired.t, all = TRUE, by = c("Gene", "PATIENT"))
paired.nt[, FC := RNA.y - RNA.x]
p1 <- ggboxplot(data = paired.nt[Gene == "NUP93" & CNVStatus.y != "Unknown"], 
                x = "Cancer.y", y = "FC", 
                color = "CNVStatus.y", fill = "CNVStatus.y", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(T/N)", outlier.size = 0.2,
                add = "jitter", width = 0.8, 
                title = "NUP93 in TCGA (with matched normal tissues & normal samples > 20)") + geom_hline(yintercept = 0, linetype = 2)
p2 <- ggboxplot(data = paired.nt[Gene == "CTCF" & CNVStatus.y != "Unknown"], 
                x = "Cancer.y", y = "FC", 
                color = "CNVStatus.y", fill = "CNVStatus.y", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(T/N)", outlier.size = 0.2,
                add = "jitter", width = 0.8, 
                title = "CTCF in TCGA (with matched normal tissues & normal samples > 20)") + geom_hline(yintercept = 0, linetype = 2)
ggsave(filename = "ggboxplot.matchedNT.NUP93.CTCF.log2FC.pdf", width = 15, height = 10, 
       ggarrange(p1, p2, nrow = 2))


p1 <- ggboxplot(data = mygenes.exprs.cnv.mut[Gene == "NUP93" & Cancer == "BRCA" & STcode %in% c("01","11")], 
                x = "Type", y = "RNA", 
                color = "CNVStatus", fill = "CNVStatus", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
                add = "jitter", width = 0.8, 
                title = "NUP93 in TCGA (with matched normal tissues)")
# p11 <- ggpar(facet(p1, facet.by = "Cancer", nrow = 1), x.text.angle = 45)
p2 <- ggboxplot(data = mygenes.exprs.cnv.mut[Gene == "CTCF" & Cancer == "BRCA" & STcode %in% c("01","11")], 
                x = "Type", y = "RNA", 
                color = "CNVStatus", fill = "CNVStatus", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
                add = "jitter", width = c(0.2, 0.8, 0.8, 0.8), 
                title = "CTCF in TCGA (with matched normal tissues)")
# p22 <- ggpar(p2, x.text.angle = 45)
ggsave(filename = "ggboxplot.BRCA.NUP93.CTCF.expression.pdf", width = 5, height = 10, 
       ggarrange(p1, p2, nrow = 2))

simdata <- mygenes.exprs.cnv.mut[Gene == "NUP93" & Cancer == "BRCA" & STcode == "01" & !is.na(MutationLoad)]
p1 <- ggscatter(simdata, x = "MutationLoad", y = "RNA",
                color = "CNVStatus", fill = "CNVStatus", shape = "MutationStatus")
ggpar(p1, xscale = "log10")
setorder(simdata, MutationLoad, CNV, RNA)
simdata[,Ranks := seq(1:nrow(simdata))]
p1 <- ggscatter(simdata, x = "Ranks", y = "RNA",
                color = "CNVStatus", fill = "CNVStatus", shape = "MutationStatus")
p1 <- ggscatter(simdata, x = "Ranks", y = "MutationLoad",
                color = "CNVStatus", fill = "CNVStatus", shape = "MutationStatus", 
                size = 0.8, label = "MutationStatus", label.select = list(criteria = "`Mutation` > 0"))
ggsave(filename = "ggscatter.BRCA.NUP93.mutation.cnv.expression.pdf", width = 15, height = 6, p1)
p1 <- ggboxplot(simdata, x = "Mutation", y = "RNA",
                color = "CNVStatus", fill = "CNVStatus", palette = "aaas", alpha = 0.8,
                xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
                add = "jitter")
# simdata[, PATIENT := myinfo[match(simdata$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$`_PATIENT`]
# simdata.t <- simdata[STcode != "11"]
# simdata[STcode == "11", CNVStatus := simdata.t[match(simdata[STcode == "11"]$PATIENT, PATIENT)]$CNVStatus]
# p1 <- ggboxplot(data = simdata[PATIENT %in% unique(simdata[STcode=="11"]$PATIENT) & STcode %in% c("01","11")],
#                 x = "Type", y = "RNA",
#                 color = "CNVStatus", fill = "CNVStatus", palette = "aaas", alpha = 0.8,
#                 xlab = "", ylab = "log2(TPM+1)", outlier.size = 0.2,
#                 add = "jitter", width = 0.8,
#                 title = "NUP93 in TCGA (with matched normal tissues)")

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