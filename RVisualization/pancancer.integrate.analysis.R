setwd("~/Data/myworkData/pancancer.data/")
library(data.table)
library(ggpubr)
library(gridExtra)
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

# https://www.cell.com/cell/fulltext/S0092-8674(17)31142-X?code=cell-site
# Having validated panel-based hypermutation testing, we examined the mutation burden in 2,885 pediatric tumors. 
# Mutation frequency ranged from 0–864 Mut/Mb (Figure 1A), with a mean and median of 6.78 Mut/Mb and 2.50 Mut/Mb, respectively. 
# Using segmented linear regression analysis, we calculated 9.91 and 9.0 Mut/Mb as appropriate thresholds for hypermutation in childhood and adult cancers (Figures S1C and S1D; STAR Methods). 
# For consistency, we use 10 Mut/Mb to define hypermutation in both cohorts. We also note that this coincides with the median mutation burden of patients previously reported to respond to checkpoint inhibition
# 10*3000

# xena pancancer maf is from https://www.ncbi.nlm.nih.gov/pubmed/?term=29596782
# over 10,000 tumor-normal exome pairs across 33 different cancer
# Variants from: MuTect, MuSE, VarScan2, Radia, Pindel, Somatic Sniper, Indelocator
# Variants called in non-exonic regions, such as introns, 5' or 3' UTR are restricted to controlled-access release.
# somatic variants that occur in exonic regions in open-access files
# The exon definitions were derived from the GAF 4.0 definition, which was based on Gencode 19 Basic.
# more gencode.v19.genes.v7.patched_contigs.gtf | awk '$3=="exon"' | cut -f 1,4-5 | sort | uniq | awk 'BEGIN{sum=0}{sum+=$3-$2+1}END{print sum}'
# 107975757
hypermutatecut <- 10 * 107975757 / 1e6

# https://elifesciences.org/articles/37294
# The CNA burden of a tumor is the degree to which a tumor's genome is altered as 
# a percentage of genome length and represents a fundamental measure of genome copy number alteration level.
# or copy nunmber burden (CNB)

# We then define a hotspot as an amino acid position in a protein-coding gene mutated more frequently than would be expected in the absence of selection.
# https://www.nature.com/articles/nbt.3391
mygenes <- fread(input = "NUPs.gene.txt")
myhotspots <- fread("nbt.3391.taylorlab.hotspots/publication_hotspots.tsv", header = T)
# myhotspots[`Hugo Symbol` %in% mygenes$targets]
# Hugo Symbol Codon Alt Common Codon Usage * Variant Amino Acid  Q-value Tumor Count Tumor Type Count Validation Level [a]
# 1:       NUP93   E14                     <NA>           K:10|G:1 1.59e-10          11                6              Level-1
# 2:      POM121   N94                     N359                S:7 1.23e-08           7                4              Level-1
# 3:       NUP93   Q15                     <NA>                *:4 8.20e-03           4                2              Level-1
# Tumor Type Composition
# 1: brca:6|thca:1|luad:1|lihc:1|hnsc:1|blca:1
# 2:                paad:4|skcm:1|prad:1|gbm:1
# 3:                             thca:3|ucec:1

# TCGA-44-A4SS-01	16	56782199	56782199	G	A	NUP93	Missense_Mutation	p.E14K	0.11	deleterious(0.05)	benign(0.129)
# TCGA-BH-A209-01	16	56782199	56782199	G	A	NUP93	Missense_Mutation	p.E14K	0.17	deleterious(0.05)	benign(0.129)
# TCGA-E2-A10C-01	16	56782199	56782199	G	A	NUP93	Missense_Mutation	p.E14K	0.33	deleterious(0.05)	benign(0.129)
# TCGA-L6-A4EU-01	16	56782199	56782199	G	A	NUP93	Missense_Mutation	p.E14K	0.07	deleterious(0.05)	benign(0.129)
# 
# bcr_patient_barcode	er_status_by_ihc	pr_status_by_ihc	her2_status_by_ihc
# TCGA-BH-A209	Positive	Positive	[Not Evaluated]
# TCGA-E2-A10C	Positive	Positive	Negative
# 
# genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.1.1.0.curated.somatic.maf
# TCGA-E9-A1NC	Negative	Positive	Positive

# head -1 EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena | tr "\t" "\n" | cat -n | grep -e "TCGA-E2-A10C" -e "TCGA-44-A4SS" -e "TCGA-BH-A209" -e "TCGA-L6-A4EU"
# 6158	TCGA-BH-A209-01
# 6159	TCGA-BH-A209-11
# 6307	TCGA-E2-A10C-01
# 7521	TCGA-L6-A4EU-01
# 8174	TCGA-44-A4SS-01

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
mutationload[N>hypermutatecut, hypermutation := "Y"]
mutationload[is.na(hypermutation), hypermutation := "N"]
mutationload[, Cancer := myinfo[match(mutationload$V1, Tumor_Sample_Barcode)]$`cancer type abbreviation`]
mutationload[, MedianN := median(N), by = "Cancer"]
sort(unique(mutationload$MedianN))
# max(mutationload$N) / 107975757 * 1e6
# 254
mutationload.count <- as.data.table(table(mutationload[,c(3:4)]))
mutationload.count[, sums := sum(N), by = "Cancer"]
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
mymuts[, Cancer := myinfo[match(mymuts$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$`cancer type abbreviation`]
mymuts[, mutationcount := mutationload[match(mymuts$Tumor_Sample_Barcode, V1)]$N]
mymuts[mutationcount > hypermutatecut, hypermutation := "HM"]
mymuts[is.na(hypermutation), hypermutation := "NotHM"]
fwrite(mymuts[hypermutation == "NotHM", c(1:10)], file = "NotHM.formatted.mc3.v0.2.8.PUBLIC.xena.maf", sep = "\t")
mymuts.samples <- unique(mymuts[,c(9,12,13)])
nup93.muts1 <- as.data.table(table(mymuts[Hugo_Symbol=="NUP93" & Amino_Acid_Change != ""]$Amino_Acid_Change))
nup93.muts1 <- nup93.muts1[N > 0]
setorder(nup93.muts1, -N)
nup93.muts1[N > 1]
# nup93.muts <- as.data.table(melt(table(mymuts[Hugo_Symbol=="NUP93" & Amino_Acid_Change != "", c(10,11)])))
nup93.muts <- as.data.table(table(mymuts[Hugo_Symbol=="NUP93" & Amino_Acid_Change != "", c(10,11)]))
nup93.muts <- nup93.muts[N > 0]
setorder(nup93.muts, -N)
nup93.muts[Amino_Acid_Change %in% nup93.muts1[N > 1]$V1]

###
# oncoplot --------------------------------------------------------------
library(maftools)
tcga.muts <- read.maf(maf = "NotHM.formatted.mc3.v0.2.8.PUBLIC.xena.maf", 
                      clinicalData = "Survival_SupplementalTable_S1_20171025_xena_sp")
getSampleSummary(tcga.muts)
getGeneSummary(tcga.muts)
getClinicalData(tcga.muts)
getFields(tcga.muts)
pdf(file ="NUPgenes.mc3.v0.2.8.PUBLIC.maf.oncoplot.NotHM.pdf", width = 16, height = 12)
oncoplot(maf = tcga.muts, genes = mygenes$targets)
dev.off()

png(filename ="NUPgenes.mc3.v0.2.8.PUBLIC.maf.oncoplot.NotHM.png", width = 1200, height = 800)
oncoplot(maf = tcga.muts, genes = mygenes$targets)
dev.off()

png(filename ="NUPgenes.mc3.v0.2.8.PUBLIC.maf.oncoplot.NotHM.1.png", width = 1600, height = 1000)
oncoplot(maf = tcga.muts, genes = mygenes$targets)
dev.off()

# CNV --------------------------------------------------------------
# https://xenabrowser.net/datapages/?dataset=broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# log(tumor/normal)
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
# Genes with focal CNV values smaller than -0.3 are categorized as a "loss" (-1)
# Genes with focal CNV values larger than 0.3 are categorized as a "gain" (+1)
# Genes with focal CNV values between and including -0.3 and 0.3 are categorized as "neutral" (0).
mygenes.cnvs <- mycnvs[Sample %in% mygenes$targets]
mygenes.cnvs.mt <- melt.data.table(mygenes.cnvs, id.vars = "Sample")
mysegs <- fread("broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.xena.seg", header = T)
myAneuploidyScore <- fread("10.1016_j.ccell.2018.03.007.mmc2.tableS2.csv", header = T, sep = ",")
mygenecoords <- fread("NUPs.gene.bed", header = F)
mygenecoords[,Gene := matrix(unlist(strsplit(V4, split = "\\|")), byrow = T, ncol = 2)[,2]]
# mysegs[value > log2(2.5/2)]
mysegs.sub <- mysegs[value > log2(2.5/2) | value < log2(1.5/2)]
mysegs.sub[, totalCNA := sum(floor((end-start)/1000)), by = "sampleID"]
mysegs.sub[, CNAburden := totalCNA / 3e6 * 100]
mysegs.sub.uniq <- unique(mysegs.sub[,c(1,7)])
# mysegs.sub.uniq[sampleID=="TCGA-OR-A5KO-01"]
# ggdensity(mysegs.sub.uniq, x = "CNAburden")
# ggdensity(mysegs.sub.uniq[sampleID %in% mymuts.samples[hypermutation=="HM"]$Tumor_Sample_Barcode], x = "CNAburden")
mysegs[, CNAburden := mysegs.sub.uniq[match(mysegs$sampleID, sampleID)]$CNAburden]
mysegs[, mutationcount := mymuts.samples[match(mysegs$sampleID, Tumor_Sample_Barcode)]$mutationcount]
mysegs[, hypermutation := mymuts.samples[match(mysegs$sampleID, Tumor_Sample_Barcode)]$hypermutation]
# quantile(mysegs.sub.uniq$CNAburden, na.rm = T, probs = 0.75)
# 18.65953
length(unique(mysegs[CNAburden < 20 & hypermutation == "NotHM"]$sampleID))
# length(unique(mysegs[CNAburden < 20]$sampleID))
# 8385
# tmpcnv <- unique(mygenes.exprs.cnv.mut[!is.na(CNV),c(1,2,4,8)])
# tmpcnv[grep("Normal", Type)]

# mygenes.cnvs.mat <- dcast(data = mygenes.exprs.cnv.mut[!is.na(CNV) & hypermutation == "NotHM", c(1,2,8)], 
# formula = Gene ~ Tumor_Sample_Barcode, mean, value.var = "CNV")
# mygenes.cnvs.mat[1:5,1:5]
# mygenes.cnvs.mat <- data.matrix(dcast(data = mygenes.cnvs.mt[variable %in% mysegs.sub.uniq[CNAburden < 20]$sampleID], formula = Gene ~ variable, fun.aggregate = mean)[,-1])
mycolidx <- na.omit(match(mysegs.sub.uniq[CNAburden < 20]$sampleID, colnames(mygenes.cnvs)))
# length(mycolidx)
# 8262
mygenes.cnvs.mat <- data.matrix(mygenes.cnvs[,.SD, .SDcols = mycolidx])
tmp <- colSums(mygenes.cnvs.mat)
length(tmp[tmp!=0])
rownames(mygenes.cnvs.mat) <- mygenes.cnvs$Sample
# nasums <- apply()
row.annos <- data.frame(geneCoords = mygenecoords[match(rownames(mygenes.cnvs.mat), Gene)]$V1,
                        geneClass = mygenes[match(rownames(mygenes.cnvs.mat), targets)]$Class)
rownames(row.annos) <- rownames(mygenes.cnvs.mat)
tmp <- myinfo[match(colnames(mygenes.cnvs.mat),Tumor_Sample_Barcode)]$`cancer type abbreviation`
tmp[is.na(tmp)] <- "Unkown"
col.annos <- data.frame(cancerType = tmp)
rownames(col.annos) <- colnames(mygenes.cnvs.mat)
# col.annos[col.annos$cancerType == " ", 1] <- "Unkown"
# mycols <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(5)
# mycols <- colorRampPalette(c("blue", "white", "red"))(5)
# display.brewer.pal(5,"RdYlGn")
mycols <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlGn")))(5)
a1 <- colorRampPalette(brewer.pal(n = 5, name = "PuBuGn"))(length(unique(row.annos$geneClass)))
names(a1) <- sort(levels(row.annos$geneClass))
a2 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(row.annos$geneCoords)))
# names(a2) <- sort(levels(row.annos$geneCoords))
names(a2) <- paste("chr", c(seq(1,13,1), seq(16,20,1),"22","X"), sep = "")
a3 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(col.annos$cancerType)))
names(a3) <- levels(col.annos$cancerType)
ann_colors = list(
    geneClass = a1,
    geneCoords = a2,
    cancerType = a3
    )
pdf(file = "pheatmap.CNV.allraw.nohyperCNB.1.pdf", width = 15, height = 10)
pheatmap(mygenes.cnvs.mat, annotation_row = row.annos, annotation_col = col.annos, 
         show_colnames = F, color = mycols, fontsize = 7, fontsize_row = 10, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D",
         annotation_colors = ann_colors, 
         main = "NUPs/Cohesin/Condensin altered in 6885/8262 of tumor samples with CNA burden < 20%")
dev.off()

###
mygenes.cnvs.mt[, CNAburden := mysegs.sub.uniq[match(mygenes.cnvs.mt$variable, sampleID)]$CNAburden]
mygenes.cnvs.mt[, mutationcount := mymuts.samples[match(mygenes.cnvs.mt$variable, Tumor_Sample_Barcode)]$mutationcount]
mygenes.cnvs.mt[, hypermutation := mymuts.samples[match(mygenes.cnvs.mt$variable, Tumor_Sample_Barcode)]$hypermutation]
mygenes.cnvs.mt[, cancerType := myinfo[match(mygenes.cnvs.mt$variable,Tumor_Sample_Barcode)]$`cancer type abbreviation`]
mygenes.cnvs.mt[value < 0, CNAtype := "Del"]
mygenes.cnvs.mt[value == 0, CNAtype := "WT"]
mygenes.cnvs.mt[value > 0, CNAtype := "Amp"]
# length(unique(mygenes.cnvs.mt$variable))
# 10845
# length(unique(mygenes.cnvs.mt[CNAburden < 20]$variable))
# 8262
colnames(mygenes.cnvs.mt)[1] <- "Gene"
# cnacounts <- as.data.table(table(mygenes.cnvs.mt[CNAburden < 20, c(1,7,8)]))
cnacounts <- as.data.table(table(mygenes.cnvs.mt[CNAburden < 20, c(1,8)]))
tmp <- cnacounts[CNAtype=="WT"]
setorder(tmp, N)
cnacounts$Gene <- factor(cnacounts$Gene, levels = tmp$Gene)
p1 <- ggbarplot(cnacounts, x = "Gene", y = "N", 
                color = "CNAtype", fill = "CNAtype", palette = "npg",
                title = "NUPs/Cohesin/Condensin altered in 6885/8262 of tumor samples with CNA burden < 20%",
                xlab = "", ylab = "Counts")
p11 <- ggpar(p1, x.text.angle = 45, font.x = 10)
# ggsave(filename = "ggbarplot.CNV.nohyperCNB.counts.pdf", p11, width = 15, height = 6)
ggsave(filename = "ggbarplot.CNV.nohyperCNB.counts.1.pdf", p11, width = 12, height = 6)

lowCNAburden.cancercount <- as.data.table(table(mygenes.cnvs.mt[CNAburden < 20 & Gene == "CTCF"]$cancerType))
# lowCNAburden.count <- as.data.table(table(mygenes.cnvs.mt[CNAburden < 20, c(1,3,7)]))
lowCNAburden.count0 <- as.data.table(table(mygenes.cnvs.mt[CNAburden < 20, c(1,8)]))
lowCNAburden.count <- as.data.table(table(mygenes.cnvs.mt[CNAburden < 20, c(1,8,7)]))
lowCNAburden.count[, cancerTotal := lowCNAburden.cancercount[match(lowCNAburden.count$cancerType, V1)]$N]
lowCNAburden.count[, percentage := N/cancerTotal*100]
lowCNAburden.count[, ID := paste(Gene, CNAtype)]
lowCNAburden.count0[, ID := paste(Gene, CNAtype)]
lowCNAburden.count[, totalN := lowCNAburden.count0[match(lowCNAburden.count$ID, ID)]$N]
lowCNAburden.count[, hyper.p := phyper(N, totalN, 8262 - totalN, cancerTotal, lower.tail = FALSE)]
# unique(lowCNAburden.count[hyper.p < 1e-3 & CNAtype != "WT" & cancerTotal > 50]$Gene)
sigcount <- as.data.table(table(lowCNAburden.count[hyper.p < 1e-3 & CNAtype != "WT" & cancerTotal > 50, c(1,2)]))
# hypergeometric test is useless
# delplot1 <- ggplot(data = lowCNAburden.count[value == -2 & cancerTotal > 100], aes(Gene, cancerType, fill = percentage)) +
#     geom_tile(color = "white") +
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                          # midpoint = 50, limit = c(0,100), 
#                          space = "Lab", name = "2 copy loss%") +
#     theme_minimal()+ 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                      size = 12, hjust = 1))+
#     coord_fixed()
# delplot2 <- ggplot(data = lowCNAburden.count[value == -1 & cancerTotal > 100], aes(Gene, cancerType, fill = percentage)) +
#     geom_tile(color = "white") +
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                          # midpoint = 50, limit = c(0,100), 
#                          space = "Lab", name = "1 copy loss%") +
#     theme_minimal()+ 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                      size = 12, hjust = 1))+
#     coord_fixed()
delplot <- ggplot(data = lowCNAburden.count[CNAtype == "Del" & cancerTotal > 50], aes(Gene, cancerType, fill = percentage)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", high = "blue", 
                         space = "Lab", name = "% samples") +
    labs(title = "Deleted in sample with CNA burden < 20") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 10, hjust = 1))+
    coord_fixed()
ampplot <- ggplot(data = lowCNAburden.count[CNAtype == "Amp" & cancerTotal > 50], aes(Gene, cancerType, fill = percentage)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", high = "red", 
                         space = "Lab", name = "% samples") +
    labs(title = "Amplified in sample with CNA burden < 20") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 10, hjust = 1))+
    coord_fixed()
ggsave(filename = "ggplot.CNV.percentage.nohyperCNB.pdf", width = 20, height = 8, 
       plot = grid.arrange(delplot, ampplot, nrow = 1))
# sigcount1 <- as.data.table(table(lowCNAburden.count[percentage > 50 & CNAtype != "WT"]))
# p1 <- ggbarplot(sigcount, x = "Gene", y = "N", 
#                 color = "CNAtype", fill = "CNAtype", palette = "npg",
#                 position = position_dodge(0.6))
# ggpar(p1, x.text.angle = 45)
lowCNAburden.count$Gene <- factor(lowCNAburden.count$Gene, levels = sort(unique(lowCNAburden.count$Gene)))
p1 <- ggboxplot(data = lowCNAburden.count[cancerTotal > 50 & CNAtype != "WT"], x = "Gene", y = "percentage",
                color = "CNAtype", fill = "CNAtype", palette = "npg", alpha = 0.6, 
                outlier.shape = NA,
                label = "cancerType", repel = T, 
                label.select = list(criteria = "`percentage` > 60 & `cancerTotal` > 200"),
                add = "jitter", add.params = list(size = 0.6)) + geom_hline(yintercept = 60, linetype = 2)
p11 <- ggpar(p1, x.text.angle = 45)
ggsave(filename = "ggboxplot.CNV.percentage.nohyperCNB.pdf", width = 20, height = 8, 
       plot = p11)
# p2 <- ggdotplot(data = lowCNAburden.count[cancerTotal > 50 & CNAtype != "WT"], x = "Gene", y = "percentage",
#                 color = "CNAtype", fill = "CNAtype", palette = "npg", 
#                 binwidth = 0.5, size = 0.6, 
#                 add = "boxplot", #add.params = list(outlier.shape = NA),
#                 label = "cancerType", repel = T, 
#                 label.select = list(criteria = "`percentage` > 50 & `cancerTotal` > 400"))
# ggpar(p2, x.text.angle = 45)

# NUP93 CTCF co-delete --------------------------------------------------------------
library(survival)
library(survminer)
# mygenes.cnvs.2genes <- data.table(Sample = colnames(mygenes.cnvs)[-1],
#                                   NUP93 = as.vector(t(mygenes.cnvs[Sample == "NUP93",-1])), 
#                                   CTCF = as.vector(t(mygenes.cnvs[Sample == "CTCF",-1])),
#                                   NUP88 = as.vector(t(mygenes.cnvs[Sample == "NUP88",-1])))
# mygenes.cnvs.2genes[NUP93 < 0 & CTCF < 0]
mygenes.cnvs.t <- as.data.table(t(mygenes.cnvs[,-1]))
colnames(mygenes.cnvs.t) <- mygenes.cnvs$Sample
mygenes.cnvs.t <- data.table(Sample = colnames(mygenes.cnvs)[-1], mygenes.cnvs.t)
mygenes.cnvs.t[,Cancer := myinfo[match(mygenes.cnvs.t$Sample, Tumor_Sample_Barcode)]$`cancer type abbreviation`]
mygenes.cnvs.t[, OS := myinfo[match(mygenes.cnvs.t$Sample, Tumor_Sample_Barcode)]$OS]
mygenes.cnvs.t[, OS.time := myinfo[match(mygenes.cnvs.t$Sample, Tumor_Sample_Barcode)]$OS.time]
codel <- mygenes.cnvs.t[NUP93 < 0 & CTCF < 0]
cowt <- mygenes.cnvs.t[NUP93 == 0 & CTCF == 0]
coamp <- mygenes.cnvs.t[NUP93 > 0 & CTCF > 0]
# codel <- mygenes.cnvs.t[NUP93 < 0 & CTCF < 0 & NUP88 < 0]
# cowt <- mygenes.cnvs.t[NUP93 == 0 & CTCF == 0 & NUP88 == 0]
# coamp <- mygenes.cnvs.t[NUP93 > 0 & CTCF > 0 & NUP88 > 0]
cocnv <- rbindlist(list(codel, cowt, coamp))
cocnv[, Group := c(rep("co-del", nrow(codel)), rep("co-wt", nrow(cowt)), rep("co-amp", nrow(coamp)))]
cancers.count <- table(cocnv$Cancer)
cancers <- names(cancers.count > 100)
# cocnv.segs <- mysegs[sampleID %in% cocnv[Group!="co-wt"]$Sample]
# fwrite(cocnv.segs, file = "NUP93.CTCF.coCNV.seg", quote = F, sep = "\t")
for(i in 1:length(cancers)){
    mygenes.cnv.cancer <- cocnv[Cancer == cancers[i]]
    fit <- survfit(Surv(OS.time, OS) ~ Group, data = mygenes.cnv.cancer)
    p1 <- ggsurvplot(fit,
                     pval = TRUE, conf.int = FALSE,
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "strata", 
                     # surv.median.line = "hv",
                     # ggtheme = theme_bw(),
                     palette = "npg", title = paste("NUP93 & CTCF", cancers[i], sep = " in "))
    # pdf(file = paste("batch.ggsurvplot.coCNV", cancers[i], "pdf", sep = "."), width = 10, height = 8)
    # print(p1)
    # dev.off()
    png(filename = paste("batch.ggsurvplot.coCNV", cancers[i], "png", sep = "."), width = 1000, height = 800)
    print(p1)
    dev.off()
}

# mygenes.cnvs.t[, Status := paste(NUP93, CTCF, sep = "/")]
# chr16   56764017        56878797        ENSG00000102900.8|NUP93 .       +
# chr16   67596310        67673086        ENSG00000102974.10|CTCF .       +
# chr17   5264258 5323480 ENSG00000108559.7|NUP88 .       -
mygenes.cnvs.t[NUP93 < 0, NUP93.cnv := "Del"]
mygenes.cnvs.t[NUP93 == 0, NUP93.cnv := "WT"]
mygenes.cnvs.t[NUP93 > 0, NUP93.cnv := "Amp"]
mygenes.cnvs.t[CTCF < 0, CTCF.cnv := "Del"]
mygenes.cnvs.t[CTCF == 0, CTCF.cnv := "WT"]
mygenes.cnvs.t[CTCF > 0, CTCF.cnv := "Amp"]
mygenes.cnvs.t[, Status := paste(NUP93.cnv, CTCF.cnv, sep = "/")]
status.count <- sort(table(mygenes.cnvs.t$Status), decreasing = TRUE)
status.count2 <- as.data.table(table(mygenes.cnvs.t[,c(50,53)]))
status.count2[, sums := sum(N), by = Cancer]
status.count2[, Percentage := N/sums*100]
status.count2$Status <- factor(status.count2$Status, 
                               levels = c("Amp/Amp", "Amp/Del", "Amp/WT", "WT/Amp", "WT/WT", "WT/Del", "Del/Amp", "Del/Del", "Del/WT"))
p2 <- ggbarplot(data = status.count2[sums > 100], x = "Cancer", y = "Percentage", 
                color = "Status", fill = "Status", palette = "jco",
                width = 0.9, xlab = "", ylab = "") + 
    geom_text(aes(x = Cancer, y = 102, label = sums))
p2
ggsave(filename = "ggbarplot.NUP93.CTCF.cocnv.pdf", width = 15, height = 6, p2)

# CNV & expression --------------------------------------------------------------
# mygenes.cnvs <- mycnvs[Sample %in% mygenes$targets]
# mygenes.cnvs.mt <- melt.data.table(mygenes.cnvs, id.vars = "Sample")
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
mygenes.exprs.cnv.mut[MutationLoad > hypermutatecut, hypermutation := "HM"]
mygenes.exprs.cnv.mut[is.na(hypermutation), hypermutation := "NotHM"]

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

##### NUP93-/CTCF- NUP93-/CTCFwt NUP93wt/CTCF- NUP93wt/CTCFwt normal
library(survival)
library(survminer)
myinfo.sim <- myinfo[,c(1,2,3,7,26,27)]
ctcf.tmp <- mygenes.exprs.cnv.mut[Gene == "CTCF"]
nup93.tmp <- mygenes.exprs.cnv.mut[Gene == "NUP93"]
cnaburden.tmp <- unique(mygenes.cnvs.mt[,c(2,4)])
myinfo.sim[, CTCF.CNV := ctcf.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$CNVStatus]
myinfo.sim[, CTCF.mut := ctcf.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$MutationStatus]
myinfo.sim[, CTCF.expr := ctcf.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$RNA]
myinfo.sim[, NUP93.CNV := nup93.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$CNVStatus]
myinfo.sim[, NUP93.mut := nup93.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$MutationStatus]
myinfo.sim[, NUP93.expr := nup93.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$RNA]
myinfo.sim[, STcode := substr(Tumor_Sample_Barcode, 14,15)]
myinfo.sim[STcode == "01", Type := "Primary Solid Tumor"]
myinfo.sim[STcode == "02", Type := "Recurrent Solid Tumor"]
myinfo.sim[STcode == "03", Type := "Primary Blood Derived Cancer"]
myinfo.sim[STcode == "05", Type := "Additional - New Primary"]
myinfo.sim[STcode == "06", Type := "Metastatic"]
myinfo.sim[STcode == "07", Type := "Additional Metastatic"]
myinfo.sim[STcode == "11", Type := "Solid Tissue Normal"]
myinfo.sim[!is.na(CTCF.CNV), CNVgroup := paste(CTCF.CNV, NUP93.CNV, sep = "/")]
myinfo.sim[, CNAburden := cnaburden.tmp[match(myinfo.sim$Tumor_Sample_Barcode, variable)]$CNAburden]
myinfo.sim <- data.table(myinfo.sim, myinfo[,c(28:33)])
table(myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11"]$CNVgroup)
table(myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11" & CNAburden < 20]$CNVgroup)
table(myinfo.sim[STcode != "11" & CNAburden < 20, c(3,15)])
myinfo.sim.brca <- myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11" & CNAburden < 20]
# myinfo.sim.brca <- myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11"]
# myinfo.sim.brca.1 <- myinfo.sim.brca[CNVgroup %in% c("Amp/Amp", "Del/Del", "WT/WT")]
myinfo.sim.brca.1 <- myinfo.sim.brca[STcode == "01" & CNVgroup %in% c("Del/Del", "WT/WT")]
fit <- survfit(Surv(OS.time, OS) ~ CNVgroup, data = myinfo.sim.brca.1)
# fit <- survfit(Surv(OS.time) ~ CNVgroup, data = myinfo.sim.brca.1)
# fit <- survfit(Surv(PFI.time, PFI) ~ CNVgroup, data = myinfo.sim.brca.1)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 # surv.median.line = "hv",
                 # ggtheme = theme_bw(),
                 palette = "npg",
                 # palette = get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 28),
                 title = "BRCA CTCF/NUP93")
pdf(file = "ggsurvplot.CNV.codel.BRCA.pdf", width = 10, height = 8)
print(p1)
dev.off()

myinfo.sim.kirp <- myinfo.sim[`cancer type abbreviation` == "KIRP" & STcode != "11" & CNAburden < 20 & CNVgroup %in% c("Amp/Amp", "Del/Del", "WT/WT")]
fit <- survfit(Surv(OS.time, OS) ~ CNVgroup, data = myinfo.sim.kirp)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 palette = "npg",
                 title = "KIRP CTCF/NUP93")
pdf(file = "ggsurvplot.CNV.codel.KIRP.pdf", width = 10, height = 8)
print(p1)
dev.off()

myinfo.sim.luad <- myinfo.sim[`cancer type abbreviation` == "LUAD" & STcode != "11" & CNAburden < 20 & CNVgroup %in% c("Amp/Amp", "Del/Del", "WT/WT")]
fit <- survfit(Surv(OS.time, OS) ~ CNVgroup, data = myinfo.sim.luad)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 palette = "npg",
                 title = "LUAD CTCF/NUP93")
pdf(file = "ggsurvplot.CNV.codel.LUAD.pdf", width = 10, height = 8)
print(p1)
dev.off()

myinfo.sim.lusc <- myinfo.sim[`cancer type abbreviation` == "LUSC" & STcode != "11" & CNAburden < 20 & CNVgroup %in% c("Amp/Amp", "Del/Del", "WT/WT")]
fit <- survfit(Surv(OS.time, OS) ~ CNVgroup, data = myinfo.sim.lusc)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 palette = "npg",
                 title = "LUSC CTCF/NUP93")
pdf(file = "ggsurvplot.CNV.codel.LUSC.pdf", width = 10, height = 8)
print(p1)
dev.off()

####
myinfo.sim.brca.expr.n <- myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode == "11" & !is.na(NUP93.expr)]
myinfo.sim.brca.expr.t <- myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11" & CNAburden < 20 & !is.na(NUP93.expr)]
myinfo.sim.brca.expr.nt <- myinfo.sim[`cancer type abbreviation` == "BRCA" & 
                                          `_PATIENT` %in% myinfo.sim.brca.expr.n$`_PATIENT` & 
                                          `_PATIENT` %in% myinfo.sim.brca.expr.t$`_PATIENT`]

##### expression
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
                title = "NUP93 mutation/CNV in BRCA (with both mutation and CNV available)",
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

# NUP93 E14K mutation --------------------------------------------------------------
NUP93.E14K.mut <- mymuts[Hugo_Symbol == "NUP93" & Amino_Acid_Change == "p.E14K"]
NUP93.E14K.dt <- mygenes.exprs.cnv.mut[Tumor_Sample_Barcode %in% NUP93.E14K.mut$Tumor_Sample_Barcode & Gene == "NUP93"]
mygenes.exprs.cnv.mut[PATIENT %in% NUP93.E14K.dt$PATIENT & Gene == "NUP93"]
# NUP93.E14K.NT.paired <- "TCGA-BH-A209"
patients.ids <- which(names(myexprs) %in% c("TCGA-BH-A209-01", "TCGA-BH-A209-11"))
patients.table <- myexprs[, .SD, .SDcols = c(1, patients.ids)]
fwrite(patients.table, file = "NUP93.E14K.TCGA-BH-A209.NT.expr.tsv", sep = "\t")
# NUP93 E14K mutation co-delete --------------------------------------------------------------
library(survival)
library(survminer)
# myexprs.t <- as.data.table(t(myexprs[,-1]))
# colnames(myexprs.t) <- myexprs$sample
# myexprs.t <- data.table(Sample = colnames(myexprs)[-1], myexprs.t)
# myexprs.t[,Cancer := myinfo[match(myexprs.t$Sample, Tumor_Sample_Barcode)]$`cancer type abbreviation`]
# myexprs.t[, OS := myinfo[match(myexprs.t$Sample, Tumor_Sample_Barcode)]$OS]
# myexprs.t[, OS.time := myinfo[match(myexprs.t$Sample, Tumor_Sample_Barcode)]$OS.time]
mcf7 <- fread("NUP93.E14K/MCF.E2.E14K.RNAseq", header = T)
patients <- c("TCGA-BH-A209-01","TCGA-BH-A209-11")
patients.muts <- mymuts[Tumor_Sample_Barcode %in% patients]
patients.muts.count <- as.data.table(table(patients.muts[,c(1,9)]))
patients.cnvs <- mycnvs[, .SD, .SDcols = c(1, which(names(mycnvs) %in% patients))]
patients.ids <- which(names(myexprs) %in% patients)
patients.table <- myexprs[, .SD, .SDcols = c(1, patients.ids)]
patients.table[, MutationStatus := patients.muts.count[match(patients.table$sample, Hugo_Symbol)]$N]
patients.table[,CNVStatus := patients.cnvs[match(patients.table$sample, Sample)]$`TCGA-BH-A209-01`]
colnames(patients.table)[2:3] <- c("Tumor","Normal")
patients.table[!is.na(MutationStatus), Mutation := "Mut"]
patients.table[is.na(MutationStatus), Mutation := "WT"]
patients.table[CNVStatus < 0, CNV := "Del"]
patients.table[CNVStatus == 0, CNV := "WT"]
patients.table[CNVStatus > 0, CNV := "Amp"]
patients.table[is.na(CNVStatus), CNV := "Unkown"]
patients.table$CNV <- factor(patients.table$CNV, levels = c("Del", "WT", "Amp", "Unkown"))
patients.table$Mutation <- factor(patients.table$Mutation, levels = c("Mut", "WT"))
patients.table[, Log2FC := Tumor - Normal]
patients.table[, E2.expr := mcf7[match(patients.table$sample, gene)]$E2.expr]
p1 <- ggscatter(patients.table, x = "Normal", y = "Tumor",
                color = "Mutation", fill = "Mutation",
                shape = "CNV", size = 0.9, 
                label = "sample", repel = TRUE, font.label = c(8, "plain"),
                label.select = list(criteria = "(`Log2FC` > 5 & `Tumor` > 10) | (`Log2FC` < -5 & `Normal` > 10)"),
                title = "expresion for BRCA TCGA-BH-A209 patient with NUP93 p.E14K mutation & NUP93/CTCF deletion")
# p1
ggsave(filename = "NUP93.E14K/TCGA-BH-A209.NT.RNAseq.pdf", plot = p1, width = 12, height = 12)
p1 <- ggscatter(patients.table, x = "Normal", y = "Tumor",
                color = "Mutation", fill = "Mutation",
                shape = "CNV", size = 0.9, 
                label = "sample", repel = TRUE, font.label = c(12, "plain","blue"),
                label.select = list(criteria = "(`Log2FC` > 0.58 & `E2.expr` == 'up') | (`Log2FC` < -0.58 & `E2.expr` == 'down')"),
                title = "expresion for BRCA TCGA-BH-A209 patient with NUP93 p.E14K mutation & NUP93/CTCF deletion")
# p1
ggsave(filename = "NUP93.E14K/TCGA-BH-A209.NT.RNAseq.cmp2MCF7mut.pdf", plot = p1, width = 12, height = 12)

# expression --------------------------------------------------------------
myexprs.sample <- data.table(sampleID = colnames(myexprs), 
                             cancerType = myinfo[match(colnames(myexprs), Tumor_Sample_Barcode)]$`cancer type abbreviation`)
myexprs.sample[, STcode := substr(sampleID, 14,15)]
myexprs.sample[STcode == "01", Type := "Primary Solid Tumor"]
myexprs.sample[STcode == "02", Type := "Recurrent Solid Tumor"]
myexprs.sample[STcode == "03", Type := "Primary Blood Derived Cancer"]
myexprs.sample[STcode == "05", Type := "Additional - New Primary"]
myexprs.sample[STcode == "06", Type := "Metastatic"]
myexprs.sample[STcode == "07", Type := "Additional Metastatic"]
myexprs.sample[STcode == "11", Type := "Solid Tissue Normal"]
myexprs.sample[, Patient := myinfo[match(myexprs.sample$sampleID, Tumor_Sample_Barcode)]$`_PATIENT`]
myexprs.sample.n <- myexprs.sample[STcode == "11" & !is.na(cancerType)]$Patient
myexprs.sample.t <- myexprs.sample[STcode != "11" & !is.na(cancerType)]$Patient
myexprs.sample.nt <- intersect(myexprs.sample.n, myexprs.sample.t)
mycolidx <- which(myexprs.sample$Patient %in% myexprs.sample.nt)
myexprs.sample.nt.mat <- data.matrix(myexprs[sample %in% mygenes$targets, .SD, .SDcols = mycolidx])
rownames(myexprs.sample.nt.mat) <- myexprs[sample %in% mygenes$targets]$sample
mycolidx1 <- which(myexprs.sample$Patient %in% myexprs.sample.nt & myexprs.sample$STcode == "11")
mycolidx2 <- which(myexprs.sample$Patient %in% myexprs.sample.nt & myexprs.sample$STcode != "11")
mydt1 <- myexprs[sample %in% mygenes$targets, .SD, .SDcols = mycolidx1]
mydt2 <- myexprs[sample %in% mygenes$targets, .SD, .SDcols = mycolidx2]
myexprs.sample.nt.mat.ordered <- data.matrix(data.table(mydt1, mydt2))
rownames(myexprs.sample.nt.mat.ordered) <- myexprs[sample %in% mygenes$targets]$sample
row.annos <- data.frame(geneCoords = mygenecoords[match(rownames(myexprs.sample.nt.mat), Gene)]$V1,
                        geneClass = mygenes[match(rownames(myexprs.sample.nt.mat), targets)]$Class)
rownames(row.annos) <- rownames(myexprs.sample.nt.mat)
col.annos <- data.frame(Tissues = myexprs.sample[match(colnames(myexprs.sample.nt.mat), sampleID)]$Type, 
                        cancerType = myexprs.sample[match(colnames(myexprs.sample.nt.mat), sampleID)]$cancerType)
rownames(col.annos) <- colnames(myexprs.sample.nt.mat)
mycols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
a1 <- colorRampPalette(brewer.pal(n = 5, name = "PuBuGn"))(length(unique(row.annos$geneClass)))
names(a1) <- sort(levels(row.annos$geneClass))
a2 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(row.annos$geneCoords)))
# names(a2) <- sort(levels(row.annos$geneCoords))
names(a2) <- paste("chr", c(seq(1,13,1), seq(16,20,1),"22","X"), sep = "")
a3 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(col.annos$cancerType)))
names(a3) <- levels(col.annos$cancerType)
a4 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(col.annos$Tissues)))
names(a4) <- levels(col.annos$Tissues)
ann_colors = list(
    geneClass = a1,
    geneCoords = a2,
    cancerType = a3, 
    Tissues = a4
)
pdf(file = "pheatmap.NT.expression.pdf", width = 15, height = 10)
pheatmap(myexprs.sample.nt.mat, annotation_row = row.annos, annotation_col = col.annos, 
         show_colnames = F, color = mycols, fontsize = 7, fontsize_row = 10, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D",
         annotation_colors = ann_colors, 
         main = "NUPs/Cohesin/Condensin expression levels in 1427 (712 N-T paired) samples")
dev.off()

mybreaks <- c(seq(-10, -3, 0.5), seq(-2.9, 2.9, 0.1), seq(3, 10, 0.5))
# mycols <- colorRampPalette(rev(brewer.pal(n = 11, name = "BrBG")))(length(mybreaks))
mycols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(mybreaks))
pdf(file = "pheatmap.NT.expression.scalebyrow.pdf", width = 15, height = 10)
pheatmap(myexprs.sample.nt.mat, annotation_row = row.annos, annotation_col = col.annos, 
         show_colnames = F, fontsize = 7, fontsize_row = 10, 
         scale = "row",color = mycols,  breaks = mybreaks,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_colors = ann_colors, 
         main = "NUPs/Cohesin/Condensin expression levels in 1427 (712 N-T paired) samples")
dev.off()

col.annos <- data.frame(Tissues = myexprs.sample[match(colnames(myexprs.sample.nt.mat.ordered), sampleID)]$Type, 
                        cancerType = myexprs.sample[match(colnames(myexprs.sample.nt.mat.ordered), sampleID)]$cancerType)
rownames(col.annos) <- colnames(myexprs.sample.nt.mat.ordered)
row.annos <- data.frame(geneCoords = mygenecoords[match(rownames(myexprs.sample.nt.mat.ordered), Gene)]$V1,
                        geneClass = mygenes[match(rownames(myexprs.sample.nt.mat.ordered), targets)]$Class)
rownames(row.annos) <- rownames(myexprs.sample.nt.mat.ordered)
a1 <- colorRampPalette(brewer.pal(n = 5, name = "PuBuGn"))(length(unique(row.annos$geneClass)))
names(a1) <- sort(levels(row.annos$geneClass))
a2 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(row.annos$geneCoords)))
# names(a2) <- sort(levels(row.annos$geneCoords))
names(a2) <- paste("chr", c(seq(1,13,1), seq(16,20,1),"22","X"), sep = "")
a3 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(col.annos$cancerType)))
names(a3) <- levels(col.annos$cancerType)
a4 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(col.annos$Tissues)))
names(a4) <- levels(col.annos$Tissues)
ann_colors = list(
    geneClass = a1,
    geneCoords = a2,
    cancerType = a3, 
    Tissues = a4
)
mybreaks <- c(seq(-10, -2, 0.5), seq(-1.9, 1.9, 0.1), seq(2, 10, 0.5))
mycols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(mybreaks))
pdf(file = "pheatmap.NT.expression.scalebyrow.ordered.pdf", width = 11, height = 8)
pheatmap(myexprs.sample.nt.mat.ordered, annotation_row = row.annos, annotation_col = col.annos, 
         show_colnames = F, show_rownames = T, fontsize = 6, fontsize_row = 10, cluster_cols = F, 
         scale = "row",color = mycols,  breaks = mybreaks,
         clustering_distance_rows = "euclidean",
         # clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_colors = ann_colors, 
         gaps_col = c(712), #gaps_row = c(0,16,35), #kmeans_k = 3,
         main = "NUPs/Cohesin/Condensin expression levels in 1427 (712 N-T paired) samples")
dev.off()

myexprs.sample.nt.dt <- melt.data.table(myexprs[sample %in% mygenes$targets, .SD, .SDcols = c(1,mycolidx)],
                                        id.vars = "sample")
myexprs.sample.nt.dt[, STcode := substr(variable, 14,15)]
myexprs.sample.nt.dt[STcode == "01", Type := "Primary Solid Tumor"]
myexprs.sample.nt.dt[STcode == "02", Type := "Recurrent Solid Tumor"]
myexprs.sample.nt.dt[STcode == "03", Type := "Primary Blood Derived Cancer"]
myexprs.sample.nt.dt[STcode == "05", Type := "Additional - New Primary"]
myexprs.sample.nt.dt[STcode == "06", Type := "Metastatic"]
myexprs.sample.nt.dt[STcode == "07", Type := "Additional Metastatic"]
myexprs.sample.nt.dt[STcode == "11", Type := "Solid Tissue Normal"]
myexprs.sample.nt.dt[, cancerType := myinfo[match(myexprs.sample.nt.dt$variable, Tumor_Sample_Barcode)]$`cancer type abbreviation`]
myexprs.sample.nt.dt[, meanExpr := mean(value), by = c("sample", "STcode", "cancerType")]
Ncounts <- table(myexprs.sample.nt.dt[STcode == "11" & sample == "NUP93"]$cancerType)
Cancers <- unique(myexprs.sample.nt.dt[cancerType %in% names(Ncounts[Ncounts > 10])]$cancerType)
Genes <- unique(myexprs.sample.nt.dt$sample)
length(Genes)
# no WAPL
pvals <- matrix(nrow = length(Genes), ncol = length(Cancers))
for(i in 1:length(Genes)){
    for(j in 1:length(Cancers)){
        pvals[i, j ] <- t.test(myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "01"]$value,
                               myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "11"]$value)$p.value
    }
}
rownames(pvals) <- Genes
colnames(pvals) <- Cancers
sigcounts <- apply(pvals, MARGIN = 1, function(x){length(x[x<1e-5])})
sort(sigcounts, decreasing = T)
# NCAPD2  NCAPH  NCAPG NUP107  NUP62  NUP85  LMNB2 NCAPG2  LMNB1 NCAPD3 NUP205  NUP93  NUP37   RAE1 NUP155 NUP188  NUP88 POM121 
# 12     12     11     11     11     11     10     10      9      8      8      8      7      7      6      6      6      6 
# SMC2   SMC4 NUP210  NUP35  NUP50  PDS5B  RAD21  SMC1B   SMC6   LMNA NUP133  NUP54  SEC13  SMC1A NUP160  NUP43  SEH1L NCAPH2 
# 6      6      5      5      5      5      5      5      5      4      4      4      4      4      3      3      3      2 
# NUP153 NUP214  STAG3   CTCF  PDS5A   SMC3  STAG1    TPR  NUP98   SMC5  STAG2 
# 2      2      2      1      1      1      1      1      0      0      0 
sigcounts.dt <- data.table(Gene = names(sigcounts), ttest.p = sigcounts, 
                           Class = mygenes[match(names(sigcounts), targets)]$Class)
p0 <- ggdotchart(sigcounts.dt, x = "Gene", y = "ttest.p", size = 2,
                 color = "Class", add = "segment", sorting = "descending", rotate = T,
                 ylab = "Count of Cancners with t test p value < 1e-5",
                 title = "15 Cancers with Normal samples > 10")
ggsave(filename = "ggdotchart.ttestsig.pdf", plot = p0, width = 8, height = 8)
pdf(file = "pheatmap.NT.expression.ttest.pdf", width = 15, height = 10)
pheatmap(-10 * log10(pvals), breaks = seq(0,20,0.5), 
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(41),
         show_colnames = T, show_rownames = T, fontsize = 7, fontsize_row = 10, cluster_cols = F, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         main = "NUPs/Cohesin/Condensin expression levels in 1427 (712 N-T paired) samples")
dev.off()

tmp <- myexprs.sample.nt.dt[STcode %in% c("01", "11") & sample %in% c("NCAPD2","NCAPH","NUP93","SEC13","TPR","STAG2") & cancerType %in% names(Ncounts[Ncounts > 10])]
tmp$sample <- factor(tmp$sample, levels = c("NCAPH","NCAPD2","NUP93","SEC13","TPR","STAG2"))
tmp$cancerType <- factor(tmp$cancerType, levels = names(sort(pvals[rownames(pvals)=="NUP93",])))
p1 <- ggboxplot(tmp,
                x = "cancerType", y = "value", 
                color = "Type", fill = "Type", palette = "npg", alpha = 0.6,
                xlab = "", ylab = "log2(TPM+1)",
                add = "jitter", add.params = list(size = 0.2), outlier.shape = NA)
p11 <- ggpar(facet(p1, facet.by = "sample", nrow = 3, scales = "free_y"), x.text.angle = 30)
ggsave(filename = "ggboxplot.NTpaired.NUP93.withConrol.pdf", plot = p11, width = 12, height = 8)

##### test
log2FC <- matrix(nrow = length(Genes), ncol = length(Cancers))
for(i in 1:length(Genes)){
    for(j in 1:length(Cancers)){
        log2FC[i, j ] <- mean(myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "01"]$value, na.rm = T) - 
            mean(myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "11"]$value, na.rm = T)
    }
}
rownames(log2FC) <- Genes
colnames(log2FC) <- Cancers

t.test.pvals <- matrix(nrow = length(Genes), ncol = length(Cancers))
for(i in 1:length(Genes)){
    for(j in 1:length(Cancers)){
        t.test.pvals[i, j ] <- t.test(myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "01"]$value,
                                      myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "11"]$value)$p.value
    }
}
rownames(t.test.pvals) <- Genes
colnames(t.test.pvals) <- Cancers
sigcounts.t.test.pvals <- apply(t.test.pvals, MARGIN = 1, function(x){length(x[x<1e-3])})
sort(sigcounts.t.test.pvals, decreasing = T)

wilcox.test.pvals <- matrix(nrow = length(Genes), ncol = length(Cancers))
for(i in 1:length(Genes)){
    for(j in 1:length(Cancers)){
        wilcox.test.pvals[i, j ] <- wilcox.test(myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "01"]$value,
                                                myexprs.sample.nt.dt[sample == Genes[i] & cancerType == Cancers[j] & STcode == "11"]$value)$p.value
    }
}
rownames(wilcox.test.pvals) <- Genes
colnames(wilcox.test.pvals) <- Cancers
sigcounts.wilcox.test.pvals <- apply(wilcox.test.pvals, MARGIN = 1, function(x){length(x[x<1e-3])})
sort(sigcounts.wilcox.test.pvals, decreasing = T)

t.test.pvals.fc <- -1 * log10(t.test.pvals) * sign(log2FC)
wilcox.test.pvals.fc <- -1 * log10(wilcox.test.pvals) * sign(log2FC)

row.annos <- data.frame(geneCoords = mygenecoords[match(rownames(wilcox.test.pvals.fc), Gene)]$V1,
                        geneClass = mygenes[match(rownames(wilcox.test.pvals.fc), targets)]$Class)
rownames(row.annos) <- rownames(wilcox.test.pvals.fc)
a1 <- colorRampPalette(brewer.pal(n = 5, name = "PuBuGn"))(length(unique(row.annos$geneClass)))
names(a1) <- sort(levels(row.annos$geneClass))
a2 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(row.annos$geneCoords)))
# names(a2) <- sort(levels(row.annos$geneCoords))
names(a2) <- paste("chr", c(seq(1,13,1), seq(16,20,1),"22","X"), sep = "")
ann_colors = list(
    geneClass = a1,
    geneCoords = a2
)
# mybreaks <- c(seq(-6, -3, 1), seq(-2.9, 2.9, 0.1), seq(3, 6, 1))
mybreaks <- seq(-5,5,1)
mycols <- colorRampPalette(c("blue", "white", "red"))(length(mybreaks)-1)
pdf(file = "pheatmap.NT.expression.wilcox.test.pvals.fc.pdf", width = 12, height = 8)
pheatmap(wilcox.test.pvals.fc, annotation_row = row.annos, 
         show_colnames = T, fontsize = 10, fontsize_row = 10, 
         scale = "none", color = mycols, breaks = mybreaks,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D",
         display_numbers = T,
         number_format = "%.1f",
         annotation_colors = ann_colors, 
         main = "wilcox.test.pvals of NUPs/Cohesin/Condensin expression in 1427 (712 N-T paired) samples")
dev.off()
log2FC.1 <- log2FC
log2FC.1[wilcox.test.pvals > 1e-2] <- 0
mybreaks <- seq(-5,5,1)
mycols <- colorRampPalette(c("blue", "white", "red"))(length(mybreaks)-1)
pdf(file = "pheatmap.NT.expression.wilcox.test.pvals.fc.1.pdf", width = 12, height = 8)
pheatmap(log2FC.1, annotation_row = row.annos, 
         show_colnames = T, fontsize = 10, fontsize_row = 10, 
         scale = "none", color = mycols, breaks = mybreaks,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D",
         display_numbers = T,
         number_format = "%.1f",
         annotation_colors = ann_colors, 
         main = "log2FC of NUPs/Cohesin/Condensin expression in 1427 (712 N-T paired) samples (p>=0.01 were labeled as 0)")
dev.off()

wilcox.test.pvals.fc.dt <- melt.data.table(data.table(Genes, wilcox.test.pvals.fc), id.vars = c("Genes"))
wilcox.test.pvals.fc.dt[value > 2, Type := "Up"]
wilcox.test.pvals.fc.dt[value < -2, Type := "Down"]
wilcox.test.pvals.fc.dt[is.na(Type), Type := "NoChange"]
wilcox.test.pvals.fc.dt.count <- as.data.table(table(wilcox.test.pvals.fc.dt[,c(1,4)]))
wilcox.test.pvals.fc.dt.count.1 <- dcast.data.table(wilcox.test.pvals.fc.dt.count, formula = Genes ~ Type)
setorder(wilcox.test.pvals.fc.dt.count.1, -Up, -Down, NoChange, -Genes)
setorder(wilcox.test.pvals.fc.dt.count, Type, -N)

wilcox.test.pvals.fc.dt.count$Type <- factor(wilcox.test.pvals.fc.dt.count$Type, levels = c("Up", "Down", "NoChange"))
wilcox.test.pvals.fc.dt.count$Genes <- factor(wilcox.test.pvals.fc.dt.count$Genes, levels = wilcox.test.pvals.fc.dt.count.1$Genes)
p1 <- ggbarplot(data = wilcox.test.pvals.fc.dt.count, x = "Genes", y = "N",
                color = "Type", fill = "Type", palette = "npg", 
                xlab = "", ylab = "Count of Cancer type",
                title = "wilcox.test of NUPs/Cohesin/Condensin expression in 1427 (712 N-T paired) samples (p<0.01)")
p11 <- ggpar(p1, x.text.angle = 45, font.x = 10)
ggsave(filename = "ggbarplot.wilcox.test.pvals.fc.count.pdf", width = 12, height = 6, p11)

# expr & survival --------------------------------------------------------------
myinfo.sim <- myinfo[,c(1,2,3,7,26,27)]
ctcf.tmp <- mygenes.exprs.cnv.mut[Gene == "CTCF" & STcode != "11"]
nup93.tmp <- mygenes.exprs.cnv.mut[Gene == "NUP93" & STcode != "11"]
ncapg.tmp <- mygenes.exprs.cnv.mut[Gene == "NCAPG" & STcode != "11"]
ctcf.tmp[, medianExpr := median(RNA), by = "Cancer"]
nup93.tmp[, medianExpr := median(RNA), by = "Cancer"]
ncapg.tmp[, medianExpr := median(RNA), by = "Cancer"]
nup93.tmp[, medianExpr := median(RNA), by = "Cancer"]
ctcf.tmp[RNA > medianExpr, RNA.group := "High"]
ctcf.tmp[RNA <= medianExpr, RNA.group := "Low"]
nup93.tmp[RNA > medianExpr, RNA.group := "High"]
nup93.tmp[RNA <= medianExpr, RNA.group := "Low"]
ncapg.tmp[RNA > medianExpr, RNA.group := "High"]
ncapg.tmp[RNA <= medianExpr, RNA.group := "Low"]
# ctcf.tmp <- mygenes.exprs.cnv.mut[Gene == "CTCF" & STcode != "11"]
# nup93.tmp <- mygenes.exprs.cnv.mut[Gene == "NUP93" & STcode != "11"]
# ctcf.tmp[, highcut := quantile(RNA, probs = 0.66), by = "Cancer"]
# nup93.tmp[, highcut := quantile(RNA, probs = 0.66), by = "Cancer"]
# ctcf.tmp[, lowcut := quantile(RNA, probs = 0.33), by = "Cancer"]
# nup93.tmp[, lowcut := quantile(RNA, probs = 0.33), by = "Cancer"]
# ctcf.tmp[RNA > highcut, RNA.group := "High"]
# ctcf.tmp[RNA < lowcut, RNA.group := "Low"]
# nup93.tmp[RNA > highcut, RNA.group := "High"]
# nup93.tmp[RNA < lowcut, RNA.group := "Low"]
# ctcf.tmp[is.na(RNA.group), RNA.group := "Median"]
# nup93.tmp[is.na(RNA.group), RNA.group := "Median"]
cnaburden.tmp <- unique(mygenes.cnvs.mt[,c(2,4)])
myinfo.sim[, CTCF.CNV := ctcf.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$CNVStatus]
myinfo.sim[, CTCF.mut := ctcf.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$MutationStatus]
myinfo.sim[, CTCF.expr := ctcf.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$RNA.group]
myinfo.sim[, NUP93.CNV := nup93.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$CNVStatus]
myinfo.sim[, NUP93.mut := nup93.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$MutationStatus]
myinfo.sim[, NUP93.expr := nup93.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$RNA.group]
myinfo.sim[, NCAPG.CNV := ncapg.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$CNVStatus]
myinfo.sim[, NCAPG.mut := ncapg.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$MutationStatus]
myinfo.sim[, NCAPG.expr := ncapg.tmp[match(myinfo.sim$Tumor_Sample_Barcode, Tumor_Sample_Barcode)]$RNA.group]
myinfo.sim[, STcode := substr(Tumor_Sample_Barcode, 14,15)]
myinfo.sim[STcode == "01", Type := "Primary Solid Tumor"]
myinfo.sim[STcode == "02", Type := "Recurrent Solid Tumor"]
myinfo.sim[STcode == "03", Type := "Primary Blood Derived Cancer"]
myinfo.sim[STcode == "05", Type := "Additional - New Primary"]
myinfo.sim[STcode == "06", Type := "Metastatic"]
myinfo.sim[STcode == "07", Type := "Additional Metastatic"]
myinfo.sim[STcode == "11", Type := "Solid Tissue Normal"]
myinfo.sim[!is.na(CTCF.CNV), CNVgroup := paste(CTCF.CNV, NUP93.CNV, sep = "/")]
myinfo.sim[, CNAburden := cnaburden.tmp[match(myinfo.sim$Tumor_Sample_Barcode, variable)]$CNAburden]
myinfo.sim <- data.table(myinfo.sim, myinfo[,c(28:33)])
colnames(myinfo.sim)[3] <- "CancerType"
fit <- survfit(Surv(OS.time, OS) ~ NCAPG.expr, data = myinfo.sim)
p1 <- ggsurvplot_facet(fit, data = myinfo.sim, facet.by = "CancerType",
                       pval = TRUE, conf.int = FALSE, 
                       scales = "free_x", nrow = 5,
                       xlab = "OSS (Days)",
                       palette = "npg", title = "NCAPG expression")
pdf(file = "ggsurvplot_facet.all.by.NCAPG.expression.OS.pdf", width = 12, height = 15)
p1
dev.off()
fit <- survfit(Surv(OS.time, OS) ~ NUP93.expr, data = myinfo.sim)
p1 <- ggsurvplot_facet(fit, data = myinfo.sim, facet.by = "CancerType",
                       pval = TRUE, conf.int = FALSE, 
                       scales = "free_x", nrow = 5,
                       xlab = "OSS (Days)",
                       palette = "npg", title = "NUP93 expression")
pdf(file = "ggsurvplot_facet.all.by.NUP93.expression.OS.pdf", width = 12, height = 15)
p1
dev.off()
fit <- survfit(Surv(OS.time, OS) ~ CTCF.expr, data = myinfo.sim)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 palette = "npg",
                 # palette = get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 28),
                 title = "pancancer CTCF expression")
# fit <- survfit(Surv(OS.time, OS) ~ NUP93.expr + NUP93.mut, data = myinfo.sim.1)
# p1 <- ggsurvplot(fit,
#                  pval = TRUE, conf.int = FALSE,
#                  risk.table = TRUE, 
#                  risk.table.col = "strata", 
#                  linetype = "strata", 
#                  palette = "npg",
#                  title = "pancancer NUP93 expression")
# # ggsurvplot_group_by(fit = fit, data = myinfo.sim, group.by = "CancerType")
# p1 <- ggsurvplot_facet(fit, data = myinfo.sim.1, facet.by = "CancerType",
#                        pval = TRUE, conf.int = FALSE, 
#                        scales = "free_x", nrow = 3,
#                        palette = "npg", title = "CTCF")
# pdf(file = paste("ggsurvplot.grouped.expression", "CTCF", "pdf", sep = "."), width = 33, height = 9)
# print(p1)
# dev.off()
myinfo.sim.1 <- myinfo.sim[NUP93.mut == "WT"]
fit <- survfit(Surv(OS.time, OS) ~ NUP93.expr, data = myinfo.sim.1)
p1 <- ggsurvplot_facet(fit, data = myinfo.sim, facet.by = "CancerType",
                       pval = TRUE, conf.int = FALSE, 
                       scales = "free_x", nrow = 3,
                       palette = "npg", title = "NUP93 expression (WT)")
pdf(file = paste("ggsurvplot.grouped.expression", "NUP93", "pdf", sep = "."), width = 33, height = 9)
print(p1)
dev.off()

a <- sort(table(myinfo.sim[STcode != "11" & CNAburden < 20 & Tumor_Sample_Barcode %in% mymuts.samples[hypermutation=="NotHM"]$Tumor_Sample_Barcode]$CancerType), decreasing = T)
myinfo.sim.largecohort <- myinfo.sim[STcode != "11" & 
                                         Tumor_Sample_Barcode %in% mymuts.samples[hypermutation=="NotHM"]$Tumor_Sample_Barcode & 
                                         CNAburden < 20 &
                                         CancerType %in% names(a[a>100]) & NUP93.mut == "WT"]
fit <- survfit(Surv(OS.time, OS) ~ NUP93.expr, data = myinfo.sim.largecohort)
p1 <- ggsurvplot_facet(fit, data = myinfo.sim.largecohort, facet.by = "CancerType",
                       pval = TRUE, conf.int = FALSE, 
                       scales = "free_x", nrow = 3,
                       palette = "npg", title = "NUP93 expression (WT)")
pdf(file = paste("ggsurvplot.grouped.expression.noHM.nohyperCNB", "NUP93", "pdf", sep = "."), width = 33, height = 9)
print(p1)
dev.off()
fit <- survfit(Surv(OS.time, OS) ~ NUP93.expr + CTCF.expr, data = myinfo.sim.largecohort)
p1 <- ggsurvplot_facet(fit, data = myinfo.sim.largecohort, facet.by = "CancerType",
                       pval = TRUE, conf.int = FALSE, 
                       scales = "free_x", nrow = 3,
                       palette = "npg", title = "NUP93 expression (WT)")
pdf(file = paste("ggsurvplot.grouped.expression.noHM.nohyperCNB", "NUP93.CTCF", "pdf", sep = "."), width = 33, height = 9)
print(p1)
dev.off()

myinfo.sim.largecohort.1 <- myinfo.sim.largecohort[CancerType %in% c("KIRC", "GBM")]
fit <- survfit(Surv(OS.time, OS) ~ NUP93.expr, data = myinfo.sim.largecohort.1)
p1 <- ggsurvplot_facet(fit, data = myinfo.sim.largecohort.1, facet.by = "CancerType",
                       pval = TRUE, conf.int = FALSE, 
                       scales = "free_x", nrow = 1,
                       palette = "npg", title = "NUP93 expression (WT) p < 0.05")
pdf(file = paste("ggsurvplot.grouped.expression.noHM.nohyperCNB", "NUP93.sig", "pdf", sep = "."), 
    width = 15, height = 6)
print(p1)
dev.off()


# table(myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11"]$CNVgroup)
# table(myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11" & CNAburden < 20]$CNVgroup)
# table(myinfo.sim[STcode != "11" & CNAburden < 20, c(3,15)])
myinfo.sim.brca <- myinfo.sim[CancerType == "BRCA" & STcode != "11"]
# myinfo.sim.brca <- myinfo.sim[`cancer type abbreviation` == "BRCA" & STcode != "11"]
# myinfo.sim.brca.1 <- myinfo.sim.brca[CNVgroup %in% c("Amp/Amp", "Del/Del", "WT/WT")]
# myinfo.sim.brca.1 <- myinfo.sim.brca[STcode == "01" & CNVgroup %in% c("Del/Del", "WT/WT")]
# fit <- survfit(Surv(OS.time, OS) ~ CNVgroup + NUP93.expr, data = myinfo.sim.brca.1)
fit <- survfit(Surv(OS.time, OS) ~ NUP93.expr, data = myinfo.sim.brca)
# fit <- survfit(Surv(PFI.time, PFI) ~ CNVgroup, data = myinfo.sim.brca.1)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 # linetype = "strata", 
                 palette = "npg",
                 # palette = get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), 28),
                 title = "BRCA NUP93")

myinfo.sim.tmp <- myinfo.sim[CancerType == "COAD" & STcode != "11" & CNAburden < 20]
fit <- survfit(Surv(OS.time) ~ NUP93.expr + NUP93.CNV, data = myinfo.sim.tmp)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = FALSE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 palette = "npg",
                 title = "NUP93")


# NUPs mutation stats --------------------------------------------------------------
Nucleoporin.patients.count <- length(unique(mymuts[Hugo_Symbol %in% mygenes[Class=="Nucleoporin"]$targets]$Tumor_Sample_Barcode))
total.patients.count <- length(unique(mymuts$Tumor_Sample_Barcode))
# 2065/9104 0.2268234
Nucleoporin.patients.count/total.patients.count
NUP93.patients.count <- length(unique(mymuts[Hugo_Symbol == "NUP93"]$Tumor_Sample_Barcode))


NUP93.patients.count2 <- as.data.table(table(unique(mymuts[Hugo_Symbol == "NUP93", c(9,11)])$Cancer))
Nucleoporin.patients.count2 <- as.data.table(table(unique(mymuts[Hugo_Symbol %in% mygenes[Class=="Nucleoporin"]$targets, c(9,11)])$Cancer))
total.patients.count2 <- as.data.table(table(unique(mymuts[, c(9,11)])$Cancer))
# unique(mymuts[, c(9,11)])[is.na(Cancer)]
# 24
total.patients.count2[, Nucleoporin.mut.patients := Nucleoporin.patients.count2[match(total.patients.count2$V1, V1)]$N]
total.patients.count2[, Nucleoporin.mut.percentage := Nucleoporin.mut.patients/N*100]
total.patients.count2[, Mutated := paste(Nucleoporin.mut.patients, N, sep = "/")]
total.patients.count2[, NUP93.mut.patients := NUP93.patients.count2[match(total.patients.count2$V1, V1)]$N]
total.patients.count2[is.na(NUP93.mut.patients), NUP93.mut.patients := 0]
total.patients.count2[, NUP93.mut.percentage := NUP93.mut.patients/N*100]
total.patients.count2[, NUP93Mutated := paste(NUP93.mut.patients, N, sep = "/")]
total.patients.count2[, hyperMutated.patients := mutationload.count[hypermutation=="Y"][match(total.patients.count2$V1, Cancer)]$N]
total.patients.count2[, hyperMutated.percentage := hyperMutated.patients/N*100]
total.patients.count2[, hyperMutated := paste(hyperMutated.patients, N, sep = "/")]

hypermNUP93.patients.count <- as.data.table(table(unique(mymuts[Hugo_Symbol == "NUP93", c(9,11,13)])[,c(2:3)]))
hypermNUP93.patients.count[, sums := sum(N), by = "Cancer"]
hypermNUP93.patients.count[, percentage := N/sums * 100]
hypermNUP93.patients.count[, hyperMutated := paste(N, sums, sep = "/")]

p1 <- ggbarplot(data = total.patients.count2, x = "V1", y = "Nucleoporin.mut.percentage",
                label = "Mutated", lab.vjust = 0.5, lab.hjust = 0.5, # ylim = c(0, 60), 
                sort.val = "asc", color = "skyblue", fill = "skyblue", 
                xlab = "", ylab = "Nucleoporin mutated samples %",
                title = "Nucleoporin mutated in 2065/9104 of PAN-CANCER samples")
p11 <- ggpar(p1, rotate = T)
ggsave(filename = "ggbarplot.Nucleoporin.mutation.count.pdf", plot = p11, width = 10, height = 7)

p2 <- ggbarplot(data = total.patients.count2, x = "V1", y = "NUP93.mut.percentage",
                label = "NUP93Mutated", lab.vjust = 0.5, lab.hjust = 0.5, # ylim = c(0, 60), 
                sort.val = "asc", color = "coral", fill = "coral", 
                xlab = "", ylab = "NUP93 mutated samples %",
                title = "NUP93 mutated in 134/9104 of PAN-CANCER samples")
p22 <- ggpar(p2, rotate = T)
ggsave(filename = "ggbarplot.Nucleoporin.NUP93.mutation.count.pdf", 
       plot = ggarrange(p11, p22, nrow = 1), width = 15, height = 7)

p3 <- ggbarplot(data = hypermNUP93.patients.count[hypermutation == "HM"], x = "Cancer", y = "percentage",
                label = "hyperMutated", lab.vjust = 0.5, lab.hjust = 0.5, 
                sort.val = "asc", color = "brown", fill = "brown", 
                xlab = "", ylab = "hypermutation % out of NUP93 mutated samples",
                title = "67 hypermutated in 134/9104 of PAN-CANCER samples")
p33 <- ggpar(p3, rotate = T)
ggsave(filename = "ggbarplot.Nucleoporin.NUP93.hypermutation.count.pdf", 
       plot = ggarrange(p11, p22, p33, nrow = 1), width = 21, height = 7)

p4 <- ggbarplot(data = total.patients.count2, x = "V1", y = "hyperMutated.percentage",
                label = "hyperMutated", lab.vjust = 0.5, lab.hjust = 0.5, # ylim = c(0, 60), 
                sort.val = "asc", color = "seagreen", fill = "seagreen", 
                xlab = "", ylab = "hypermutated samples %",
                title = "hypermutated in 456/9104 of PAN-CANCER samples")
p44 <- ggpar(p4, rotate = T)
ggsave(filename = "ggbarplot.hypermutation.count.pdf", plot = p44, width = 10, height = 7)

# unique(mymuts[Hugo_Symbol %in% mygenes[Class=="Nucleoporin"]$targets, c(9,11,13)][hypermutation == "HM"])
# 429
# unique(mymuts[, c(9,11,13)][hypermutation == "HM"])
# 456
# unique(mymuts[Hugo_Symbol %in% mygenes[Class=="Nucleoporin"]$targets, c(9,11,13)][hypermutation == "NotHM"])
# 1636
# unique(mymuts[, c(9,11,13)][hypermutation == "NotHM"])
# 8648
# unique(mymuts[, c(9,11)])[is.na(Cancer)]
# 24
# > sum(total.patients.count3$N)
# [1] 8624
# > total.patients.count - 456
# [1] 8648
# > total.patients.count - 456 - 24
# [1] 8624
Nucleoporin.patients.count3 <- as.data.table(table(unique(mymuts[Hugo_Symbol %in% mygenes[Class=="Nucleoporin"]$targets, c(9,11,13)])[hypermutation == "NotHM"]$Cancer))
NUP93.patients.count3 <- as.data.table(table(unique(mymuts[Hugo_Symbol == "NUP93", c(9,11,13)])[hypermutation == "NotHM"]$Cancer))
total.patients.count3 <- as.data.table(table(unique(mymuts[hypermutation == "NotHM", c(9,11)])$Cancer))
total.patients.count3[, Nucleoporin.mut.patients := Nucleoporin.patients.count3[match(total.patients.count3$V1, V1)]$N]
total.patients.count3[, Nucleoporin.mut.percentage := Nucleoporin.mut.patients/N*100]
total.patients.count3[, Mutated := paste(Nucleoporin.mut.patients, N, sep = "/")]
total.patients.count3[, NUP93.mut.patients := NUP93.patients.count3[match(total.patients.count3$V1, V1)]$N]
total.patients.count3[is.na(NUP93.mut.patients), NUP93.mut.patients := 0]
total.patients.count3[, NUP93.mut.percentage := NUP93.mut.patients/N*100]
total.patients.count3[, NUP93Mutated := paste(NUP93.mut.patients, N, sep = "/")]
p1 <- ggbarplot(data = total.patients.count3, x = "V1", y = "Nucleoporin.mut.percentage",
                label = "Mutated", lab.vjust = 0.5, lab.hjust = 0.5, # ylim = c(0, 60), 
                sort.val = "asc", color = "skyblue", fill = "skyblue", 
                xlab = "", ylab = "Nucleoporin mutated samples %",
                title = "Nucleoporin mutated in 1635/8624 of non-hypermutated PAN-CANCER samples")
p11 <- ggpar(p1, rotate = T)
ggsave(filename = "ggbarplot.NotHM.Nucleoporin.mutation.count.pdf", plot = p11, width = 10, height = 7)

p2 <- ggbarplot(data = total.patients.count3, x = "V1", y = "NUP93.mut.percentage",
                label = "NUP93Mutated", lab.vjust = 0.5, lab.hjust = 0.5, # ylim = c(0, 60), 
                sort.val = "asc", color = "coral", fill = "coral", 
                xlab = "", ylab = "NUP93 mutated samples %",
                title = "NUP93 mutated in 67/8624 of non-hypermutated PAN-CANCER samples")
p22 <- ggpar(p2, rotate = T)
ggsave(filename = "ggbarplot.NotHM.Nucleoporin.NUP93.mutation.count.pdf", 
       plot = ggarrange(p11, p22, nrow = 1), width = 15, height = 7)


