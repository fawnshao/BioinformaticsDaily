setwd("~/Data/myworkData/pancancer.data/")
# Tumor_Sample_Barcode column not found in provided clinical data. Rename column name containing sample names to Tumor_Sample_Barcode if necessary.
library(maftools)
library(data.table)
library(ggpubr)
# maftitles <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", 
#                "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "Gene_Type", 
#                "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Cancer_Type")
mygenes <- fread(input = "NUPs.gene.txt")
# mymuts <- fread(input = "mc3.v0.2.8.PUBLIC.xena")
# myinfo <- fread(input = "Survival_SupplementalTable_S1_20171025_xena_sp")

# mygenes.muts <- mymuts[gene %in% mygenes$targets]
# Samples   10295 
# notice: match Tumor_Sample_Barcode
# tcga.muts <- read.maf(maf = "mc3.v0.2.8.PUBLIC.maf.gz", clinicalData = "formatted.clinical_PANCAN_patient_with_followup.tsv")
# -Finished in 00:05:27 elapsed (00:06:14 cpu) 
tcga.muts <- read.maf(maf = "formatted.mc3.v0.2.8.PUBLIC.xena.maf", 
                      clinicalData = "Survival_SupplementalTable_S1_20171025_xena_sp")
getSampleSummary(tcga.muts)
getGeneSummary(tcga.muts)
getClinicalData(tcga.muts)
getFields(tcga.muts)
pdf(file ="NUPgenes.mc3.v0.2.8.PUBLIC.maf.oncoplot.pdf", width = 16, height = 12)
oncoplot(maf = tcga.muts, genes = mygenes$targets)
dev.off()

png(filename ="NUPgenes.mc3.v0.2.8.PUBLIC.maf.oncoplot.png", width = 1200, height = 800)
oncoplot(maf = tcga.muts, genes = mygenes$targets)
dev.off()


lollipopPlot(maf = tcga.muts, gene = 'TPR', AACol = 'Amino_Acid_Change', showMutationRate = TRUE)
lollipopPlot(maf = tcga.muts, gene = 'NUP93', AACol = 'Amino_Acid_Change', showMutationRate = TRUE, proteinID = "NP_055484")

#Survival analysis based on grouping of mutation status
mafSurvival(maf = tcga.muts, genes = 'NUP93', time = 'OS.time', Status = 'OS', isTCGA = FALSE)




# stats --------------------------------------------------------------------
mutsamples <- fread(input = "NUP93.mut.samples", header = F)
mutsamples.df <- dcast(data = mutsamples, formula = V1 ~ V3, value.var = "V2")
mutsamples.df[, percent := NUP93mut / All * 100]
mutsamples.df[, labs := paste(NUP93mut, All, sep = "/")]
p1 <- ggbarplot(mutsamples.df[!is.na(NUP93mut)], x = "V1", y = "percent", 
                label = "labs", color = "steelblue", fill = "steelblue",
                xlab = "", ylab = "NUP93 mutated samples%")
ggsave(filename = "ggbarplot.NUP93.mut.samples.pdf", width = 12, height = 6, p1)

# UCEC --------------------------------------------------------------------
ucec.muts <- read.maf(maf = "UCEC.formatted.mc3.v0.2.8.PUBLIC.xena.maf", 
                      clinicalData = "UCEC.Survival_SupplementalTable_S1_20171025_xena_sp")
pdf(file ="UCEC.NUPgenes.mc3.v0.2.8.PUBLIC.maf.mafSurvival.pdf", width = 10, height = 5)
mafSurvival(maf = ucec.muts, genes = 'NUP93', time = 'OS.time', Status = 'OS', isTCGA = FALSE)
dev.off()

pdf(file ="UCEC.TP53.mc3.v0.2.8.PUBLIC.maf.mafSurvival.pdf", width = 10, height = 5)
mafSurvival(maf = ucec.muts, genes = 'TP53', time = 'OS.time', Status = 'OS', isTCGA = FALSE)
dev.off()

mafSurvGroup(maf = ucec.muts, geneSet = c("NUP93", "CTCF","TPR", "RAD21"), time = "OS.time", Status = "OS")

# mafSurvival(maf = ucec.muts, genes = 'NUP93', time = 'PFI.time', Status = 'PFI', isTCGA = FALSE)
pdf(file ="UCEC.NUPgenes.mc3.v0.2.8.PUBLIC.maf.oncoplot.pdf", width = 16, height = 12)
oncoplot(maf = ucec.muts, genes = mygenes$targets)
dev.off()



# prog_geneset <- survGroup(maf = tcga.muts, top = 20, geneSetSize = 2, time = "OS.time", Status = "OS", verbose = FALSE)
# print(prog_geneset)
# mafSurvGroup(maf = tcga.muts, geneSet = mygenes$targets, time = "OS.time", Status = "OS")












# plotmafSummary(maf = laml.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# lollipopPlot(maf = tcga.muts, gene = 'NUP93', AACol = 'HGVSp_Short', showMutationRate = TRUE, proteinID = "NP_055484")
# lollipopPlot(maf = laml.maf, gene = 'NUP93', showMutationRate = TRUE, refSeqID = "NM_014669")
# lollipopPlot(maf = tcga.muts, gene = 'TPR', AACol = 'HGVSp_Short', showMutationRate = TRUE)


# ## ------------------------------------------------------------------------
# #path to TCGA LAML MAF file
# laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
# #clinical information containing survival information and histology. This is optional
# laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
# 
# laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
# 
# ## ------------------------------------------------------------------------
# #Typing laml shows basic summary of MAF file.
# laml
# #Shows sample summry.
# getSampleSummary(laml)
# #Shows gene summary.
# getGeneSummary(laml)
# #shows clinical data associated with samples
# getClinicalData(laml)
# #Shows all fields in MAF
# getFields(laml)
# #Writes maf summary to an output file with basename laml.
# write.mafSummary(maf = laml, basename = 'laml')
# 
# ## ----fig.height=5, fig.width=6-------------------------------------------
# plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# 
# ## ---- fig.align='left',fig.height=5,fig.width=10, fig.align='left'-------
# #oncoplot for top ten mutated genes.
# oncoplot(maf = laml, top = 10)
# 
# ## ---- fig.height=2.4,fig.width=8,fig.align='left'------------------------
# oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'))
# 
# ## ---- fig.height=5, fig.width=6, eval = T, fig.align='left'--------------
# laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
# #plot titv summary
# plotTiTv(res = laml.titv)
# 
# ## ----fig.align='left', fig.width=6, fig.height=3-------------------------
# #lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
# lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)
# 
# ## ----fig.align='left', fig.width=6, fig.height=3-------------------------
# lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = 816, refSeqID = 'NM_000222')
# 
# ## ---- results='hide', message=FALSE--------------------------------------
# brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
# brca = read.maf(maf = brca, verbose = FALSE)
# 
# ## ---- fig.height=5,fig.width=12,fig.align='center'-----------------------
# rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.6)
# 
# ## ---- fig.align='left', fig.height=5, fig.width=12, message=FALSE, results='hide'----
# laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML')
# 
# ## ---- fig.align='left', fig.height=4, fig.width=4------------------------
# plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')
# 
# ## ---- fig.align='left',fig.width=7, fig.height=5, eval=T-----------------
# geneCloud(input = laml, minMut = 3)
# 
# ## ------------------------------------------------------------------------
# all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
# amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
# del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
# scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")
# 
# laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
# 
# #GISTIC object
# laml.gistic
# 
# ## ---- fig.width=6, fig.height=4, fig.align='left'------------------------
# gisticChromPlot(gistic = laml.gistic, markBands = "all")
# 
# ## ---- fig.width=5, fig.height=4, fig.align='left'------------------------
# gisticBubblePlot(gistic = laml.gistic)
# 
# ## ---- fig.align='left',fig.width=7, fig.height=5, eval=T-----------------
# gisticOncoPlot(gistic = laml.gistic, clinicalData = getClinicalData(x = laml), clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10)
# 
# ## ---- fig.height=3,fig.width=8,fig.align='center'------------------------
# tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.3009.hg19.seg.txt", package = "maftools")
# plotCBSsegments(cbsFile = tcga.ab.009.seg)
# 
# ## ---- message=FALSE------------------------------------------------------
# #exclusive/co-occurance event analysis on top 10 mutated genes. 
# somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
# 
# ## ---- fig.align='default', fig.width=7,fig.height=5, message=F,results='hide', eval=T----
# laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
# head(laml.sig)
# 
# ## ---- fig.align='left', fig.width=5, fig.height=4------------------------
# plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
# 
# ## ---- fig.align='left', fig.width=5, fig.height=4------------------------
# laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)
# #Protein summary (Printing first 7 columns for display convenience)
# laml.pfam$proteinSummary[,1:7, with = FALSE]
# #Domain summary (Printing first 3 columns for display convenience)
# laml.pfam$domainSummary[,1:3, with = FALSE]
# 
# ## ---- fig.width=5, fig.height=4------------------------------------------
# #MutsigCV results for TCGA-AML
# laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
# pancanComparison(mutsigResults = laml.mutsig, qval = 0.1, cohortName = 'LAML', inputSampleSize = 200, label = 1)
# 
# ## ---- fig.width=5, fig.height=5------------------------------------------
# #Survival analysis based on grouping of DNMT3A mutation status
# mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)
# 
# ## ------------------------------------------------------------------------
# #Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
# prog_geneset = survGroup(maf = laml, top = 20, geneSetSize = 2, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)
# 
# print(prog_geneset)
# 
# ## ---- fig.width=5, fig.height=5------------------------------------------
# mafSurvGroup(maf = laml, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", Status = "Overall_Survival_Status")
# 
# ## ----results='hide', message=FALSE---------------------------------------
# #Primary APL MAF
# primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
# primary.apl = read.maf(maf = primary.apl)
# #Relapse APL MAF
# relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
# relapse.apl = read.maf(maf = relapse.apl)
# 
# ## ---- fig.align='left'---------------------------------------------------
# #Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
# pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
# print(pt.vs.rt)
# 
# ## ---- fig.width=5, fig.height=5, fig.align='left'------------------------
# forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
# 
# ## ---- fig.height=3,fig.width=11, eval=T, fig.align='left'----------------
# genes = c("PML", "RARA", "RUNX1", "ARID1B", "FLT3")
# coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', genes = genes, removeNonMutated = TRUE)
# 
# ## ---- warning=FALSE, message=FALSE,fig.align='left', results='hide'------
# lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "PML", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")
# 
# ## ------------------------------------------------------------------------
# fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'FAB_classification')
# 
# #Results are returned as a list. Significant associations p-value < 0.05
# fab.ce$groupwise_comparision[p_value < 0.05]
# 
# ## ---- fig.width=6, fig.height=4------------------------------------------
# plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)
# 
# ## ---- fig.height=4, fig.width=8------------------------------------------
# dgi = drugInteractions(maf = laml, fontSize = 0.75)
# 
# ## ------------------------------------------------------------------------
# dnmt3a.dgi = drugInteractions(genes = "DNMT3A", drugs = TRUE)
# #Printing selected columns.
# dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
# 
# ## ---- fig.width=4, fig.height=4------------------------------------------
# OncogenicPathways(maf = laml)
# 
# ## ---- fig.width=10, fig.height=3-----------------------------------------
# PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")
# 
# ## ---- echo = TRUE, fig.align='left', fig.height=4, fig.width=6, eval=T----
# #Heterogeneity in sample TCGA.AB.2972
# tcga.ab.2972.het = inferHeterogeneity(maf = laml, tsb = 'TCGA-AB-2972', vafCol = 'i_TumorVAF_WU')
# print(tcga.ab.2972.het$clusterMeans)
# #Visualizing results
# plotClusters(clusters = tcga.ab.2972.het)
# 
# ## ---- fig.align='left', fig.height=4, fig.width=6, eval=T----------------
# seg = system.file('extdata', 'TCGA.AB.3009.hg19.seg.txt', package = 'maftools')
# tcga.ab.3009.het = inferHeterogeneity(maf = laml, tsb = 'TCGA-AB-3009', segFile = seg, vafCol = 'i_TumorVAF_WU')
# #Visualizing results. Highlighting those variants on copynumber altered variants.
# plotClusters(clusters = tcga.ab.3009.het, genes = 'CN_altered', showCNvars = TRUE)
# 
# ## ---- eval=TRUE----------------------------------------------------------
# #Requires BSgenome object
# library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
# laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
# 
# ## ---- eval=TRUE, fig.height=4, fig.width=7-------------------------------
# plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)
# 
# ## ---- echo=FALSE---------------------------------------------------------
# par(mar = c(2, 2, 2, 1))
# plot(NA, xlim = c(1, 10), ylim = c(0, 30), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
# rect(xleft = 3, ybottom = 28, xright = 7, ytop = 30, col = grDevices::adjustcolor("gray70", alpha.f = 0.6), lwd = 1.2, border = "maroon")
# text(x = 5, y = 29, labels = "MAF", font = 2)
# arrows(x0 = 5, y0 = 28, x1 = 5, y1 = 26, length = 0.1, lwd = 2)
# text(x = 5, y = 25, labels = "trinucleotideMatrix()", font = 3)
# arrows(x0 = 5, y0 = 24, x1 = 5, y1 = 21, length = 0.1, lwd = 2)
# text(x = 5, y = 20, labels = "estimateSignatures()", font = 3)
# arrows(x0 = 5, y0 = 19, x1 = 5, y1 = 16, length = 0.1, lwd = 2)
# text(x = 5, y = 15, labels = "plotCophenetic()", font = 3)
# arrows(x0 = 5, y0 = 14, x1 = 5, y1 = 11, length = 0.1, lwd = 2)
# text(x = 5, y = 10, labels = "extractSignatures()", font = 3)
# arrows(x0 = 5, y0 = 9, x1 = 5, y1 = 6, length = 0.1, lwd = 2)
# text(x = 5, y = 5, labels = "compareSignatures()", font = 3)
# arrows(x0 = 5, y0 = 4, x1 = 5, y1 = 1, length = 0.1, lwd = 2)
# text(x = 5, y = 0, labels = "plotSignatures()", font = 3)
# 
# ## ---- fig.height=5, fig.width=5, eval=FALSE, message=FALSE---------------
# #  library('NMF')
# #  laml.sign = estimateSignatures(mat = laml.tnm, nTry = 4)
# 
# ## ---- fig.height=3, fig.width=3, eval=TRUE, message=FALSE, echo=FALSE, include=FALSE----
# #Run main function with maximum 6 signatures. 
# library('NMF')
# laml.sign = estimateSignatures(mat = laml.tnm, nTry = 4, pConstant = 0.1, plotBestFitRes = FALSE)
# 
# ## ---- fig.width=3, fig.height=3------------------------------------------
# plotCophenetic(res = laml.sign)
# 
# ## ---- echo=FALSE---------------------------------------------------------
# laml.sig = extractSignatures(mat = laml.tnm, n = 3, pConstant = 0.1)
# 
# ## ------------------------------------------------------------------------
# #Compate against original 30 signatures 
# laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
# #Compate against updated version3 60 signatures 
# laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")
# 
# ## ---- fig.width=7, fig.height=2.5, fig.align='center'--------------------
# library('pheatmap')
# pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
# 
# ## ---- fig.width=6, fig.height=4, fig.align='center', eval = T------------
# maftools::plotSignatures(nmfRes = laml.sig, title_size = 0.8)
# 
# ## ---- echo=FALSE,eval=FALSE----------------------------------------------
# #  colnames(laml.sign$contributions) = as.character(getSampleSummary(x = laml)[,Tumor_Sample_Barcode])[1:187]
# 
# ## ---- fig.height=3.5, fig.width=5, warning=FALSE-------------------------
# laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig)
# 
# ## ---- fig.height=4, fig.width=6------------------------------------------
# plotEnrichmentResults(enrich_res = laml.se, pVal = 0.05)
# 
# ## ---- eval=T-------------------------------------------------------------
# var.annovar = system.file("extdata", "variants.hg19_multianno.txt", package = "maftools")
# var.annovar.maf = annovarToMaf(annovar = var.annovar, Center = 'CSI-NUS', refBuild = 'hg19', 
#                                tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene')
# 
# 
# ## ------------------------------------------------------------------------
# #Read sample ICGC data for ESCA
# esca.icgc <- system.file("extdata", "simple_somatic_mutation.open.ESCA-CN.sample.tsv.gz", package = "maftools")
# esca.maf <- icgcSimpleMutationToMAF(icgc = esca.icgc, addHugoSymbol = TRUE)
# #Printing first 16 columns for display convenience.
# print(esca.maf[1:5,1:16, with = FALSE])
# 
# ## ---- eval=FALSE---------------------------------------------------------
# #  laml.mutsig.corrected = prepareMutSig(maf = laml)
# #  # Converting gene names for 1 variants from 1 genes
# #  #    Hugo_Symbol MutSig_Synonym N
# #  # 1:    ARHGAP35          GRLF1 1
# #  # Original symbols are preserved under column OG_Hugo_Symbol.
# 
# ## ------------------------------------------------------------------------
# #Extract data for samples 'TCGA.AB.3009' and 'TCGA.AB.2933'  (Printing just 5 rows for display convenience)
# subsetMaf(maf = laml, tsb = c('TCGA-AB-3009', 'TCGA-AB-2933'), mafObj = FALSE)[1:5]
# ##Same as above but return output as an MAF object (Default behaviour)
# subsetMaf(maf = laml, tsb = c('TCGA-AB-3009', 'TCGA-AB-2933'))
# 
# ## ------------------------------------------------------------------------
# #Select all Splice_Site mutations from DNMT3A and NPM1
# subsetMaf(maf = laml, genes = c('DNMT3A', 'NPM1'), mafObj = FALSE,query = "Variant_Classification == 'Splice_Site'")
# 
# #Same as above but include only 'i_transcript_name' column in the output
# subsetMaf(maf = laml, genes = c('DNMT3A', 'NPM1'), mafObj = FALSE, query = "Variant_Classification == 'Splice_Site'", fields = 'i_transcript_name')
# 
# ## ---- eval=FALSE---------------------------------------------------------
# #  devtools::install_github(repo = "PoisonAlien/TCGAmutations")
# 
# ## ------------------------------------------------------------------------
# sessionInfo()
# 



#path to TCGA LAML MAF file
# laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
# laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

