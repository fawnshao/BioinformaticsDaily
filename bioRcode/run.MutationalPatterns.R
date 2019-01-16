# whole mut, highTIP Mut, lowTIP Mut, whle eMut, highTIP eMut, lowTIP eMut
myvcf <- commandArgs(TRUE)
# myvcf <- "hg19.BRCA-US.enhancer.mut.vcf"
# myvcf <- c("hg19.BRCA-US.all.mut.vcf", "hg19.BRCA-US.all.highTIP60.mut.vcf", "hg19.BRCA-US.all.lowTIP60.mut.vcf", 
#            "hg19.BRCA-US.enhancer.mut.vcf", "hg19.BRCA-US.enhancer.highTIP60.mut.vcf", 
#            "hg19.BRCA-US.enhancer.lowTIP60.mut.vcf")

###################################################
### code chunk number 1: convert bed file to vcf
###################################################
# library(bedr)

###################################################
### code chunk number 3: loading_reference_data
###################################################
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", version = "3.8")
library(BSgenome)
head(available.genomes())


###################################################
### code chunk number 4: loading_reference_data
###################################################
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)


###################################################
### code chunk number 5: load_package
###################################################
library(MutationalPatterns)


###################################################
### code chunk number 6: locate_files
###################################################
# vcf_files <- list.files(system.file("extdata", package="MutationalPatterns"),
#                         pattern = ".vcf", full.names = TRUE)


###################################################
### code chunk number 7: set_sample_names
###################################################
sample_names <- c("all.MUT", "all.highTIP60.MUT", "all.lowTIP60.MUT",
                  "all.eMUT", "all.highTIP60.eMUT", "all.lowTIP60.eMUT")

###################################################
### code chunk number 8: read_vcfs_as_granges
###################################################
vcfs <- read_vcfs_as_granges(myvcf, sample_names, ref_genome)
summary(vcfs)


###################################################
### code chunk number 9: store_tissue_variable
###################################################
tissue <- sample_names
tissue <- factor(tissue, levels = unique(tissue))


###################################################
### code chunk number 10: mutations_from_vcf
###################################################
muts = mutations_from_vcf(vcfs[[1]])
# head(muts, 12)


###################################################
### code chunk number 11: mut_type
###################################################
types = mut_type(vcfs[[1]])
# head(types, 12)


###################################################
### code chunk number 12: mut_context
###################################################
context = mut_context(vcfs[[1]], ref_genome)
# head(context, 12)


###################################################
### code chunk number 13: type_context
###################################################
type_context = type_context(vcfs[[1]], ref_genome)
# lapply(type_context, head, 12)


###################################################
### code chunk number 14: mut_type_occurrences
###################################################
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
# type_occurrences


###################################################
### code chunk number 15: plot_spectrum
###################################################
p1 <- plot_spectrum(type_occurrences)
# pdf(file = paste("MutationalPatterns.plot_spectrum", myvcf[1], "pdf", sep = "."), width = 6, height = 6)
# print(p1)
# dev.off()


###################################################
### code chunk number 16: plot_spectrum_2
###################################################
p2 <- plot_spectrum(type_occurrences, CT = TRUE)
# pdf(file = paste("MutationalPatterns.plot_spectrum.CT", myvcf[1], "pdf", sep = "."), width = 6, height = 6)
# print(p2)
# dev.off()

###################################################
### code chunk number 17: plot_spectrum_3
###################################################
p3 <- plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)


###################################################
### code chunk number 18: combine_plot_spectrum_noeval (eval = FALSE)
###################################################
## library("gridExtra")
## grid.arrange(p1, p2, p3, ncol=3, widths=c(3,3,1.75))


###################################################
### code chunk number 19: combine_plot_spectrum
###################################################
library(gridExtra)
library(ggplot2)
ggsave(paste("MutationalPatterns.plot_spectrum", myvcf[1], "pdf", sep = "."),
       grid.arrange(p1, p2, p3, ncol = 3, widths = c(3,3,1.75)),
       width = 10,
       height = 3)


###################################################
### code chunk number 20: plot_spectrum_4
###################################################
p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)
ggsave(paste("MutationalPatterns.plot_spectrum.byTissue", myvcf[1], "pdf", sep = "."),
       p4, width = 10, height = 6)

###################################################
### code chunk number 21: plot_spectrum_5
###################################################
# palette <- c("pink", "orange", "blue", "lightblue", "green", "red", "purple")
# p5 <- plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE, colors = palette)


###################################################
### code chunk number 22: combine_plot_spectrum_2_noeval
###################################################
# ggsave(paste("MutationalPatterns.plot_spectrum.byTissue", myvcf[1], "pdf", sep = "."),
#        grid.arrange(p4, p5, ncol = 2, widths = c(4,2.3)), 
#        width = 10, 
#        height = 3)


###################################################
### code chunk number 23: combine_plot_spectrum_2 (eval = FALSE)
###################################################
## grid.arrange(p4, p5, ncol=2, widths=c(4,2.3))


###################################################
### code chunk number 24: mut_matrix
###################################################
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
# head(mut_mat)


###################################################
### code chunk number 25: plot_96_profile_2_e (eval = FALSE)
###################################################
## plot_96_profile(mut_mat[,c(1,7)])


###################################################
### code chunk number 26: plot_96_profile_2
###################################################
# plot_96_profile(mut_mat[,c(1,7)])


###################################################
### code chunk number 27: plot_96_profile_3_e (eval = FALSE)
###################################################
## plot_96_profile(mut_mat[,c(1,7)], condensed = TRUE)


###################################################
### code chunk number 28: plot_96_profile_3
###################################################
pdf(file = paste("MutationalPatterns.plot_96_profile", myvcf[1], "pdf", sep = "."), width = 12, height = 8)
plot_96_profile(mut_mat, condensed = TRUE)
dev.off()

###################################################
### code chunk number 29: psuedo_count
###################################################
mut_mat <- mut_mat + 0.0001


###################################################
### code chunk number 30: use_nmf
###################################################
library("NMF")
estimate <- nmf(mut_mat, rank = 2:5, method = "brunet", nrun = 10, seed = 123456)


###################################################
### code chunk number 31: estimate_rank_e (eval = FALSE)
###################################################
## plot(estimate)


###################################################
### code chunk number 32: estimate_rank
###################################################
pdf(file = paste("MutationalPatterns.NMF.estimate", myvcf[1], "pdf", sep = "."), width = 10, height = 8)
plot(estimate)
dev.off()

###################################################
### code chunk number 33: extract_signatures
###################################################
nmf_res <- extract_signatures(mut_mat, rank = 5, nrun = 10)


###################################################
### code chunk number 34: add_column_names
###################################################
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")


###################################################
### code chunk number 35: plot_96_profile_e (eval = FALSE)
###################################################
## plot_96_profile(nmf_res$signatures, condensed = TRUE)


###################################################
### code chunk number 36: plot_96_profile
###################################################
pdf(file = paste("MutationalPatterns.NMF.plot_96_profile", myvcf[1], "pdf", sep = "."), width = 12, height = 8)
plot_96_profile(nmf_res$signatures, condensed = TRUE)
dev.off()
pdf(file = paste("MutationalPatterns.NMF.plot_96_profile.contriubution", myvcf[1], "pdf", sep = "."), width = 12, height = 8)
plot_96_profile(nmf_res$contribution, condensed = TRUE)
dev.off()

###################################################
### code chunk number 37: plot_contribution
###################################################
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature,
                         mode = "relative")


###################################################
### code chunk number 38: plot_contribution_2
###################################################
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature,
                         mode = "absolute")


###################################################
### code chunk number 39: plot_contribution_2_fig_e (eval = FALSE)
###################################################
## grid.arrange(pc1, pc2)


###################################################
### code chunk number 40: plot_contribution_2_fig
###################################################
# grid.arrange(pc1, pc2)
ggsave(paste("MutationalPatterns.plot_contribution", myvcf[1], "pdf", sep = "."),
       grid.arrange(pc1, pc2),
       width = 10,
       height = 8)


###################################################
### code chunk number 41: plot_contribution_4_e (eval = FALSE)
###################################################
## plot_contribution(nmf_res$contribution, nmf_res$signature,
##                          mode = "absolute", coord_flip = TRUE)


###################################################
### code chunk number 42: plot_contribution_4
###################################################
p <- plot_contribution(nmf_res$contribution, nmf_res$signature,
                       mode = "absolute", coord_flip = TRUE)
ggsave(paste("MutationalPatterns.plot_contribution.absolute", myvcf[1], "pdf", sep = "."),
       grid.arrange(p),
       width = 10,
       height = 8)

###################################################
### code chunk number 43: plot_contribution_heatmap
###################################################
pch1 <- plot_contribution_heatmap(nmf_res$contribution,
                                  sig_order = c("Signature E", "Signature D", "Signature C", "Signature B", "Signature A"))


###################################################
### code chunk number 44: plot_contribution_heatmap2_e
###################################################
pch2 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=FALSE)


###################################################
### code chunk number 45: plot_contribution_heatmap2_e2 (eval = FALSE)
###################################################
## grid.arrange(pch1, pch2, ncol = 2, widths = c(2,1.6))


###################################################
### code chunk number 46: plot_contribution_heatmap2
###################################################
ggsave(paste("MutationalPatterns.plot_contribution_heatmap", myvcf[1], "pdf", sep = "."),
       grid.arrange(pch1, pch2, ncol = 2, widths = c(2,1.6)),
       width = 10,
       height = 8)

###################################################
### code chunk number 47: plot_compare_profiles_e (eval = FALSE)
###################################################
## plot_compare_profiles(mut_mat[,1],
##                         nmf_res$reconstructed[,1],
##                         profile_names = c("Original", "Reconstructed"),
##                         condensed = TRUE)


###################################################
### code chunk number 48: plot_compare_profiles
###################################################
p <- plot_compare_profiles(mut_mat[,1], 
                           nmf_res$reconstructed[,1],
                           profile_names = c("Original", "Reconstructed"),
                           condensed = TRUE)
ggsave(paste("MutationalPatterns.plot_compare_profiles", myvcf[1], "pdf", sep = "."),
       grid.arrange(p),
       width = 10,
       height = 8)

# ###################################################
# ### code chunk number 49: download_cancer_signatures
# ###################################################
# sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
#                 "signatures_probabilities.txt", sep = "")
# 
# cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# # Match the order of the mutation types to MutationalPatterns standard
# new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# # Reorder cancer signatures dataframe
# cancer_signatures = cancer_signatures[as.vector(new_order),]
# # Add trinucletiode changes names as row.names
# row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# # Keep only 96 contributions of the signatures in matrix
# cancer_signatures = as.matrix(cancer_signatures[,4:33])
# 
# 
# ###################################################
# ### code chunk number 50: plot_96_profile_COSMIC_e (eval = FALSE)
# ###################################################
# ## plot_96_profile(cancer_signatures[,1:2], condensed = TRUE, ymax = 0.3)
# 
# 
# ###################################################
# ### code chunk number 51: plot_96_profile_COSMIC
# ###################################################
# plot_96_profile(cancer_signatures[,1:2], condensed = TRUE, ymax = 0.3)
# 
# 
# ###################################################
# ### code chunk number 52: cluster_COSMIC
# ###################################################
# hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
# # store signatures in new order
# cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
# plot(hclust_cosmic)
# 
# 
# ###################################################
# ### code chunk number 53: cos_sim
# ###################################################
# cos_sim(mut_mat[,1], cancer_signatures[,1])
# 
# 
# ###################################################
# ### code chunk number 54: cos_sim_cosmic_samples
# ###################################################
# cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
# # Plot heatmap with specified signature order
# plot_cosine_heatmap(cos_sim_samples_signatures,
#                     col_order = cosmic_order,
#                     cluster_rows = TRUE)
# 
# 
# ###################################################
# ### code chunk number 55: fit_to_signatures
# ###################################################
# fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
# 
# 
# ###################################################
# ### code chunk number 56: plot_contribution_3_noeval (eval = FALSE)
# ###################################################
# ## # Select signatures with some contribution
# ## select <- which(rowSums(fit_res$contribution) > 10)
# ## # Plot contribution barplot
# ## plot_contribution(fit_res$contribution[select,],
# ##                     cancer_signatures[,select],
# ##                     coord_flip = FALSE,
# ##                     mode = "absolute")
# 
# 
# ###################################################
# ### code chunk number 57: plot_contribution_3
# ###################################################
# # Select signatures with some contribution
# select <- which(rowSums(fit_res$contribution) > 0)
# # Plot contribution
# ggsave("plot_contribution_3.pdf",
#        plot_contribution(fit_res$contribution[select,],
#                          cancer_signatures[,select],
#                          coord_flip = FALSE,
#                          mode = "absolute"),
#        width=9,
#        height=5)
# 
# 
# ###################################################
# ### code chunk number 58: plot_contribution_heatmap3_e (eval = FALSE)
# ###################################################
# ## plot_contribution_heatmap(fit_res$contribution,
# ##                           cluster_samples = TRUE,
# ##                           method = "complete")
# 
# 
# ###################################################
# ### code chunk number 59: plot_contribution_heatmap3
# ###################################################
# plot_contribution_heatmap(fit_res$contribution,
#                           cluster_samples = TRUE,
#                           method = "complete")
# 
# 
# ###################################################
# ### code chunk number 60: plot_compare_profiles_2_e (eval = FALSE)
# ###################################################
# ## plot_compare_profiles(mut_mat[,1], fit_res$reconstructed[,1],
# ##                         profile_names = c("Original", "Reconstructed"),
# ##                         condensed = TRUE)
# 
# 
# ###################################################
# ### code chunk number 61: plot_compare_profiles_2
# ###################################################
# plot_compare_profiles(mut_mat[,1], fit_res$reconstructed[,1],
#                       profile_names = c("Original", "Reconstructed"),
#                       condensed = TRUE)
# 
# 
# ###################################################
# ### code chunk number 62: cos_sim_ori_rec
# ###################################################
# # calculate all pairwise cosine similarities
# cos_sim_ori_rec <- cos_sim_matrix(mut_mat, fit_res$reconstructed)
# # extract cosine similarities per sample between original and reconstructed
# cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
# 
# 
# ###################################################
# ### code chunk number 63: cos_sim_ori_rec
# ###################################################
# # Adjust data frame for plotting with gpplot
# colnames(cos_sim_ori_rec) = "cos_sim"
# cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
# 
# 
# ###################################################
# ### code chunk number 64: plot_cos_sim_ori_rec_e (eval = FALSE)
# ###################################################
# ## # Load ggplot2
# ## library(ggplot2)
# ## # Make barplot
# ## ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
# ##   geom_bar(stat="identity", fill = "skyblue4") +
# ##   coord_cartesian(ylim=c(0.8, 1)) +
# ##   # coord_flip(ylim=c(0.8,1)) +
# ##   ylab("Cosine similarity\n original VS reconstructed") +
# ##   xlab("") +
# ##   # Reverse order of the samples such that first is up
# ##   # xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
# ##   theme_bw() +
# ##   theme(panel.grid.minor.y=element_blank(),
# ##         panel.grid.major.y=element_blank()) +
# ##   # Add cut.off line
# ##   geom_hline(aes(yintercept=.95))
# 
# 
# ###################################################
# ### code chunk number 65: plot_cos_sim_ori_rec
# ###################################################
# # Load ggplot2
# library(ggplot2)
# # Make barplot
# ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
#     geom_bar(stat="identity", fill = "skyblue4") +
#     coord_cartesian(ylim=c(0.8, 1)) +
#     # coord_flip(ylim=c(0.8,1)) +
#     ylab("Cosine similarity\n original VS reconstructed") +
#     xlab("") +
#     # Reverse order of the samples such that first is up
#     # xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
#     theme_bw() +
#     theme(panel.grid.minor.y=element_blank(),
#           panel.grid.major.y=element_blank()) +
#     # Add cut.off line
#     geom_hline(aes(yintercept=.95))
# 

###################################################
### code chunk number 66: get_genes
###################################################
# For example get known genes table from UCSC for hg19 using 
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes_hg19


###################################################
### code chunk number 67: mut_strand
###################################################
strand = mut_strand(vcfs[[1]], genes_hg19)
head(strand, 10)


###################################################
### code chunk number 68: mut_matrix_stranded
###################################################
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
mut_mat_s[1:5,1:5]


###################################################
### code chunk number 69: strand_occurrences
###################################################
strand_counts <- strand_occurrences(mut_mat_s, by=tissue)
head(strand_counts)


###################################################
### code chunk number 70: strand_bias_test
###################################################
strand_bias <- strand_bias_test(strand_counts)
strand_bias


###################################################
### code chunk number 71: plot_strand
###################################################
ps1 <- plot_strand(strand_counts, mode = "relative")


###################################################
### code chunk number 72: plot_strand_bias_3
###################################################
ps2 <- plot_strand_bias(strand_bias)


###################################################
### code chunk number 73: plot_strand_bias_e (eval = FALSE)
###################################################
## grid.arrange(ps1, ps2)


###################################################
### code chunk number 74: plot_strand_bias
###################################################
# grid.arrange(ps1, ps2)
ggsave(paste("MutationalPatterns.plot_strand_bias", myvcf[1], "pdf", sep = "."),
       grid.arrange(ps1, ps2),
       width = 10,
       height = 8)

###################################################
### code chunk number 75: repli_file
###################################################
repli_file = system.file("extdata/ReplicationDirectionRegions.bed",
                         package = "MutationalPatterns")
repli_strand = read.table(repli_file, header = TRUE)
# Store in GRanges object
repli_strand_granges = GRanges(seqnames = repli_strand$Chr,
                               ranges = IRanges(start = repli_strand$Start + 1,
                                                end = repli_strand$Stop),
                               strand_info = repli_strand$Class)
# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) = "UCSC"
repli_strand_granges


###################################################
### code chunk number 76: strand_from_vcf_rep
###################################################
strand_rep <- mut_strand(vcfs[[1]], repli_strand_granges, mode = "replication")
head(strand_rep, 10)


###################################################
### code chunk number 77: mut_matrix_stranded_rep
###################################################
mut_mat_s_rep <- mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges,
                                     mode = "replication")
mut_mat_s_rep[1:5, 1:5]


###################################################
### code chunk number 78: specify_levels
###################################################
repli_strand_granges$strand_info <- factor(repli_strand_granges$strand_info,
                                           levels = c("right", "left"))
mut_mat_s_rep2 <- mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges,
                                      mode = "replication")
mut_mat_s_rep2[1:5, 1:5]


###################################################
### code chunk number 79: strand_occurrences_rep
###################################################
strand_counts_rep <- strand_occurrences(mut_mat_s_rep, by=tissue)
head(strand_counts)


###################################################
### code chunk number 80: strand_bias_test_rep
###################################################
strand_bias_rep <- strand_bias_test(strand_counts_rep)
strand_bias_rep


###################################################
### code chunk number 81: plot_strand_rep
###################################################
ps1 <- plot_strand(strand_counts_rep, mode = "relative")


###################################################
### code chunk number 82: plot_strand_bias_rep2
###################################################
ps2 <- plot_strand_bias(strand_bias_rep)


###################################################
### code chunk number 83: plot_strand_bias_rep_e (eval = FALSE)
###################################################
## grid.arrange(ps1, ps2)


###################################################
### code chunk number 84: plot_strand_bias_rep
###################################################
# grid.arrange(ps1, ps2)
ggsave(paste("MutationalPatterns.plot_strand_bias.2", myvcf[1], "pdf", sep = "."),
       grid.arrange(ps1, ps2),
       width = 10,
       height = 8)

###################################################
### code chunk number 85: extract_signatures
###################################################
nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)

# Provide signature names
colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")


###################################################
### code chunk number 86: plot_192
###################################################
a <- plot_192_profile(nmf_res_strand$signatures, condensed = TRUE)


###################################################
### code chunk number 87: plot_strand_bias
###################################################
b <- plot_signature_strand_bias(nmf_res_strand$signatures)


###################################################
### code chunk number 88: plot_192_profile_noeval (eval = FALSE)
###################################################
## grid.arrange(a, b, ncol = 2, widths = c(5, 1.8))


###################################################
### code chunk number 89: plot_192_profile
###################################################
ggsave(paste("MutationalPatterns.plot_192_profile", myvcf[1], "pdf", sep = "."),
       grid.arrange(a, b, ncol = 2, widths = c(5, 2)),
       width = 12,
       height = 5)


###################################################
### code chunk number 90: plot_rainfall_noeval (eval = FALSE)
###################################################
## # Define autosomal chromosomes
## chromosomes <- seqnames(get(ref_genome))[1:22]
## 
## # Make a rainfall plot
## plot_rainfall(vcfs[[1]], title = names(vcfs[1]),
##                 chromosomes = chromosomes, cex = 1.5, ylim = 1e+09)


###################################################
### code chunk number 91: plot_rainfall
###################################################
# Define autosomal chromosomes
chromosomes <- seqnames(get(ref_genome))[1:22]

# Make a rainfall plot
ggsave(paste("MutationalPatterns.plot_rainfall", myvcf[1], "pdf", sep = "."),
       plot_rainfall(vcfs[[1]], title = names(vcfs[1]),
                     chromosomes = chromosomes, cex = 1.5, ylim = 1e+09),
       width = 9,
       height = 3)
ggsave(paste("MutationalPatterns.plot_rainfall", myvcf[2], "pdf", sep = "."),
       plot_rainfall(vcfs[[2]], title = names(vcfs[2]),
                     chromosomes = chromosomes, cex = 1.5, ylim = 1e+09),
       width = 9,
       height = 3)
ggsave(paste("MutationalPatterns.plot_rainfall", myvcf[3], "pdf", sep = "."),
       plot_rainfall(vcfs[[3]], title = names(vcfs[3]),
                     chromosomes = chromosomes, cex = 1.5, ylim = 1e+09),
       width = 9,
       height = 3)
ggsave(paste("MutationalPatterns.plot_rainfall", myvcf[4], "pdf", sep = "."),
       plot_rainfall(vcfs[[4]], title = names(vcfs[4]),
                     chromosomes = chromosomes, cex = 1.5, ylim = 1e+09),
       width = 9,
       height = 3)
ggsave(paste("MutationalPatterns.plot_rainfall", myvcf[5], "pdf", sep = "."),
       plot_rainfall(vcfs[[5]], title = names(vcfs[5]),
                     chromosomes = chromosomes, cex = 1.5, ylim = 1e+09),
       width = 9,
       height = 3)
ggsave(paste("MutationalPatterns.plot_rainfall", myvcf[6], "pdf", sep = "."),
       plot_rainfall(vcfs[[6]], title = names(vcfs[6]),
                     chromosomes = chromosomes, cex = 1.5, ylim = 1e+09),
       width = 9,
       height = 3)


# ###################################################
# ### code chunk number 92: install_biomaRt (eval = FALSE)
# ###################################################
# ## if (!requireNamespace("BiocManager", quietly=TRUE))
# ##     install.packages("BiocManager")
# ## BiocManager::install("biomaRt")
# 
# 
# ###################################################
# ### code chunk number 93: load_biomart
# ###################################################
# library(biomaRt)
# 
# 
# ###################################################
# ### code chunk number 94: download_using_biomaRt
# ###################################################
# # regulatory <- useEnsembl(biomart="regulation",
# #                          dataset="hsapiens_regulatory_feature",
# #                          GRCh = 37)
# 
# ## Download the regulatory CTCF binding sites and convert them to
# ## a GRanges object.
# # CTCF <- getBM(attributes = c('chromosome_name',
# #                             'chromosome_start',
# #                             'chromosome_end',
# #                             'feature_type_name',
# #                             'cell_type_name'),
# #              filters = "regulatory_feature_type_name", 
# #              values = "CTCF Binding Site", 
# #              mart = regulatory)
# #
# # CTCF_g <- reduce(GRanges(CTCF$chromosome_name,
# #                 IRanges(CTCF$chromosome_start,
# #                 CTCF$chromosome_end)))
# 
# CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
#                               package="MutationalPatterns"))
# 
# ## Download the promoter regions and convert them to a GRanges object.
# 
# # promoter = getBM(attributes = c('chromosome_name', 'chromosome_start',
# #                                 'chromosome_end', 'feature_type_name'),
# #                  filters = "regulatory_feature_type_name", 
# #                  values = "Promoter", 
# #                  mart = regulatory)
# # promoter_g = reduce(GRanges(promoter$chromosome_name,
# #                     IRanges(promoter$chromosome_start,
# #                             promoter$chromosome_end)))
# 
# promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
#                                   package="MutationalPatterns"))
# 
# ## Download the promoter flanking regions and convert them to a GRanges object.
# 
# # flanking = getBM(attributes = c('chromosome_name',
# #                                 'chromosome_start',
# #                                 'chromosome_end',
# #                                 'feature_type_name'),
# #                  filters = "regulatory_feature_type_name", 
# #                  values = "Promoter Flanking Region", 
# #                  mart = regulatory)
# # flanking_g = reduce(GRanges(
# #                        flanking$chromosome_name,
# #                        IRanges(flanking$chromosome_start,
# #                        flanking$chromosome_end)))
# 
# flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
#                                   package="MutationalPatterns"))
# 
# 
# ###################################################
# ### code chunk number 95: combine_genomic_regions
# ###################################################
# regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
# 
# names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
# 
# 
# ###################################################
# ### code chunk number 96: combine_genomic_regions_2
# ###################################################
# seqlevelsStyle(regions) <- "UCSC"
# 
# 
# ###################################################
# ### code chunk number 97: download_bed_data
# ###################################################
# ## Get the filename with surveyed/callable regions
# surveyed_file <- system.file("extdata/callableloci-sample.bed",
#                              package = "MutationalPatterns")
# 
# ## Import the file using rtracklayer and use the UCSC naming standard
# library(rtracklayer)
# surveyed <- import(surveyed_file)
# seqlevelsStyle(surveyed) <- "UCSC"
# 
# ## For this example we use the same surveyed file for each sample.
# surveyed_list <- rep(list(surveyed), 9)
# 
# 
# ###################################################
# ### code chunk number 98: genomic_distribution
# ###################################################
# ## Calculate the number of observed and expected number of mutations in
# ## each genomic regions for each sample.
# distr <- genomic_distribution(vcfs, surveyed_list, regions)
# 
# 
# ###################################################
# ### code chunk number 99: enrichment_depletion_test
# ###################################################
# ## Perform the enrichment/depletion test by tissue type.
# distr_test <- enrichment_depletion_test(distr, by = tissue)
# head(distr_test)
# 
# 
# ###################################################
# ### code chunk number 100: plot_enrichment_depletion_e (eval = FALSE)
# ###################################################
# ## plot_enrichment_depletion(distr_test)
# 
# 
# ###################################################
# ### code chunk number 101: plot_enrichment_depletion
# ###################################################
# plot_enrichment_depletion(distr_test)
# 
# 
# ###################################################
# ### code chunk number 102: sessionInfo
# ###################################################
# toLatex(sessionInfo())

