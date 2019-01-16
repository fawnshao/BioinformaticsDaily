# whole mut, highTIP Mut, lowTIP Mut, whle eMut, highTIP eMut, lowTIP eMut
myvcf <- commandArgs(TRUE)
# myvcf <- "hg19.BRCA-US.enhancer.mut.vcf"
# myvcf <- c("hg19.BRCA-US.all.mut.vcf", "hg19.BRCA-US.all.highTIP60.mut.vcf", "hg19.BRCA-US.all.lowTIP60.mut.vcf", 
#            "hg19.BRCA-US.enhancer.mut.vcf", "hg19.BRCA-US.enhancer.highTIP60.mut.vcf", 
#            "hg19.BRCA-US.enhancer.lowTIP60.mut.vcf")

## ----results='hide', echo=FALSE, message=FALSE, warning=FALSE-------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("SomaticSignatures", version = "3.8")
set.seed(1)

options(width = 70)

library(knitr)

inlineCode <- function(file, format = c) {
    file_exist = sapply(file, file.exists)
    file = file[file_exist]
    if(length(file) == 0)
        return("")
    style = sapply(file,
                   function(file) {
                       paste(readLines(file), collapse = "\n")},
                   USE.NAMES = FALSE)
    style = sapply(style, format)
    style = paste(style, "\n", collapse = "\n\n")
    return(style)
}

knitrHeader <- function(css, js) {
    header = opts_knit$get("header")
    if(!missing(css) && !identical(css, character())) {
        header["highlight"] = inlineCode(css)
    }
    if(!missing(js) && !identical(js, character())) {
        header["js"] = inlineCode(js, formatInlineJS)
    }
    return(header)
}

base_dir = system.file(package = "SomaticSignatures")
css_path = file.path(base_dir, "css", "bioc.css")

opts_knit$set(self.contained = TRUE,
              upload.fun = image_uri,
              header = knitrHeader(css = css_path))

opts_chunk$set(comment = "  ",
               fig.path = "",
               fig.align = "center",
               out.width = "65%",
               dpi = 300,
               indent = 10,
               cache = FALSE,
               cache.path = "../cache")

knit_hooks$set(fig.cap = function(before, options, envir) {
    if(!before) {
        paste('<p class="caption">',options$fig.cap,"</p>",sep="")
    }
})

## ----load_ss, results='hide',message=FALSE--------------------------
library(SomaticSignatures)

## ----load_data_package, results='hide',message=FALSE----------------
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

## ----sca_metadata---------------------------------------------------
sca_metadata = scaMetadata()

sca_metadata

## ----sca_to_vranges-------------------------------------------------
sca_data = unlist(scaLoadDatasets())

sca_data$study = factor(gsub("(.*)_(.*)", "\\1", toupper(names(sca_data))))
sca_data = unname(subset(sca_data, Variant_Type %in% "SNP"))
sca_data = keepSeqlevels(sca_data, hsAutosomes(), pruning.mode = "coarse")

sca_vr = VRanges(
    seqnames = seqnames(sca_data),
    ranges = ranges(sca_data), 
    ref = sca_data$Reference_Allele,
    alt = sca_data$Tumor_Seq_Allele2, 
    sampleNames = sca_data$Patient_ID,
    seqinfo = seqinfo(sca_data), 
    study = sca_data$study)

sca_vr

## ----sca_study_table------------------------------------------------
sort(table(sca_vr$study), decreasing = TRUE)

## ----sca_vr_to_motifs-----------------------------------------------
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.1000genomes.hs37d5)
head(sca_motifs)

## ----sca_motif_occurrence-------------------------------------------
sca_mm = motifMatrix(sca_motifs, group = "study", normalize = TRUE)

head(round(sca_mm, 4))

## ----sca_mutation_spectrum, fig.cap='Mutation spectrum over studies'----
plotMutationSpectrum(sca_motifs, "study")

## ----sca_nmf_pca----------------------------------------------------
n_sigs = 5

sigs_nmf = identifySignatures(sca_mm, n_sigs, nmfDecomposition)

sigs_pca = identifySignatures(sca_mm, n_sigs, pcaDecomposition)

## ----sca_explore_nmf------------------------------------------------
sigs_nmf

## ----sca_explore_pca------------------------------------------------
sigs_pca

## -------------------------------------------------------------------
n_sigs = 2:8

gof_nmf = assessNumberSignatures(sca_mm, n_sigs, nReplicates = 5)

gof_pca = assessNumberSignatures(sca_mm, n_sigs, pcaDecomposition)

## ----fig.cap='Summary statistics for selecting the number of signatures in the NMF decomposition.'----
plotNumberSignatures(gof_nmf)

## ----fig.cap='Summary statistics for selecting the number of signatures in the PCA decomposition.'----
plotNumberSignatures(gof_pca)

## ----load_ggplot2, results='hide',message=FALSE---------------------
library(ggplot2)

## ----sca_plot_nmf_signatures_map, fig.cap='Composition of somatic signatures estimated with the NMF, represented as a heatmap.'----
plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")

## ----sca_plot_nmf_signatures, fig.cap='Composition of somatic signatures estimated with the NMF, represented as a barchart.'----
plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Barchart")

## -------------------------------------------------------------------
plotObservedSpectrum(sigs_nmf)

## -------------------------------------------------------------------
plotFittedSpectrum(sigs_nmf)

## ----sca_plot_nmf_samples_map, fig.cap='Occurrence of signatures estimated with the NMF, represented as a heatmap.'----
plotSampleMap(sigs_nmf)

## ----sca_plot_nmf_samples, fig.cap='Occurrence of signatures estimated with the NMF, represented as a barchart.'----
plotSamples(sigs_nmf)

## ----sca_plot_pca_signatures_map, fig.cap='Composition of somatic signatures estimated with the PCA, represented as a heatmap.'----
plotSignatureMap(sigs_pca) + ggtitle("Somatic Signatures: PCA - Heatmap")

## ----sca_plot_pca_signatures, fig.cap='Composition of somatic signatures estimated with the PCA, represented as a barchart.'----
plotSignatures(sigs_pca) + ggtitle("Somatic Signatures: PCA - Barchart")

## -------------------------------------------------------------------
plotFittedSpectrum(sigs_pca)

## -------------------------------------------------------------------
plotObservedSpectrum(sigs_pca)

## ----load_ggplot2_again, results='hide',message=FALSE---------------
library(ggplot2)

## ----sca_plot_nmf_samples_mod, results='hide',message=FALSE---------
p = plotSamples(sigs_nmf)

## (re)move the legend
p = p + theme(legend.position = "none")
## (re)label the axis
p = p + xlab("Studies")
## add a title
p = p + ggtitle("Somatic Signatures in TGCA WES Data")
## change the color scale
p = p + scale_fill_brewer(palette = "Blues")
## decrease the size of x-axis labels
p = p + theme(axis.text.x = element_text(size = 9))

## ----sca_plot_nmf_samples_mod_print, fig.cap='Occurrence of signatures estimated with the NMF, customized plot. See the original plot above for comparisons.'----
p

## -------------------------------------------------------------------
clu_motif = clusterSpectrum(sca_mm, "motif")

## ----fig.cap='Hierachical clustering of the mutational spectrum, according to motif.'----
library(ggdendro)

p = ggdendrogram(clu_motif, rotate = TRUE)
p

## ----sva_load, results='hide',message=FALSE-------------------------
library(sva)

## ----sva_batch------------------------------------------------------
sca_anno = as.data.frame(lapply(sca_metadata, unlist))

model_null = model.matrix(~ 1, sca_anno)

sca_mm_batch = ComBat(sca_mm, batch = sca_anno$Sequence_Source, mod = model_null)

## ----kmer_hs_chrs, results='hide',message=FALSE---------------------
k = 3
n = 1e4

hs_chrs = as(seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5), "GRanges")
hs_chrs = keepSeqlevels(hs_chrs, c(1:22, "X", "Y"), pruning.mode = "coarse")

k3_hs_chrs = kmerFrequency(BSgenome.Hsapiens.1000genomes.hs37d5, n, k, hs_chrs)
k3_hs_chrs

## ----kmer_exons, eval=FALSE-----------------------------------------
#  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#  
#  k = 3
#  n = 1e4
#  
#  hs_exons = reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
#  hs_exons = ncbi(keepStandardChromosomes(hs_exons))
#  
#  k3_exons = kmerFrequency(BSgenome.Hsapiens.1000genomes.hs37d5, n, k, hs_exons)

## ----normalize_motifs-----------------------------------------------
data(kmers)
norms = k3wg / k3we
head(norms)

sca_mm_norm = normalizeMotifs(sca_mm, norms)

## -------------------------------------------------------------------
showMethods("getSeq")

## ----eval=FALSE-----------------------------------------------------
#  ## Somatic variant calls
#  vr_A = readVcfAsVRanges(vcf_A_path, "GenomeA")
#  vr_B = readVcfAsVRanges(vcf_B_path, "GenomeB")
#  
#  ## Genomic sequences
#  fa_A = FaFile(fasta_A_path)
#  fa_B = FaFile(fasta_B_path)
#  
#  ## Somatic motifs
#  vr_A = mutationContext(vr_A, fa_A)
#  vr_B = mutationContext(vr_B, fa_B)
#  
#  ## Combine for further analysis
#  vr = c(vr_A, vr_B)

## -------------------------------------------------------------------
citation("SomaticSignatures")

## ----eval=FALSE-----------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("SomaticSignatures")

## ----eval=FALSE-----------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite()

## ----results='hide', message=FALSE----------------------------------
library(VariantAnnotation)

## -------------------------------------------------------------------
vr = VRanges(
    seqnames = "chr1",
    ranges = IRanges(start = 1000, width = 1),
    ref = "A",
    alt = "C")

vr

## ----echo=FALSE, results='markup'-----------------------------------
sessionInfo()
