library(monocle3)
library(dplyr)
setwd("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/NS27.NUP93.scRNA")
cds1 <- load_cellranger_data(pipestance_path = "X1", genome = "hg19")
cds2 <- load_cellranger_data(pipestance_path = "X2", genome = "hg19")
cds <- new_cell_data_set(,
                     cell_metadata = cell_metadata,
                     gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, reduction_method = "UMAP") 
cds <- learn_graph(cds)
plot_cells(cds)
plot_cells(cds, genes = c("MYC", "KITLG", "NUP93"))


#################
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

# Make the CDS object
testcds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

##################
# test
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
