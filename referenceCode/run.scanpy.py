import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.logging.print_versions()
results_file = './test.h5ad'
# sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
# adata = sc.read_10x_mtx(
#     'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
#     var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
#     cache=True)                              # write a cache file for faster subsequent reading

mydata=sc.read("../ScienceDirect_files_19May2020_03-59-36.683/sim.GSE111113_Table_S1_FilterNormal10xExpMatrix.txt", delimiter='\t', cache=True)
sc.pp.recipe_zheng17(mydata)
sc.tl.pca(mydata, svd_solver='arpack')
sc.pp.neighbors(mydata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(mydata)
# sc.pl.draw_graph(mydata, color='paul15_clusters', legend_loc='on data')

sc.tl.umap(mydata)
sc.tl.leiden(mydata, resolution=1.0)
sc.tl.paga(mydata, groups='leiden')
sc.pl.paga(mydata, color=['leiden', 'Hba-a2', 'Elane', 'Irf8'])