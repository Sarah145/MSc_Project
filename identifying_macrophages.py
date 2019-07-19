# load libraries
import numpy as np
import scipy as sp
import scanpy as sc
import pandas as pd
from bbknn import bbknn
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150)
sc.logging.print_version_and_date()

# read in filtered data
adata_filt = sc.read_h5ad('filtered_human_data.h5ad')

# get list of tissues
tissues = list(set(adata_filt.obs['Tissue']))
tissues.sort()

# identify macrophages in each tissue individually
counter = 0
mac_codes = list()
adata=adata_filt[adata_filt.obs.Tissue==tissues[counter]]
print(adata.obs.Tissue[0], "with", adata.n_obs, "cells from", len(set(adata.obs.Dataset)), "datasets")

# filter out uninformative genes to reduce size of dataset
sc.pp.filter_genes(adata, min_cells=3)

# normalise to CPM and transform to log scale and save as raw
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata

# find highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=3000)

# subset to only highly variable genes
adata = adata[:, adata.var['highly_variable']]

# compute pca
sc.tl.pca(adata, svd_solver='arpack')

# compute connectivities - if tissue contains more than 1 dataset use bbknn - if not use neighbors
# bbknn(adata, neighbors_within_batch = 3, trim = 100, n_pcs = 30, batch_key = 'Dataset')
sc.pp.neighbors(adata, n_pcs=30)

# compute UMAP coordinates
sc.tl.umap(adata)

# cluster
sc.tl.leiden(adata, resolution=0.5)

# visualise
sc.pl.umap(adata, color=['leiden'], legend_loc='on data')

# look for clusters expressing known macrophage markers
sc.pl.umap(adata, color=['CD68', 'CD14', 'CD163', 'C1QA', 'FCGR3A', 'LYZ'], ncols = 3)
sc.pl.violin(adata, ['CD68', 'CD14', 'CD163', 'C1QA', 'FCGR3A', 'LYZ'], groupby = 'leiden', rotation = 90)

# take subset with only macrophages
tmp = adata[(adata.obs.leiden== 'macrophage_cluster' ) ]
print(adata.obs.Tissue[0], "contains ", tmp.n_obs, "macrophages which is ", (tmp.n_obs/adata.n_obs)*100, "% of the cells in this tissue.")

# save barcodes of macrophages and repeat for each tissue
mac_codes += list(tmp.obs.Barcodes)
counter += 1

# take subset of adata_filt with only high quality cells and write out to macrophages.h5ad
adata = adata_filt[mac_codes, :]
adata.write("macrophages.h5ad")





























