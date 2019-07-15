# load libraries
import numpy as np
import scipy as sp
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150)
sc.logging.print_version_and_date()


# read in raw data
adata_raw = sc.read_h5ad('raw_human_data.h5ad')

# calculate qc metrics for each cell
adata_raw.obs['n_counts'] = adata_raw.X.sum(1)
adata_raw.obs['n_genes'] = (adata_raw.X > 0).sum(1)
mito_genes = adata_raw.var_names.str.startswith('MT-')
adata_raw.obs['percent_mito'] = np.sum(adata_raw[:, mito_genes].X, axis=1).A1 / np.sum(adata_raw.X, axis=1).A1

# visualise distributions of QC metrics for different datasets
sc.pl.violin(adata_raw, 'n_counts', groupby='Dataset', size=2, log=True, cut=0, rotation=90)
sc.pl.violin(adata_raw, 'n_genes', groupby='Dataset', size=2, log=True, cut=0, rotation=90)
sc.pl.violin(adata_raw, 'percent_mito', groupby='Dataset', rotation=90)

# get list of datasets
data = list(set(adata_raw.obs['Dataset']))
data.sort()

# perform QC on each dataset individually
counter = 0
filtered_codes = list()
adata=adata_raw[adata_raw.obs.Dataset==data[counter]]
print(adata.obs.Dataset[0], "with", adata.n_obs, "cells from", set(adata.obs.Tissue))

# visualise QC metrics for this dataset
sc.pl.scatter(adata, 'n_counts', 'n_genes', color='percent_mito')

# if QC information for this dataset is available from the original publication then apply the same thresholds here 
# if not - cut out cells in top and bottom percentile for n_counts and n_genes and top percentile for percent_mito
high_count_thresh = np.percentile(adata.obs['n_counts'], 99)
low_count_thresh = np.percentile(adata.obs['n_counts'], 1)

high_gene_thresh = np.percentile(adata.obs['n_genes'], 99)
low_gene_thresh = np.percentile(adata.obs['n_genes'], 1)

mito_thresh = np.percentile(adata.obs['percent_mito'], 98)

adata = adata[(adata.obs['n_counts'] <= high_count_thresh)]
adata = adata[(adata.obs['n_genes'] <= high_gene_thresh)]
adata = adata[adata.obs['percent_mito'] <= mito_thresh]

print(adata.n_obs, "cells left after filtering")

# visualise again
sc.pl.scatter(adata, 'n_counts', 'n_genes', color='percent_mito')

# save barcodes of high quality cells and repeat for each dataset
filtered_codes += list(adata.obs['Barcodes'])
counter += 1

# take subset of adata_raw with only high quality cells and write out to filtered_human_data.h5ad
adata = adata_raw[filtered_codes, :]
adata.write("filtered_human_data.h5ad")


