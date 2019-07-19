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
adata = sc.read_h5ad('macrophages_con.h5ad')

# filter out uninformative genes to reduce size of dataset
sc.pp.filter_genes(adata, min_cells=5)

# remove mitochondrial genes because they've already been removed from some datasets
gene_list = list()
for i in adata.var_names:
    if not i.startswith('MT-'):
        gene_list.append(i)
adata = adata[:, gene_list]

# normalise to CPM and transform to log scale and save as raw
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata

# find highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=4000)

# subset to only highly variable genes
adata = adata[:, adata.var['highly_variable']]

# compute pca
sc.tl.pca(adata, svd_solver='arpack')

# find neighbors and cluster (without batch alignment)
nbr = sc.pp.neighbors(adata, n_pcs=30, copy=True)
sc.tl.umap(nbr)
sc.tl.leiden(nbr)
sc.pl.umap(nbr, color=['leiden'], legend_loc='on data')
sc.pl.umap(nbr, color=['Tissue'])
sc.pl.umap(nbr, color=['Dataset'])

# find neighbors and cluster (with batch alignment)
bbknn(bb, n_pcs=30, batch_key='Protocol', neighbors_within_batch=6, trim=40)
sc.tl.umap(bb, min_dist=0.6, spread=1.4)
sc.tl.leiden(bb, resolution = 0.27)
sc.pl.umap(bb, color=['leiden'], legend_loc='on data', size=3.5)
sc.pl.umap(bb, color=['Tissue'], legend_loc='on data',size = 3.5)
sc.pl.umap(bb, color=['Dataset'], size=3.5)

# change colors
bb.uns['leiden_colors'] = [ 
	'#FF230A', #1
	'#BA2223', #2
	'#993636', #3
	'#7D4D43', #4
	'#DC5BFC', #5
	'#DBD767', #6
	'#98df8a', #7
	'#279e68', #8
	'#17becf', #9
	'#aec7e8', #10
	'#005EFF', #11
	'#FF70DE', #12
	'#ffbb78', #13
	'#c5b0d5', #14
	'#65069C', #15
	'#f7b6d2' #16                               
]
bb.uns['Tissue_colors'] = [
     '#d62728', # blood
     '#279e68', # brain
     '#964B00', # colon
     '#ff00ef', # decidua
     '#aa40fc', # fetal liver
     '#E8E848', # fetal skin
     '#aec7e8', # kidney
     '#1406D6', # lung
     '#b5bd61', # microglia
     '#ff7f0e', # pancreas
     '#ff9896', # placenta
     '#00fff3', # prostate
     '#00ff89', # testis
     '#f7b6d2', # upper airway
     '#ff7f0e'] # mLN

# write out final adata object
bb.write("adata_final.h5ad")

# write out cluster/tissue info for plotting R
clusters = adata.obs['leiden1']
tissues = adata.obs['Tissue']
pd.DataFrame(zip(clusters, tissues), columns = ["clusters", "tissues"]).to_csv("clus_tis_df.csv", index = False)

# write out data to analyse in R
adata = bb

# write mtx file
adataT = adata.raw.X.transpose()
sp.io.mmwrite("../project_R/matrix.mtx", adataT)

# write barcodes to tsv
j = list()
for i in adata.obs.Barcode.values:
    k = i.split("-", 1)[1]
    j.append(k)
h = list()
for i in adata.obs.Barcode.values:
    k = i.split("-", 1)[0]
    h.append(k)
barcodes = pd.DataFrame(h)
barcodes['b'] = pd.Series(j, index=barcodes.index)
barcodes.index = list(range(1,len(barcodes)+1))
barcodes.to_csv('../project_R/barcodes.tsv', sep='\t', header=False)

# write genes to tsv
j= adata.raw.var_names.to_list()
genes = pd.DataFrame(j)
genes.index = list(range(1,len(genes)+1))
genes.to_csv('../project_R/genes.tsv', sep='\t',  header=False)

# write metadata to tsv
j = list(adata.obs['Tissue'])
k = list(adata.obs['Dataset'])
l = list(adata.obs['leiden'])
meta = pd.DataFrame(list(zip(j,k,l)), columns=['Tissue', 'Dataset', 'leiden'])
meta.index = list(range(1,len(meta)+1))
meta.to_csv('../project_R/meta.tsv', sep='\t', index=False)
