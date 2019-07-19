# load libraries
import numpy as np
import scipy as sp
import scanpy as sc
import pandas as pd
from gprofiler import gprofiler
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150)
sc.logging.print_version_and_date()

# read in data
adata = sc.read_h5ad('adata_final.h5ad')

# create custom background gene list
bg = adata.raw.var_names.tolist()

# read in tissue DEGs calculated in R
de_genes_tissues = pd.read_csv("../de_genes_tissues.csv", sep='\t')

# perform functional enrichment for each tissue
for i in list(set(de_genes_tissues['cluster'])):
    print("Calculating", i)
    tmp = de_genes_tissues[de_genes_tissues['cluster']==i]
    tmp = tmp[tmp['avg_logFC'] > 1] # take only upregulated genes
    tmp = tmp[tmp['p_val_adj'] <= 0.05]
    # find GO terms
    enrichment = gprofiler(tmp['gene'], custom_bg=bg , organism='hsapiens', correction_method='fdr', src_filter=['GO:BP'])
    try:     # if enriched terms are found, write them to a csv
    	enrichment.sort_values('p.value').iloc[:,[2,3,4,5,6,11,13]].to_csv("BP_enrichment/tissues/" + str(tmp['cluster'].iloc[0]) + "_enrichment_BP.csv")
    # find Reactome pathways
    enrichment = gprofiler(tmp['gene'], custom_bg=bg , organism='hsapiens', correction_method='fdr', src_filter=['REAC'])
    try:
	    enrichment.sort_values('p.value').iloc[:,[2,3,4,5,6,11,13]].to_csv("REAC_enrichment/tissues/" + str(tmp['cluster'].iloc[0]) + "_enrichment_REAC.csv")
 
 # read in cluster DEGs calculated in R
de_genes_tissues = pd.read_csv("../de_genes_clusters.csv", sep='\t')

# perform functional enrichment for each cluster
for i in list(set(de_genes_clusters['cluster'])):
    print("Calculating", i)
    tmp = de_genes_clusters[de_genes_clusters['cluster']==i]
    tmp = tmp[tmp['avg_logFC'] > 1] # take only upregulated genes
    tmp = tmp[tmp['p_val_adj'] <= 0.05]
    # find GO terms
    enrichment = gprofiler(tmp['gene'], custom_bg=bg , organism='hsapiens', correction_method='fdr', src_filter=['GO:BP'])
    try:     # if enriched terms are found, write them to a csv
    	enrichment.sort_values('p.value').iloc[:,[2,3,4,5,6,11,13]].to_csv("BP_enrichment/clusters/" + str(tmp['cluster'].iloc[0]) + "_enrichment_BP.csv")
    # find Reactome pathways
    enrichment = gprofiler(tmp['gene'], custom_bg=bg , organism='hsapiens', correction_method='fdr', src_filter=['REAC'])
    try:
	    enrichment.sort_values('p.value').iloc[:,[2,3,4,5,6,11,13]].to_csv("REAC_enrichment/clusters/" + str(tmp['cluster'].iloc[0]) + "_enrichment_REAC.csv")
 