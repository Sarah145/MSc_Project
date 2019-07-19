# load libraries
import numpy as np
import scipy as sp
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
%load_ext rpy2.ipython

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150)
sc.logging.print_version_and_date()

# read in macrophage data
adata = sc.read_h5ad('macrophages.h5ad')
aliases = adata.var_names.to_list()


# use R to convert all genes to official symbols
%%R  -i aliases -o symbols 
library(limma)
symbols <- alias2SymbolTable(aliases, species='Hs')

# make dataframe with aliases and symbols, remove genes with no symbol
genes_to_symbols = pd.DataFrame(zip(aliases, symbols), columns = ['alias', 'symbol'])
no_symbol = genes_to_symbols['alias'][genes_to_symbols['symbol']== 'NA']
genes_to_symbols = genes_to_symbols[genes_to_symbols['symbol'] != 'NA']

# get list of datasets
data = list(set(adata_raw.obs['Dataset']))
data.sort()

# take subset with just one dataset
tmp = adata[adata.obs['Dataset']==data[0]]

# filter out genes that are not expressed in this dataset
sc.pp.filter_genes(tmp, min_cells=1)

# get aliases
alias = tmp.var_names

# get symbols if they're available
symbol = list()
keep_aliases = list()
for i in alias:
    if i not in no_symbol:
        symbol.append(genes_to_symbols.loc[i,'symbol'])
        keep_aliases.append(i)
        
# subset to only genes with symbols available
tmp = tmp[:,keep_aliases]

# change aliases to symbols and save object as adata_con
tmp.var_names = symbol
adata_con = tmp

# loop over all other datasets and concatenate into one object
for i in range(1, len(data)+1):
    tmp = adata[adata.obs['Dataset']==data[i]]
    sc.pp.filter_genes(tmp, min_cells=1)
    aliases = tmp.var_names
    symbols = list()
    keep_aliases = list()
    for i in aliases:
        if i not in no_symbol:
            symbols.append(genes_to_symbols.loc[i,'symbol'])
            keep_aliases.append(i)
    tmp = tmp[:,keep_aliases]
    tmp.var_names = symbols
    adata_con = sc.AnnData.concatenate(adata_con, tmp,  join='outer')
    
# write out converted adata object
adata_con.write("macrophages_con.h5ad")


