#### MSc. Research Project

##### "Characterising transcriptional variation between tissue-resident macrophage subsets using single-cell RNA-seq"



This directory contains the scripts used to process data for my masters research project which involved a cross-tissue comparison of scRNA-seq data from macrophages in 15 different tissues. The final manuscript can be found [here](https://github.com/Sarah145/MSc_Project/blob/master/final_manuscript.pdf) and the macrophage dataset can be explored through an interactive web app [here](https://tiny.cc/cross-tissue-macrophages).



These scripts were run interactively (mainly from Jupyter notebooks) to analyze cells in the following order:

1. [`qc.py`](https://github.com/Sarah145/MSc_Project/blob/master/qc.py) - To remove low-quality cells in the raw dataset.
2. [`identifying_macrophages.py`](https://github.com/Sarah145/MSc_Project/blob/master/identifying_macrophages.py) - To identify macrophages among cells from each tissue.
3. [`converting_genes.py`](https://github.com/Sarah145/MSc_Project/blob/master/converting_genes.py) - To make sure all gene names were official symbols.
4. [`analysing_macrophages.py`](https://github.com/Sarah145/MSc_Project/blob/master/analysing_macrophages.py) - To analyse the macrophage dataset (normalisation, dimensionality reduction, visualisation).
5. [`seurat_DEGs.R`](https://github.com/Sarah145/MSc_Project/blob/master/seurat_DEGs.R) - To compute differential expression analysis with Seurat's implementation of MAST.
6. [`functional_enrichment.py`](https://github.com/Sarah145/MSc_Project/blob/master/functional_enrichment.py) - To find enrichment of Gene Ontology terms and Reactome pathways among differentially expressed genes.

-----

`Supplementary_table_2_extended.xls` contains the gene symbols, names, Uniprot descriptions, logFC and adj. p-value for the top 10 DE genes in each tissue. 

`plots.R` contains the ggplot scripts that were used to create some of the figures in the manuscript. 

