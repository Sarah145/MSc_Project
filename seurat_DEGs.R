library(Seurat)
library(dplyr)
library(cowplot)
library(Matrix)
library(scales)
library(MAST)

setwd("~")

# Read in the dataset
counts <- Read10X(data.dir = "project_R")

# Create Seurat object
mac <- CreateSeuratObject(counts)

# Attach metadata
meta <- read.delim("project_R/meta.tsv", header=T, sep="\t")
mac@meta.data <- meta

Idents(mac) <- "Tissue"
mac_tissue.markers <- FindAllMarkers(mac, min.pct = 0.3, test.use = "MAST")
write.table(mac_tissue.markers, file = "de_genes_tissues.csv", col.names=T, row.names = F, sep="\t")

Idents(mac) <- "leiden"
mac_cluster.markers <- FindAllMarkers(mac, min.pct = 0.3, test.use = "MAST")
write.table(mac_cluster.markers, file = "de_genes_clusters.csv", col.names=T, row.names= F, sep="\t")

print("Finished!")
