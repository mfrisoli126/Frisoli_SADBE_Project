#############################################################
### Find Markers
#############################################################
library(Seurat)
library(SeuratObject)
library(SignallingSingleCell)
library(tibble)
library(dplyr)
library(umap)
library(FNN)
library(igraph)
library(stringr)
library(tidyr)
library(clusterProfiler)
library(magrittr)
library(aliases2entrez)
library(data.table)
library(Matrix)
library(ggrepel)


setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load data:
load(file = "Data/contact_derm_filtered3.Rdata")

DefaultAssay(contact_derm) <- "RNA"

all.markers <- FindAllMarkers(contact_derm, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all.markers, file = "Data/contact_derm_filtered.markers.tsv" , row.names=FALSE, sep="\t")
