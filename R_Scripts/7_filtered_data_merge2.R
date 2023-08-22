##############################################################################################
### Merge filtered cells back to a single embedding and integrate for sequencer batch effect:
##############################################################################################
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

setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load data:
load(file = "Data/contact_derm_5_9_22.Rdata")

#Load filtered cell ID files:
cell.type.file.list <- list.files(path = "Data/Basic_Celltype_Cell_IDs/filtered2/", pattern = "_filtered2_Cell_IDs.Rdata")

all.cell.IDs <- c()

for (file in cell.type.file.list) {
  
  subset.name <- str_replace(file, pattern = "_Cell_IDs.Rdata", replacement = "")
  
  cell.IDs <- get(load(paste("Data/Basic_Celltype_Cell_IDs/filtered2/", subset.name, "_Cell_IDs.Rdata", sep = "")))
  
  all.cell.IDs <- c(all.cell.IDs, cell.IDs)
}

#Inspect filtered cells:
#contact_derm@meta.data$Final.cell <- colnames(contact_derm) %in% all.cell.IDs

#plot_tsne_metadata_srt(contact_derm, color_by = "Final.cell", size = 0.25, facet_by = c("Sequence_instrument"))
#plot_tsne_metadata_srt(contact_derm, color_by = "Final.cell", size = 0.25, colors = c("#F65937", "#5A5A5A"), plot_labels = FALSE, plot_legend = TRUE)

#plot_tsne_gene_srt(contact_derm, gene = c("KRT10", "KRT1", "KRT5"), size = 0.25)

#Filter by Cell IDs:
rownames(contact_derm@meta.data) <- colnames(contact_derm)

meta.data <- contact_derm@meta.data[which(colnames(contact_derm@assays$RNA@counts) %in% all.cell.IDs),]
counts <- contact_derm@assays$RNA@counts[,which(colnames(contact_derm@assays$RNA@counts) %in% all.cell.IDs)]
embedding <- contact_derm@reductions$umap

#Reprocess filtered dataset in Seurat with integration for sequencing instrument batch effect:
contact_derm <- CreateSeuratObject(counts = counts, project = "contact_derm", min.cells = 3, min.features = 150)
contact_derm@meta.data <- meta.data
contact_derm@reductions$umap<- embedding
contact_derm@reductions$umap@cell.embeddings <- contact_derm@reductions$umap@cell.embeddings[which(rownames(contact_derm@reductions$umap@cell.embeddings) %in% all.cell.IDs),]


#Try new 2022 Seurat integration method:
# split the dataset into a list of two seurat objects based on capture and alignment method:
dataset.list <- SplitObject(contact_derm, split.by = "Sequence_instrument")

# normalize and identify variable features for each dataset independently
dataset.list <- lapply(X = dataset.list, FUN = function(x) {
  
  x <- SCTransform(x, variable.features.n = 6000, ncells = 10000)
  
})

#Run integration:
integration.features <- SelectIntegrationFeatures(object.list = dataset.list, nfeatures = 6000)

#Remove exclusion genes before integration:
load("gene_lists/established_exclusion_genes_5.17.22.Rdata")

integration.features <- integration.features[-which(integration.features %in% established_exclusion_genes)]

#to integrate all features from all_genes, add them to the PrepSCTIntegration step, rerun the PrepSCTIntegration and FindIntegrationAnchors steps
dataset.list <- PrepSCTIntegration(object.list = dataset.list, anchor.features = integration.features)
anchors <- FindIntegrationAnchors(object.list = dataset.list, normalization.method = "SCT", anchor.features = integration.features)

#now IntegrateData takes the features.to.integrate 
contact_derm <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)

#Remove exclusion genes before PCA:
load("gene_lists/established_exclusion_genes_5.17.22.Rdata")

if (all(rownames(contact_derm)[which(str_detect(rownames(contact_derm), pattern = "^RPS|^RPL|^MT-"))] %in% established_exclusion_genes) == FALSE) {
  established_exclusion_genes <- unique(c(established_exclusion_genes, rownames(contact_derm)[which(str_detect(rownames(contact_derm), pattern = "^RPS|^RPL|^MT-"))]))
}

variable.features <- VariableFeatures(contact_derm)

if (any(variable.features %in% established_exclusion_genes)) {
  
  variable.features <- variable.features[-which(variable.features %in% established_exclusion_genes)]
}

if (length(variable.features) > 4000) {
  variable.features <- variable.features[1:4000]
}

VariableFeatures(contact_derm) <- variable.features

#PCA:
contact_derm <- RunPCA(contact_derm, verbose = FALSE, npcs = 75)

contact_derm <- FindNeighbors(contact_derm, dims = 1:75, k.param = 30)
contact_derm <- FindClusters(contact_derm, resolution = 0.5)

contact_derm <- RunUMAP(contact_derm, reduction = "pca", dims = 1:30, n.neighbors = 30L)

p1 <- DimPlot(contact_derm, reduction = "umap", group.by = "Sequence_instrument")
p2 <- DimPlot(contact_derm, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
              repel = TRUE)
p1 + p2



#Save Outputs:

save(contact_derm, file = "Data/contact_derm_filtered3.Rdata")


