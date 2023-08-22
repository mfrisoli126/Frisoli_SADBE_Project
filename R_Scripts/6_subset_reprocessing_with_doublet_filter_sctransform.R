#####################################################
#Re-processing Basic Celltype Subsets with Seurat:
#####################################################
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

#Load Data
load(file = "Data/contact_derm_5_9_22.Rdata")

#Define basic cell type specific markers:
MEL_specific_markers <- c("TYRP1", "MLANA", "PMEL", "MITF", "DCT")
KRT_specific_markers <- c("KRT10", "KRT14", "KRT17", "DMKN", "KRTDAP", "MUCL1")
LYMPH_specific_markers <- c("TRAC", "TRDC", "NCAM1")
APC_specific_markers <- c("CD207", "LYZ", "HLA.DR", "CD86", "CD83", "CLEC4A", "CLEC4C", "MRC1", "NRP1", "JCHAIN")
RBC_specific_markers <- c("HBB", "HBA2")

Basic_cell_specific_markers <- c(MEL_specific_markers, KRT_specific_markers, LYMPH_specific_markers, APC_specific_markers, RBC_specific_markers)



cell.type.file.list <- list.files(path = "Data/Basic_Celltype_Cell_IDs/initial/")

#1st Filter of each Basic Celltype Subset:

    for (file in cell.type.file.list) {
      
      subset.name <- str_replace(file, pattern = "_Cell_IDs.Rdata", replacement = "")
     
      cell.IDs <- get(load(paste("Data/Basic_Celltype_Cell_IDs/initial/", subset.name, "_Cell_IDs.Rdata", sep = "")))
      
      #Filter by Cell IDs:
      subset_data <- contact_derm@assays$RNA@counts[,which(colnames(contact_derm@assays$RNA@counts) %in% cell.IDs)]
      

      #Reprocess in Seurat:
            #Create Suerat Object, normalize, and scale:
            subset_data <- CreateSeuratObject(counts = subset_data, project = "contact_derm", min.cells = 3, min.features = 150)
            
            #Try SCTransform for processing:
            subset_data <- SCTransform(subset_data, variable.features.n = 4000, ncells = 10000)
            
            #PCA, Clustering, and UMAP:
            subset_data <- RunPCA(subset_data, npcs = 75)
            
            subset_data <- FindNeighbors(subset_data, dims = 1:75)
            subset_data <- FindClusters(subset_data, resolution = 0.5)
            subset_data <- RunUMAP(subset_data, dims = 1:75)
            
            DimPlot(subset_data, reduction = "umap", label = TRUE, pt.size = 0.75) + NoLegend()
            
            all.markers <- FindAllMarkers(subset_data, min.pct = 0.25, logfc.threshold = 0.25)
                
            
                  # Filter clusters that appear to be doublets:
                  unique.cluster.labels <- as.character(unique(subset_data@meta.data$seurat_clusters))
                  
                  clusters.to.remove <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>%
                    filter(gene %in% Basic_cell_specific_markers[-which(Basic_cell_specific_markers %in% eval(as.symbol(paste(subset.name,"_specific_markers", sep = ""))))]) %>% 
                    dplyr::select(cluster) %>% unique()
                  clusters.to.remove <- as.vector(clusters.to.remove$cluster)
                  
                  if (length(clusters.to.remove) == 0) {
                    clusters.to.keep <- unique.cluster.labels
                  }
                  
                  if (length(clusters.to.remove) != 0) {
                    clusters.to.keep <- unique.cluster.labels[-which(unique.cluster.labels %in% clusters.to.remove)]
                  }
                  
                  cell.IDs.to.keep <- colnames(subset_data)[which(subset_data@meta.data$seurat_clusters %in% clusters.to.keep)]
                  
                  
      #Save Outputs:
      
      save(subset_data, file = paste("Data/Basic_Celltype_Cell_IDs/filtered1/", subset.name, "_filtered1_srt.Rdata", sep = ""))
      
      write.table(all.markers, file = paste("Data/Basic_Celltype_Cell_IDs/filtered1/", subset.name, "_filtered1_markers.tsv", sep = ""), row.names=FALSE, sep="\t")

      save(cell.IDs.to.keep, file = paste("Data/Basic_Celltype_Cell_IDs/filtered1/", subset.name, "_filtered1_Cell_IDs.Rdata", sep = ""))
      
    }


#2nd Filter of each Basic Celltype Subset:

second.round.cell.type.file.list <- list.files(path = "Data/Basic_Celltype_Cell_IDs/filtered1/")
second.round.cell.type.file.list <- second.round.cell.type.file.list[which(str_detect(string = second.round.cell.type.file.list, pattern = "_filtered1_Cell_IDs.Rdata"))]


    for (file in second.round.cell.type.file.list) {
      
            subset.name <- str_replace(file, pattern = "_filtered1_Cell_IDs.Rdata", replacement = "")
            
            cell.IDs <- get(load(paste("Data/Basic_Celltype_Cell_IDs/filtered1/", subset.name, "_filtered1_Cell_IDs.Rdata", sep = "")))
            
            #Filter by Cell IDs:
            subset_data <- contact_derm@assays$RNA@counts[,which(colnames(contact_derm@assays$RNA@counts) %in% cell.IDs)]
            
            #Reprocess in Seurat:
                #Create Suerat Object, normalize, and scale:
                subset_data <- CreateSeuratObject(counts = subset_data, project = "Multi_disease_analysis", min.cells = 3, min.features = 150)
                
                #Try SCTransform for processing:
                subset_data <- SCTransform(subset_data, variable.features.n = 4000, ncells = 10000)
                
                #PCA, Clustering, and UMAP:
                subset_data <- RunPCA(subset_data, npcs = 75)
                
                subset_data <- FindNeighbors(subset_data, dims = 1:75)
                
                if (subset.name == "LYMPH") {
                  subset_data <- FindClusters(subset_data, resolution = 2.0)
                }
                
                if (subset.name != "LYMPH") {
                  subset_data <- FindClusters(subset_data, resolution = 1.0)
                }
                
                subset_data <- RunUMAP(subset_data, dims = 1:75)
                
                DimPlot(subset_data, reduction = "umap", label = TRUE, pt.size = 0.75) + NoLegend()
                
                all.markers <- FindAllMarkers(subset_data, min.pct = 0.25, logfc.threshold = 0.25)
                
                
                      # Filter clusters that appear to be doublets:
                      unique.cluster.labels <- as.character(unique(subset_data@meta.data$seurat_clusters))
                      
                      clusters.to.remove <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>%
                        filter(gene %in% Basic_cell_specific_markers[-which(Basic_cell_specific_markers %in% eval(as.symbol(paste(subset.name,"_specific_markers", sep = ""))))]) %>% 
                        dplyr::select(cluster) %>% unique()
                      clusters.to.remove <- as.vector(clusters.to.remove$cluster)
                      
                      if (length(clusters.to.remove) == 0) {
                        clusters.to.keep <- unique.cluster.labels
                      }
                      
                      if (length(clusters.to.remove) != 0) {
                        clusters.to.keep <- unique.cluster.labels[-which(unique.cluster.labels %in% clusters.to.remove)]
                      }
                      
                cell.IDs.to.keep <- colnames(subset_data)[which(subset_data@meta.data$seurat_clusters %in% clusters.to.keep)]
                
            
      #Save Outputs:
      
      save(subset_data, file = paste("Data/Basic_Celltype_Cell_IDs/filtered2/", subset.name, "_filtered2_srt.Rdata", sep = ""))
      
      write.table(all.markers, file = paste("Data/Basic_Celltype_Cell_IDs/filtered2/", subset.name, "_filtered2_markers.tsv", sep = ""), row.names=FALSE, sep="\t")
      
      save(cell.IDs.to.keep, file = paste("Data/Basic_Celltype_Cell_IDs/filtered2/", subset.name, "_filtered2_Cell_IDs.Rdata", sep = ""))
      
    }