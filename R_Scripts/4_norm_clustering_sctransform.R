######################
#Pre-processing with Seurat:
######################

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
library(sctransform)


setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load Data
master_data <- fread(file = "Data/master_data_5_6_22.tsv", header = TRUE)

#Filter contact derm, GVHD, and SJS out of this analysis:
if (any(str_detect(string = master_data$Cell_ID, pattern = "GH|SJS")) == TRUE) {
  
  master_data <- master_data[-which(str_detect(string = master_data$Cell_ID, pattern = "GH|SJS")),]
}

#Filter Patients with less than 100 cells:
  split_cellID_to_elements <- function(input) {
    
    name.list.length <- as.numeric(lengths(str_split(input, pattern = "_")))
    
    name <- str_c(unlist(str_split(input, pattern = "_"))[1:(name.list.length-1)], collapse = "_")
    
    #Remove "A / B" sample splits
    if (str_detect(name, pattern = "_A_|_B_")) {
      name <- str_replace(name, pattern = "_A_|_B_", replacement = "_") 
    }
    
    name.list.length <- as.numeric(lengths(str_split(name, pattern = "_")))
    
    Patient <- unlist(str_split(name, pattern = "_"))[1]
    
    #remove padded zeros:
    if (str_detect(Patient, pattern = "^[:alpha:]*0")) {
      Patient <- str_replace(Patient, pattern = "0", replacement = "")
    }
    if (str_detect(Patient, pattern = "^[:alpha:]*0")) {
      Patient <- str_replace(Patient, pattern = "0", replacement = "")
    }
    
    #Process Controls here:
    if (str_detect(Patient, pattern = "CB")) {
      
      Lesion <- "NL"
      
      sequence_date <- NA
    }
    
    #Process Contact Derm here:
    if (str_detect(Patient, pattern = "HC")) {
      
      sequence_date <- str_c(unlist(str_split(name, pattern = "_"))[(name.list.length-2):(name.list.length)], collapse = "_")
      Lesion <- str_c(unlist(str_split(name, pattern = "_"))[2:(name.list.length-3)], collapse = "_")
      
      output <- list(Patient,Lesion,sequence_date)  
      names(output) <- c("Patient", "Lesion", "sequence_date")
    }
    
    output <- list(Patient,Lesion,sequence_date)  
    names(output) <- c("Patient", "Lesion", "sequence_date")
    return(output)
}
  
  cell.info <- data.frame(cellID = master_data$Cell_ID,
                            Cell_ID_lessbarcode = str_replace_all(master_data$Cell_ID, pattern = "_[:upper:]*$", replacement = ""))
  
  
  meta.data <- as.data.frame(matrix(unlist(t(sapply(cell.info$cellID, FUN = split_cellID_to_elements))), ncol = 3))
  colnames(meta.data) <- c("Patient", "Lesion", "Sequence_date")
  
  cell.info <- cbind(cell.info, meta.data)

  patient.info <- cell.info %>% group_by(Patient) %>% summarize(cell.count = n()) %>% arrange(cell.count)
  patients.to.filter <- patient.info %>% filter(cell.count < 100) %>% dplyr::select(Patient) %>% unlist()

  cellIDs.to.filter <- cell.info %>% filter(Patient %in% patients.to.filter) %>% dplyr::select(cellID) %>% unlist()
  
  if (length(cellIDs.to.filter) > 0) {
    master_data <- master_data[-which(master_data$Cell_ID %in% cellIDs.to.filter),]
  }

#Remove Summary Stat Columns:
master_data <- as.data.frame(master_data)

summary.stat.cols <- c("UMI.Sum", "gene.Sum", "mito.UMI.Sum", "percent.mt")
master_data <- master_data[,-which(colnames(master_data) %in% summary.stat.cols)]

#Prepare matrix for Seurat:
master_data <- column_to_rownames(master_data, var = "Cell_ID")
sample.ind <- grep("Sample", colnames(master_data))
master_data <- t(master_data[,-sample.ind])

#Create Suerat Object:
master_data <- CreateSeuratObject(counts = master_data, project = "Contact_derm", min.cells = 3, min.features = 150)


#Try SCTransform for processing:
master_data <- SCTransform(master_data, variable.features.n = 4000, ncells = 10000)

    #PCA, Clustering, and UMAP:
    master_data <- RunPCA(master_data, npcs = 75)
        
    master_data <- FindNeighbors(master_data, dims = 1:75)
    master_data <- FindClusters(master_data, resolution = 0.5)
    master_data <- RunUMAP(master_data, dims = 1:75)
    
    DimPlot(master_data, reduction = "umap", label = TRUE, pt.size = 0.75) + NoLegend()
        
        #top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
        
  #Save output:
    save_path <- "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/"
    
    save(master_data, file = paste(save_path, "unfiltered.seurat.Rdata", sep =""))
    
    
    all.markers <- FindAllMarkers(master_data, min.pct = 0.25, logfc.threshold = 0.25)
    
    write.table(all.markers, file = paste(save_path, "unfiltered.cluster.markers.tsv", sep = ""), row.names=FALSE, sep="\t")



