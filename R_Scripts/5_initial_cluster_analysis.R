#############################################################
### Initial Cluster Analysis
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
library(ggrepel)
library(data.table)


setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load data:
load(file = "Data/unfiltered.seurat.Rdata")

#Append MetaData:

    #Cell ID minus barcode:
    master_data@meta.data$Cell_ID_lessbarcode <- str_replace_all(rownames(master_data@meta.data), pattern = "_[:upper:]*$", replacement = "")

    #Patient:
    master_data@meta.data$Patient <- as.character(master_data@meta.data$orig.ident)
    
      #run twice to replace 2 levels of padded zeros:
      master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "^[:alpha:]*0")),"Patient"] <- master_data@meta.data$Patient[which(str_detect(string = master_data@meta.data$Patient, pattern = "^[:alpha:]*0"))] %>% str_replace(string = ., pattern = "0", replacement = "")
      master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "^[:alpha:]*0")),"Patient"] <- master_data@meta.data$Patient[which(str_detect(string = master_data@meta.data$Patient, pattern = "^[:alpha:]*0"))] %>% str_replace(string = ., pattern = "0", replacement = "")

    #Disease:
    master_data@meta.data$Disease <- NA  
  
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "CB")),"Disease"] <- "Healthy"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "HC")),"Disease"] <- "Allergy"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "CL")),"Disease"] <- "Lupus"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "DM")),"Disease"] <- "Dermatomyositis"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "VB")),"Disease"] <- "Vitiligo"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Patient, pattern = "VB81|VB92|VB93|VB97|VB98|VB102")),"Disease"] <- "Psoriasis"
    
    master_data@meta.data$Disease <- factor(master_data@meta.data$Disease, 
                                               levels = c("Healthy", "Allergy"))

    #Lesion:
    master_data@meta.data$Lesion <- NA
    
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Cell_ID_lessbarcode, pattern = "_H[:digit:]_|_H_")),"Lesion"] <- "Nonlesional"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Cell_ID_lessbarcode, pattern = "_NL_")),"Lesion"] <- "Nonlesional"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Cell_ID_lessbarcode, pattern = "Day2_SADBE|D2_SADBE")),"Lesion"] <- "Day2_Allergy"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Cell_ID_lessbarcode, pattern = "Day4_SADBE|D4_SADBE")),"Lesion"] <- "Day4_Allergy"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Cell_ID_lessbarcode, pattern = "SLS")),"Lesion"] <- "Irritant"
    master_data@meta.data[which(str_detect(string = master_data@meta.data$Cell_ID_lessbarcode, pattern = "Acetone")),"Lesion"] <- "Acetone_Vehicle"
    
    master_data@meta.data$Lesion <- factor(master_data@meta.data$Lesion, 
                                            levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
    
    #Blister Visit:
    master_data@meta.data$Visit <- NA
    
    #Sequence Date:
    split_sample_to_elements <- function(input) {
      
      name.list.length <- as.numeric(lengths(str_split(input, pattern = "_")))
      
      name <- input
      
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
    cellID_elements <- as.data.frame(matrix(unlist(t(sapply(master_data@meta.data$Cell_ID_lessbarcode, FUN = split_sample_to_elements))), ncol = 3))
    colnames(cellID_elements) <- c("Patient", "Lesion", "Sequence_date")
    
    
    master_data@meta.data$Sequence_date <- cellID_elements$Sequence_date
    master_data@meta.data$Sequence_date <- factor(master_data@meta.data$Sequence_date,
                                                  levels = c("4_16_20", "1_20_21", "12_3_21", "4_28_22"))
   
    #Sequence Instrument:
    master_data@meta.data$Sequence_instrument <- NA
    master_data@meta.data[which(master_data@meta.data$Sequence_date %in% c("4_16_20", "1_20_21")),"Sequence_instrument"] <- "NB501205"
    master_data@meta.data[which(master_data@meta.data$Sequence_date %in% c("12_3_21", "4_28_22")),"Sequence_instrument"] <- "VH00230"
    master_data@meta.data[which(master_data@meta.data$Patient %in% c("CB17", "CB20", "CB23")),"Sequence_instrument"] <- "NS500602"
    master_data@meta.data[which(master_data@meta.data$Patient %in% c("CB31")),"Sequence_instrument"] <- "NB501205"
    master_data@meta.data[which(master_data@meta.data$Patient %in% c("CB27")),"Sequence_instrument"] <- "A00439"
    master_data@meta.data[which(master_data@meta.data$Patient %in% c("CB43", "CB44", "CB46", "CB48")),"Sequence_instrument"] <- "A00197"
    
    master_data@meta.data[which(master_data@meta.data$Cell_ID_lessbarcode %in% c("CB019_H1_V1_S1")),"Sequence_instrument"] <- "NS500602"
    master_data@meta.data[which(master_data@meta.data$Cell_ID_lessbarcode %in% c("CB019_H2_V1_S1")),"Sequence_instrument"] <- "NB501205"
    
    master_data@meta.data[which(master_data@meta.data$Cell_ID_lessbarcode %in% c("CB021_H1_V1_S1")),"Sequence_instrument"] <- "NB501205"
    master_data@meta.data[which(master_data@meta.data$Cell_ID_lessbarcode %in% c("CB021_H2_V1_S1")),"Sequence_instrument"] <- "NS500602"
    
    master_data@meta.data[which(master_data@meta.data$Cell_ID_lessbarcode %in% c("CB022_H1_V1_S1", "CB022_H1_V1_S3", "CB022_H1_V1_S4", "CB022_H1_V1_S5", "CB022_H2_V1_S1", "CB022_H2_V1_S3", "CB022_H2_V1_S4", "CB022_H2_V1_S5")),"Sequence_instrument"] <- "NS500602"
    master_data@meta.data[which(master_data@meta.data$Cell_ID_lessbarcode %in% c("CB022_H1_V1_S2", "CB022_H2_V1_S2")),"Sequence_instrument"] <- "NB501205"
    
    master_data@meta.data$Sequence_instrument <- factor(master_data@meta.data$Sequence_instrument,
                                                      levels = c("NS500602", "NB501205", "VH00230", "A00439", "A00197"))
    
    #Sequence Platform:
    master_data@meta.data$Sequence_platform <- NA
    master_data@meta.data[which(master_data@meta.data$Sequence_instrument %in% c("NS500602", "NB501205")),"Sequence_platform"] <- "NextSeq550"
    master_data@meta.data[which(master_data@meta.data$Sequence_instrument %in% c("VH00230")),"Sequence_platform"] <- "NextSeq2000"
    master_data@meta.data[which(master_data@meta.data$Sequence_instrument %in% c("A00439", "A00197")),"Sequence_platform"] <- "NovaSeq"
    
    master_data@meta.data$Sequence_platform <- factor(master_data@meta.data$Sequence_platform,
                                                      levels = c("NextSeq550", "NextSeq2000", "NovaSeq"))


#Inspect data visually:
    plot_tsne_metadata_srt(master_data, color_by = "seurat_clusters", facet_by = c("Lesion", "Sequence_date"), facet_grid = TRUE, plot_legend = FALSE, plot_labels = FALSE)
    plot_tsne_metadata_srt(master_data, color_by = "seurat_clusters", size = 0.5, facet_by = c("Sequence_platform", "Sequence_instrument"), facet_grid = TRUE, plot_legend = FALSE, plot_labels = FALSE)
    plot_tsne_metadata_srt(master_data, color_by = "seurat_clusters", size = 0.5, facet_by = c("Lesion", "Sequence_instrument"), facet_grid = TRUE, plot_legend = FALSE, plot_labels = FALSE)
    plot_tsne_metadata_srt(master_data, color_by = "seurat_clusters", size = 0.1, facet_by = c("Sequence_instrument", "Patient"), facet_grid = TRUE, plot_legend = FALSE, plot_labels = FALSE)
   
    #setwd("Plots/")
    #ggsave(filename = "unfiltered_UMAP_Sequencer_vs_Patient.png", height = 40, width = 80, scale = 0.25)
   
  #Inspect Clusters:
    DimPlot(master_data, reduction = "umap", label = TRUE, pt.size = 0.75) + NoLegend()
    
    #DimPlot(master_data, reduction = "umap", label = FALSE, pt.size = 0.75, group.by = "Disease")
    
    #FeaturePlot(master_data, features = "IL2RB")
    
    #VlnPlot(master_data, features = "IL22", group.by = "Disease")
    

### Cluster Classification & Filtration of RBCS:


#Load data:
all.cluster.markers <- fread(file = "Data/unfiltered.cluster.markers.tsv", header = TRUE)
all.cluster.markers <- as.data.frame(all.cluster.markers)

top10 <- all.cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


#Create dataframe to categorize each seurat cluster:
Basic_cluster_categories <- data.frame(cluster = unique(all.cluster.markers$cluster),
                                       KRT = rep(NA, times = length(unique(all.cluster.markers$cluster))),
                                       MEL = rep(NA, times = length(unique(all.cluster.markers$cluster))),
                                       APC = rep(NA, times = length(unique(all.cluster.markers$cluster))),
                                       LYMPH = rep(NA, times = length(unique(all.cluster.markers$cluster))),
                                       RBC = rep(NA, times = length(unique(all.cluster.markers$cluster))))
    
    #Analyze top10 markers for basic cell type markers:
    for (i in 1:nrow(Basic_cluster_categories)) {
      
      #each cluster markers:
      clus <-Basic_cluster_categories[i,"cluster"]
      clus.markers <- top10 %>% filter(cluster == clus) %>% ungroup() %>% select(gene) %>% unlist()
      
      #KRT:
      if (any(str_detect(string = clus.markers, pattern = "KRT1|KRT2|KRT5|KRT6A|KRT6B|KRT10|KRT14|KRT15|KRT17|KRT77|KRTDAP|DMKN|MUCL1"))) {
        Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "KRT"] <- 1
      }
      
      #MEL:
      if (any(str_detect(string = clus.markers, pattern = "TYR|MLANA|PMEL|DCT"))) {
        Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "MEL"] <- 1
      }
      
      #APC:
      if (any(str_detect(string = clus.markers, pattern = "HLA.DR|LYZ|CD207|CD86|CD83|CLEC|CD163|MRC1|NRP1|HMOX1|JCHAIN"))) {
        Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "APC"] <- 1
      }
      
      #LYMPH:
      if (any(str_detect(string = clus.markers, pattern = "TRAC|TRDC|NCAM1"))) {
        Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "LYMPH"] <- 1
      }
      
      #RBC:
      if (any(str_detect(string = clus.markers, pattern = "HBB|HBA2"))) {
        Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "RBC"] <- 1
      }
      
    }
    
    Basic_cluster_categories$rowsum <- apply(Basic_cluster_categories[,-1], MARGIN = 1, FUN = function(x) sum(x, na.rm = TRUE))
    
        #Some clusters aren't obvious by top 10 markers, so look at top 50 markers:
        no.category <- Basic_cluster_categories %>% filter(rowsum == 0) %>% ungroup() %>% select(cluster) %>% unlist() %>% as.numeric()
        
        difficult.clus.top50 <- all.cluster.markers %>% filter(cluster %in% no.category) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
        
        for (j in 1:length(no.category)) {
          
          #each cluster markers:
          clus <- no.category[j]
          clus.markers <- difficult.clus.top50 %>% filter(cluster == clus) %>% ungroup() %>% select(gene) %>% unlist()
          
          #KRT:
          if (any(str_detect(string = clus.markers, pattern = "KRT1|KRT2|KRT5|KRT6A|KRT6B|KRT10|KRT14|KRT15|KRT17|KRT77|KRTDAP|DMKN|MUCL1"))) {
            Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "KRT"] <- 1
          }
          
          #MEL:
          if (any(str_detect(string = clus.markers, pattern = "TYR|MLANA|PMEL|DCT|MITF"))) {
            Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "MEL"] <- 1
          }
          
          #APC:
          if (any(str_detect(string = clus.markers, pattern = "HLA.DR|LYZ|CD207|CD86|CD83|CLEC|CD163|MRC1|NRP1"))) {
            Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "APC"] <- 1
          }
          
          #LYMPH:
          if (any(str_detect(string = clus.markers, pattern = "TRAC|TRDC|NCAM1"))) {
            Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "LYMPH"] <- 1
          }
          
          #RBC:
          if (any(str_detect(string = clus.markers, pattern = "HBB|HBA2"))) {
            Basic_cluster_categories[which(Basic_cluster_categories$cluster == clus), "RBC"] <- 1
          }
          
          
        }
        
        Basic_cluster_categories$rowsum <- apply(Basic_cluster_categories[,2:6], MARGIN = 1, FUN = function(x) sum(x, na.rm = TRUE))
        
        Basic_cluster_categories <- Basic_cluster_categories %>% mutate(Basic_Celltype = case_when(rowsum > 1 ~ "doublet",
                                                                                                   rowsum == 1 & KRT == 1 ~ "KRT",
                                                                                                   rowsum == 1 & MEL == 1 ~ "MEL",
                                                                                                   rowsum == 1 & APC == 1 ~ "APC",
                                                                                                   rowsum == 1 & LYMPH == 1 ~ "LYMPH",
                                                                                                   rowsum == 1 & RBC == 1 ~ "RBC",
                                                                                                   rowsum == 0 ~ "unknown"))
        
        #Some clusters still are unclear, so look at markers manually:
        unknowns <- filter(Basic_cluster_categories, Basic_Celltype == "unknown")
        unknown.clus.markers <- all.cluster.markers %>% filter(cluster %in% unknowns$cluster) %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
        
            #This cluster clearly looks like keratinocytes:
            Basic_cluster_categories$Basic_Celltype[which(Basic_cluster_categories$Basic_Celltype == "unknown")] <- "KRT"
            
        Basic_cluster_categories$cluster <- factor(Basic_cluster_categories$cluster, 
                                                   levels = unique(Basic_cluster_categories$cluster))

        master_data@meta.data <- left_join(master_data@meta.data, Basic_cluster_categories[,which(colnames(Basic_cluster_categories) %in% c("cluster", "Basic_Celltype"))], by = c("seurat_clusters" = "cluster"))
        
        plot_tsne_metadata_srt(master_data, color_by = "Basic_Celltype", plot_legend = TRUE)
  
        
    #Filter RBCs, doublets, and any cluster with fewer than 10 cells:
      cells.per.cluster <- master_data@meta.data %>% group_by(seurat_clusters) %>% summarize(n())
      colnames(cells.per.cluster) <- c("cluster", "cell_count")
      Basic_cluster_categories <- full_join(Basic_cluster_categories, cells.per.cluster)
      
      Basic_cluster_categories <- Basic_cluster_categories %>% filter(!Basic_Celltype %in% c("RBC", "doublet"))
       
      
      clusters.to.keep <- as.character(unique(Basic_cluster_categories$cluster))

      contact_derm <- subset(x = master_data, idents = clusters.to.keep)

      contact_derm@meta.data <- master_data@meta.data[which(master_data@meta.data$seurat_clusters %in% clusters.to.keep),]
      contact_derm@meta.data$seurat_clusters <- factor(contact_derm@meta.data$seurat_clusters,
                                                            levels = as.character(sort(as.numeric(unique(contact_derm@meta.data$seurat_clusters)))))
      
      #save(contact_derm, file = "Data/contact_derm_5_9_22.Rdata")

  #Inspect filtered clusters:
      load("Data/contact_derm_5_9_22.Rdata")
      
      plot_tsne_metadata_srt(contact_derm, color_by = "Basic_Celltype", plot_legend = TRUE)
      
 #Save vectors of Cell IDs for subsetted reprocessing by each Basic Celltype:
      #Filter RBCs and unknown weird cells:
      basic.celltypes <- as.character(unique(contact_derm@meta.data$Basic_Celltype))
                                      
      for (celltype in basic.celltypes) {
        
        cell.IDs.to.keep <- colnames(contact_derm)[which(contact_derm@meta.data$Basic_Celltype == celltype)]
        
        assign(paste(celltype,"_Cell_IDs", sep = ""), cell.IDs.to.keep)
      
      }
      
      #Save Cell ID vectors:
          #save(KRT_Cell_IDs, file = "Data/Basic_Celltype_Cell_IDs/initial/KRT_Cell_IDs.Rdata")
          #save(MEL_Cell_IDs, file = "Data/Basic_Celltype_Cell_IDs/initial/MEL_Cell_IDs.Rdata")
          #save(APC_Cell_IDs, file = "Data/Basic_Celltype_Cell_IDs/initial/APC_Cell_IDs.Rdata")
          #save(LYMPH_Cell_IDs, file = "Data/Basic_Celltype_Cell_IDs/initial/LYMPH_Cell_IDs.Rdata")
          
          

  