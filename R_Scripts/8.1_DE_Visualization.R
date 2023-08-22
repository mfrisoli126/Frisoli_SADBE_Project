#######################################################
# DE Analysis: Visualizing Significant Genes
#######################################################
library(Seurat)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(purrr)
library(readxl)
library(data.table)
library(SignallingSingleCell)
library(scales)
library(clusterProfiler)
library(aliases2entrez)
library(writexl)

#Load data:

setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/")

#Load data:
load(file = "Data/contact_derm_filtered.Rdata")

contact_derm@meta.data$Semi_Detailed_Celltype <- factor(contact_derm@meta.data$Semi_Detailed_Celltype, 
                                                        levels = c("KRT", "MEL", "LC","Mac1", "Mac2", "Bcell", "CD4", "Treg", "CD8", "NK"))
contact_derm@meta.data$Detailed_Celltype <- factor(contact_derm@meta.data$Detailed_Celltype, 
                                                   levels = c("KRT-b1", "KRT-b2", "KRT-b3", "KRT-sp", "KRT-wr1", "KRT-wr2", "KRT-m", "KRT-g", "MEL1", "MEL2", "LC", "Mac1", "Mac2", "Bcell", "CD4", "Treg", "CD8", "NK"))

contact_derm@meta.data$Lesion <- factor(contact_derm@meta.data$Lesion, 
                                                        levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))

      
      
      ####################################
      # Semi_Detailed_CellType DE by EdgeR
      ####################################
      
      
      # Set Path to group of DE files:
      
      file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_DE/"
      
      
      filelist = list.files(path = file.path, pattern = ".txt")
      
      DE.data.list <- c()
      
      for (i in 1:length(filelist)) {
        
        #Load Data:
        DE.data <- fread(paste(file.path, filelist[i], sep = ""))
        DE.data <- DE.data %>% dplyr::rename("Gene" = V1)
        
        #Note: We used EdgeR package to compute DE, which labels the Benjamini-Hochberg adjusted Pvalue as "FDR". DESeq2 labels the same column "padj" as it is the adjusted p value.
        DE.data <- DE.data %>% dplyr::rename("padj" = FDR)
        
        #Gather info from file name:
        string <- filelist[i]
        string <- str_replace(string, pattern = ".txt", replacement = "")
        string <- unlist(str_split(string, pattern = "_vs_"))  
        string <- unlist(str_split(string, pattern = "_", n = 2))  
        
        CellType <- string[1]
        Comparison <- string[2]
        Reference <- string[3]
        
        #Add labels:
        DE.data$CellType <- CellType
        DE.data$Comparison <- Comparison
        DE.data$Reference <- Reference
        
        
        assign(paste(CellType,"_",Comparison,".data", sep = ""), DE.data)
        
        DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,".data", sep = ""))
        
      }
      
      DE.data <- purrr::reduce(mget(DE.data.list), full_join)
      
      #Add columns for plotting:
      
      #Log-modulus transformation from limma package:
      logfc2fc <- function(logFC){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(logFC>=0))
        fc <- sgn*2^abs(logFC)
        return(fc)}
      
      
      fc_backto_logFC <- function(foldchange){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(foldchange>=0))
        logFC = sign(foldchange)*(log2(abs(foldchange)))
        return(logFC)}
      
      DE.data <- DE.data %>% mutate(FoldChange = logfc2fc(logFC))
      
      DE.data <- DE.data %>% mutate(FoldChange_Direction = ifelse(FoldChange > 1, "Up", "Down"))
      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & FoldChange >= 2^0.5 ~ "Significantly Up", 
                                                             FoldChange_Direction == "Down" & padj <= 0.01 & FoldChange <= 2^0.5 ~ "Significantly Down",
                                                             padj > 0.01 ~ "Not Significant",
                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
      
      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
      
      DE.data$Comparison <- factor(DE.data$Comparison, 
                                   levels = c("Irritant", "Day2_Allergy", "Day4_Allergy"))
      
      DE.data$CellType <- factor(DE.data$CellType, 
                                 levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))
                                 
      #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt", sep = ""), sep = "\t")
      
      
      
      ############################################################################################################
      # Semi_Detailed_CellType DE by EdgeR: filtered to patch test patients with batch covariate by EdgeR
      ############################################################################################################
      
      # Set Path to group of DE files:
      
      file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_patch_patient_batched_DE/"
      
      
      filelist = list.files(path = file.path, pattern = ".txt")
      
      DE.data.list <- c()
      
      for (i in 1:length(filelist)) {
        
        #Load Data:
        DE.data <- fread(paste(file.path, filelist[i], sep = ""))
        DE.data <- DE.data %>% dplyr::rename("Gene" = V1)
        
        #Note: We used EdgeR package to compute DE, which labels the Benjamini-Hochberg adjusted Pvalue as "FDR". DESeq2 labels the same column "padj" as it is the adjusted p value.
        DE.data <- DE.data %>% dplyr::rename("padj" = FDR)
        
        #Gather info from file name:
        string <- filelist[i]
        string <- str_replace(string, pattern = ".txt", replacement = "")
        string <- unlist(str_split(string, pattern = "_vs_"))  
        string <- unlist(str_split(string, pattern = "_", n = 2))  
        
        CellType <- string[1]
        Comparison <- string[2]
        Reference <- string[3]
        
        #Add labels:
        DE.data$CellType <- CellType
        DE.data$Comparison <- Comparison
        DE.data$Reference <- Reference
        
        
        assign(paste(CellType,"_",Comparison,".data", sep = ""), DE.data)
        
        DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,".data", sep = ""))
        
      }
      
      DE.data <- purrr::reduce(mget(DE.data.list), full_join)
      
      #Add columns for plotting:
      
      #Log-modulus transformation from limma package:
      logfc2fc <- function(logFC){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(logFC>=0))
        fc <- sgn*2^abs(logFC)
        return(fc)}
      
      
      fc_backto_logFC <- function(foldchange){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(foldchange>=0))
        logFC = sign(foldchange)*(log2(abs(foldchange)))
        return(logFC)}
      
      DE.data <- DE.data %>% mutate(FoldChange = logfc2fc(logFC))
      
      DE.data <- DE.data %>% mutate(FoldChange_Direction = ifelse(FoldChange > 1, "Up", "Down"))
      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & FoldChange >= 2^0.5 ~ "Significantly Up", 
                                                             FoldChange_Direction == "Down" & padj <= 0.01 & FoldChange <= 2^0.5 ~ "Significantly Down",
                                                             padj > 0.01 ~ "Not Significant",
                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
      
      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
      
      DE.data$Comparison <- factor(DE.data$Comparison, 
                                   levels = c("Irritant", "Day2_Allergy", "Day4_Allergy"))
      
      DE.data$CellType <- factor(DE.data$CellType, 
                                 levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))
      
      #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt", sep = ""), sep = "\t")
      
      ############################################################################################################
      # Semi_Detailed_CellType DE by EdgeR: filtered to patch test patients (without any batch covariate terms) by EdgeR
      ############################################################################################################
      
      # Set Path to group of DE files:
      
      file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_patch_patient_DE/"
      
      
      filelist = list.files(path = file.path, pattern = ".txt")
      
      DE.data.list <- c()
      
      for (i in 1:length(filelist)) {
        
        #Load Data:
        DE.data <- fread(paste(file.path, filelist[i], sep = ""))
        DE.data <- DE.data %>% dplyr::rename("Gene" = V1)
        
        #Note: We used EdgeR package to compute DE, which labels the Benjamini-Hochberg adjusted Pvalue as "FDR". DESeq2 labels the same column "padj" as it is the adjusted p value.
        DE.data <- DE.data %>% dplyr::rename("padj" = FDR)
        
        #Gather info from file name:
        string <- filelist[i]
        string <- str_replace(string, pattern = ".txt", replacement = "")
        string <- unlist(str_split(string, pattern = "_vs_"))  
        string <- unlist(str_split(string, pattern = "_", n = 2))  
        
        CellType <- string[1]
        Comparison <- string[2]
        Reference <- string[3]
        
        #Add labels:
        DE.data$CellType <- CellType
        DE.data$Comparison <- Comparison
        DE.data$Reference <- Reference
        
        
        assign(paste(CellType,"_",Comparison,".data", sep = ""), DE.data)
        
        DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,".data", sep = ""))
        
      }
      
      DE.data <- purrr::reduce(mget(DE.data.list), full_join)
      
      #Add columns for plotting:
      
      #Log-modulus transformation from limma package:
      logfc2fc <- function(logFC){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(logFC>=0))
        fc <- sgn*2^abs(logFC)
        return(fc)}
      
      
      fc_backto_logFC <- function(foldchange){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(foldchange>=0))
        logFC = sign(foldchange)*(log2(abs(foldchange)))
        return(logFC)}
      
      DE.data <- DE.data %>% mutate(FoldChange = logfc2fc(logFC))
      
      DE.data <- DE.data %>% mutate(FoldChange_Direction = ifelse(FoldChange > 1, "Up", "Down"))
      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & FoldChange >= 2^0.5 ~ "Significantly Up", 
                                                             FoldChange_Direction == "Down" & padj <= 0.01 & FoldChange <= 2^0.5 ~ "Significantly Down",
                                                             padj > 0.01 ~ "Not Significant",
                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
      
      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
      
      DE.data$Comparison <- factor(DE.data$Comparison, 
                                   levels = c("Irritant", "Day2_Allergy", "Day4_Allergy"))
      
      DE.data$CellType <- factor(DE.data$CellType, 
                                 levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))
      
      #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt", sep = ""), sep = "\t")
      
      
      #########################################################################################
      # Pseudobulk DE analysis by EdgeR: Aggergated by Semi_Detailed_CellType, Patient, Lesion
      #########################################################################################
      
      
      # Set Path to group of DE files:
      
      file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_Bulk_DE/"
      
      
      filelist = list.files(path = file.path, pattern = ".txt")
      
      DE.data.list <- c()
      
      for (i in 1:length(filelist)) {
        
        #Load Data:
        DE.data <- fread(paste(file.path, filelist[i], sep = ""))
        DE.data <- DE.data[,-1]
      
        #Note: We used EdgeR package to compute DE, which labels the Benjamini-Hochberg adjusted Pvalue as "FDR". DESeq2 labels the same column "padj" as it is the adjusted p value.
        DE.data <- DE.data %>% dplyr::rename("padj" = FDR)
        
        #Gather info from file name:
        string <- filelist[i]
        string <- str_replace(string, pattern = ".txt", replacement = "")
        string <- unlist(str_split(string, pattern = "_vs_"))  
        string <- unlist(str_split(string, pattern = "_", n = 2))  
        
        CellType <- string[1]
        Comparison <- string[2]
        Reference <- string[3]
        
        #Add labels:
        DE.data$CellType <- CellType
        DE.data$Comparison <- Comparison
        DE.data$Reference <- Reference
        
        
        assign(paste(CellType,"_",Comparison,".data", sep = ""), DE.data)
        
        DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,".data", sep = ""))
        
      }
      
      DE.data <- purrr::reduce(mget(DE.data.list), full_join)
      
      #Add columns for plotting:
      
      #Log-modulus transformation from limma package:
      logfc2fc <- function(logFC){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(logFC>=0))
        fc <- sgn*2^abs(logFC)
        return(fc)}
      
      
      fc_backto_logFC <- function(foldchange){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(foldchange>=0))
        logFC = sign(foldchange)*(log2(abs(foldchange)))
        return(logFC)}
      
      DE.data <- DE.data %>% mutate(FoldChange = logfc2fc(logFC))
      
      DE.data <- DE.data %>% mutate(FoldChange_Direction = ifelse(FoldChange > 1, "Up", "Down"))
      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & FoldChange >= 2^0.5 ~ "Significantly Up", 
                                                             FoldChange_Direction == "Down" & padj <= 0.01 & FoldChange <= 2^0.5 ~ "Significantly Down",
                                                             padj > 0.01 ~ "Not Significant",
                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
      
      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
      
      DE.data$Comparison <- factor(DE.data$Comparison, 
                                   levels = c("Irritant", "Day2_Allergy", "Day4_Allergy"))
      
      DE.data$CellType <- factor(DE.data$CellType, 
                                 levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))
      
      #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt", sep = ""), sep = "\t")

      ############################################################################################################
      # Semi_Detailed_CellType DE by EdgeR: Batch-corrected with Combat-seq before DE
      ############################################################################################################
      
      # Set Path to group of DE files:
      
      file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_ComBatseq_DE/"
      
      
      filelist = list.files(path = file.path, pattern = ".txt")
      
      DE.data.list <- c()
      
      for (i in 1:length(filelist)) {
        
        #Load Data:
        DE.data <- fread(paste(file.path, filelist[i], sep = ""))
        DE.data <- DE.data %>% dplyr::rename("Gene" = V1)
        
        #Note: We used EdgeR package to compute DE, which labels the Benjamini-Hochberg adjusted Pvalue as "FDR". DESeq2 labels the same column "padj" as it is the adjusted p value.
        DE.data <- DE.data %>% dplyr::rename("padj" = FDR)
        
        #Gather info from file name:
        string <- filelist[i]
        string <- str_replace(string, pattern = ".txt", replacement = "")
        string <- unlist(str_split(string, pattern = "_vs_"))  
        string <- unlist(str_split(string, pattern = "_", n = 2))  
        
        CellType <- string[1]
        Comparison <- string[2]
        Reference <- string[3]
        
        #Add labels:
        DE.data$CellType <- CellType
        DE.data$Comparison <- Comparison
        DE.data$Reference <- Reference
        
        
        assign(paste(CellType,"_",Comparison,".data", sep = ""), DE.data)
        
        DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,".data", sep = ""))
        
      }
      
      DE.data <- purrr::reduce(mget(DE.data.list), full_join)
      
      #Add columns for plotting:
      
      #Log-modulus transformation from limma package:
      logfc2fc <- function(logFC){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(logFC>=0))
        fc <- sgn*2^abs(logFC)
        return(fc)}
      
      
      fc_backto_logFC <- function(foldchange){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(foldchange>=0))
        logFC = sign(foldchange)*(log2(abs(foldchange)))
        return(logFC)}
      
      DE.data <- DE.data %>% mutate(FoldChange = logfc2fc(logFC))
      
      DE.data <- DE.data %>% mutate(FoldChange_Direction = ifelse(FoldChange > 1, "Up", "Down"))
      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & FoldChange >= 2^0.5 ~ "Significantly Up", 
                                                             FoldChange_Direction == "Down" & padj <= 0.01 & FoldChange <= 2^0.5 ~ "Significantly Down",
                                                             padj > 0.01 ~ "Not Significant",
                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
      
      DE.data$Comparison <- factor(DE.data$Comparison, 
                                   levels = c("Irritant", "Day2_Allergy", "Day4_Allergy"))
      
      DE.data$CellType <- factor(DE.data$CellType, 
                                 levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))
      
      #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt", sep = ""), sep = "\t")

      
      ##########################################################
      # Semi_Detailed_CellType DE by EdgeR (no min cell filter)
      ##########################################################
      
      
      # Set Path to group of DE files:
      
      file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE_mincell0/Semi_Detailed_Celltype_DE/"
      
      
      filelist = list.files(path = file.path, pattern = ".txt")
      
      DE.data.list <- c()
      
      for (i in 1:length(filelist)) {
        
        #Load Data:
        DE.data <- fread(paste(file.path, filelist[i], sep = ""))
        DE.data <- DE.data %>% dplyr::rename("Gene" = V1)
        
        #Note: We used EdgeR package to compute DE, which labels the Benjamini-Hochberg adjusted Pvalue as "FDR". DESeq2 labels the same column "padj" as it is the adjusted p value.
        DE.data <- DE.data %>% dplyr::rename("padj" = FDR)
        
        #Gather info from file name:
        string <- filelist[i]
        string <- str_replace(string, pattern = ".txt", replacement = "")
        string <- unlist(str_split(string, pattern = "_vs_"))  
        string <- unlist(str_split(string, pattern = "_", n = 2))  
        
        CellType <- string[1]
        Comparison <- string[2]
        Reference <- string[3]
        
        #Add labels:
        DE.data$CellType <- CellType
        DE.data$Comparison <- Comparison
        DE.data$Reference <- Reference
        
        
        assign(paste(CellType,"_",Comparison,".data", sep = ""), DE.data)
        
        DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,".data", sep = ""))
        
      }
      
      DE.data <- purrr::reduce(mget(DE.data.list), full_join)
      
      #Add columns for plotting:
      
      #Log-modulus transformation from limma package:
      logfc2fc <- function(logFC){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(logFC>=0))
        fc <- sgn*2^abs(logFC)
        return(fc)}
      
      
      fc_backto_logFC <- function(foldchange){
        # sign is -1 if logFC<0; 1 if logFC>=0
        sgn <- (-1)^(1+as.numeric(foldchange>=0))
        logFC = sign(foldchange)*(log2(abs(foldchange)))
        return(logFC)}
      
      DE.data <- DE.data %>% mutate(FoldChange = logfc2fc(logFC))
      
      DE.data <- DE.data %>% mutate(FoldChange_Direction = ifelse(FoldChange > 1, "Up", "Down"))
      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & FoldChange >= 2^0.5 ~ "Significantly Up", 
                                                             FoldChange_Direction == "Down" & padj <= 0.01 & FoldChange <= 2^0.5 ~ "Significantly Down",
                                                             padj > 0.01 ~ "Not Significant",
                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
      
      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
      
      DE.data$Comparison <- factor(DE.data$Comparison, 
                                   levels = c("Irritant", "Day2_Allergy", "Day4_Allergy"))
      
      DE.data$CellType <- factor(DE.data$CellType, 
                                 levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))
      
      #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt", sep = ""), sep = "\t")
      
      
##################
# Comparisons:
##################
setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/")

#Open all DE method outputs:
    sc_DE <- fread("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt")
    sc_DE <- sc_DE[,-1]
    
    sc_nofilter_DE <- fread("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE_mincell0/Semi_Detailed_Celltype_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt")
    sc_nofilter_DE <- sc_nofilter_DE[,-1]
    
    bulk_DE <- fread("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_Bulk_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt")
    bulk_DE <- bulk_DE[,-1]
    colnames(bulk_DE)[which(colnames(bulk_DE) == "pvalue")] <- "PValue"
    
    sc_patch_patient_batched_DE <- fread("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_patch_patient_batched_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt")
    sc_patch_patient_batched_DE <- sc_patch_patient_batched_DE[,-1]
    
    sc_patch_patient_DE <- fread("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_patch_patient_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt")
    sc_patch_patient_DE <- sc_patch_patient_DE[,-1]
    
    Combatseq_sc_DE <- fread("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_ComBatseq_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt")
    Combatseq_sc_DE <- Combatseq_sc_DE[,-1]

    
    #Filter all:
      sc_DE_f <- DE_filter_srt(DE_data = sc_DE, single_cell_data = contact_derm, sc_metadata_cols = c("Lesion", "Semi_Detailed_Celltype"))
      sc_nofilter_DE_f <- DE_filter_srt(DE_data = sc_nofilter_DE, single_cell_data = contact_derm, sc_metadata_cols = c("Lesion", "Semi_Detailed_Celltype"))
      bulk_DE_f <- DE_filter_srt(DE_data = bulk_DE, single_cell_data = contact_derm, sc_metadata_cols = c("Lesion", "Semi_Detailed_Celltype"))
      sc_patch_patient_batched_DE_f <- DE_filter_srt(DE_data = sc_patch_patient_batched_DE, single_cell_data = contact_derm, sc_metadata_cols = c("Lesion", "Semi_Detailed_Celltype"))
      sc_patch_patient_DE_f<- DE_filter_srt(DE_data = sc_patch_patient_DE, single_cell_data = contact_derm, sc_metadata_cols = c("Lesion", "Semi_Detailed_Celltype"))
      Combatseq_sc_DE_f <- DE_filter_srt(DE_data = Combatseq_sc_DE, single_cell_data = contact_derm, sc_metadata_cols = c("Lesion", "Semi_Detailed_Celltype"))
      
      #Add method column:
      sc_DE$method <- "default_sc"
      sc_nofilter_DE$method <- "nofilter_default_sc"
      bulk_DE$method <- "bulk"
      sc_patch_patient_batched_DE$method <- "patch_patients_only:batch_covariate"
      sc_patch_patient_DE$method <- "patch_patients_only:no_batch_covariate"
      Combatseq_sc_DE$method <- "combat_seq_sc"
      
      
      sc_DE_f$method <- "default_sc"
      sc_nofilter_DE_f$method <- "nofilter_default_sc"
      bulk_DE_f$method <- "bulk"
      sc_patch_patient_batched_DE_f$method <- "patch_patients_only:batch_covariate"
      sc_patch_patient_DE_f$method <- "patch_patients_only:no_batch_covariate"
      Combatseq_sc_DE_f$method <- "combat_seq_sc"
      
      file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE_mincell0/Semi_Detailed_Celltype_DE/"
      
      #write.table(sc_DE_f, "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged_filtered.txt", sep = "\t")
      #write.table(sc_nofilter_DE_f, "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE_mincell0/Semi_Detailed_Celltype_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged_filtered.txt", sep = "\t")

    
    
    
    
DE.method.comparison <- bind_rows(sc_DE, bulk_DE, sc_patch_patient_batched_DE, sc_patch_patient_DE, Combatseq_sc_DE, sc_nofilter_DE)
DE.method.comparison$Comparison <- str_replace_all(DE.method.comparison$Comparison, pattern = "_", replacement = " ")

DE.method.comparison$Comparison <- factor(DE.method.comparison$Comparison,
                             levels = c("Irritant" , "Day2 Allergy", "Day4 Allergy"))
DE.method.comparison$Color_factor <- factor(DE.method.comparison$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))


#Volcano Plot:
          #For x axis transformation:
          logfc2fc <- function(logFC){
            # sign is -1 if logFC<0; 1 if logFC>=0
            sgn <- (-1)^(1+as.numeric(logFC>=0))
            fc <- sgn*2^abs(logFC)
            return(fc)}
          
          
          fc_backto_logFC <- function(foldchange){
            # sign is -1 if logFC<0; 1 if logFC>=0
            sgn <- (-1)^(1+as.numeric(foldchange>=0))
            logFC = sign(foldchange)*(log2(abs(foldchange)))
            return(logFC)}
          
          tn <- trans_new("log-mod",
                          transform = function(x) fc_backto_logFC(x),
                          inverse = function(x) logfc2fc(x),
                          breaks=c(-10,-5, -2, 1, 2,5, 10,20),
                          domain = c(-Inf, Inf))
          
         
          
          p1 <- DE.method.comparison %>% filter(method == "default_sc") %>% ggplot(., aes(x = FoldChange, y = -log2(padj), color = Color_factor)) + 
                                facet_grid(CellType ~ Comparison, scales = "free") +
                                geom_point(alpha = 0.6) + 
                                theme_classic() + 
                                labs(title = "Single Cell DE (all patient data)", color = "Gene Regulation") + 
                                theme(plot.title = element_text(hjust=0.5)) +
                                scale_color_manual(values = c("coral1", "cornflowerblue", "gray")) +
                                scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn)
          
          p2 <- DE.method.comparison %>% filter(method == "bulk") %>% ggplot(., aes(x = FoldChange, y = -log2(padj), color = Color_factor)) + 
                                facet_grid(CellType ~ Comparison, scales = "free") +
                                geom_point(alpha = 0.6) + 
                                theme_classic() + 
                                labs(title = "Psuedobulk DE (all patient data, filtered at bulks >=5 cells)", color = "Gene Regulation") + 
                                theme(plot.title = element_text(hjust=0.5)) +
                                scale_color_manual(values = c("coral1", "cornflowerblue", "gray")) +
                                scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn)
                              
          
          p3 <- DE.method.comparison %>% filter(method == "patch_patients_only:no_batch_covariate") %>% ggplot(., aes(x = FoldChange, y = -log2(padj), color = Color_factor)) + 
                                facet_grid(CellType ~ Comparison, scales = "free") +
                                geom_point(alpha = 0.6) + 
                                theme_classic() + 
                                labs(title = "Single Cell DE (Patch Test Patients Only; no batch covariate)", color = "Gene Regulation") + 
                                theme(plot.title = element_text(hjust=0.5)) +
                                scale_color_manual(values = c("coral1", "cornflowerblue", "gray")) +
                                scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn)
          
          
          p4 <- DE.method.comparison %>% filter(method == "patch_patients_only:batch_covariate") %>% ggplot(., aes(x = FoldChange, y = -log2(padj), color = Color_factor)) + 
                                facet_grid(CellType ~ Comparison, scales = "free") +
                                geom_point(alpha = 0.6) + 
                                theme_classic() + 
                                labs(title = "Single Cell DE (Patch Test Patients Only; covariate = patient)", color = "Gene Regulation") + 
                                theme(plot.title = element_text(hjust=0.5)) +
                                scale_color_manual(values = c("coral1", "cornflowerblue", "gray")) +
                                scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn)
          
          
          p5 <- DE.method.comparison %>% filter(method == "combat_seq_sc") %>% ggplot(., aes(x = FoldChange, y = -log2(padj), color = Color_factor)) + 
                                facet_grid(CellType ~ Comparison, scales = "free") +
                                geom_point(alpha = 0.6) + 
                                theme_classic() + 
                                labs(title = "Single Cell DE (all data pre-processed by combat-seq)", color = "Gene Regulation") + 
                                theme(plot.title = element_text(hjust=0.5)) +
                                scale_color_manual(values = c("coral1", "cornflowerblue", "gray")) +
                                scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn)
          
          p6 <- DE.method.comparison %>% filter(method == "nofilter_default_sc") %>% ggplot(., aes(x = FoldChange, y = -log2(padj), color = Color_factor)) + 
                                facet_grid(CellType ~ Comparison, scales = "free") +
                                geom_point(alpha = 0.6) + 
                                theme_classic() + 
                                labs(title = "Single Cell DE (all patient data no filter)", color = "Gene Regulation") + 
                                theme(plot.title = element_text(hjust=0.5)) +
                                scale_color_manual(values = c("coral1", "cornflowerblue", "gray")) +
                                scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn)
                                                  
          plot(p1 + p6)
          plot(p1 + p3+ p4 + p5)
          
          
          
          
          DE.method.comparison2 <- bind_rows(sc_DE_f, sc_nofilter_DE_f)
          
          
          extra.genes <- unique(sc_nofilter_DE_f$Gene)[-which(unique(sc_nofilter_DE_f$Gene) %in% unique(sc_DE_f$Gene))]
          extra.genes[which(str_detect(extra.genes, pattern = "^CXCL|^CCL|^IL[:digit:]"))]
          
          plot_violin_srt(contact_derm, gene = "IL13RA1", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype")
          
          rownames(contact_derm)[which(str_detect(rownames(contact_derm), pattern = "IFN"))]
          
  #LogFC DotPlots:
          DE.data.wider <- DE.method.comparison %>% pivot_wider(names_from = method, values_from = c(FoldChange, FoldChange_Direction, Color_factor, logFC, logCPM, LR, PValue, padj))
          DE.data.wider <- DE.data.wider %>% arrange()
          
          
          
          ggplot(DE.data.wider, aes(x = FoldChange_sc, y = `FoldChange_sc patient batched; patch patients only`)) + 
            facet_grid(CellType ~ Comparison, scales = "free") +
            geom_point(alpha = 0.6) + 
            theme_classic() + 
            labs(title = "Gene LogFC comparison using patient batch DE covariate vs. non-batched all data") + 
            theme(plot.title = element_text(hjust=0.5)) +
            scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn) +
            scale_y_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn) +
            geom_vline(aes(xintercept = 2^0.5), color = "#BD3535", linetype = "dashed") +
            geom_vline(aes(xintercept = -2^0.5), color = "#BD3535", linetype = "dashed") +
            geom_hline(aes(yintercept = 2^0.5), color = "#BD3535", linetype = "dashed") +
            geom_hline(aes(yintercept = -2^0.5), color = "#BD3535", linetype = "dashed") 
            
          
          DE.data <- sc_nofilter_DE
          DE.data$Comparison <- str_replace_all(DE.data$Comparison, pattern = "_", replacement = " ")
          
          DE.data$Comparison <- factor(DE.data$Comparison,
                                                    levels = c("Irritant" , "Day2 Allergy", "Day4 Allergy"))
          DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
          
          
          
          ggplot(DE.data, aes(x = FoldChange, y = -log2(padj), color = Color_factor)) + 
            facet_grid(CellType ~ Comparison, scales = "free") +
            geom_point(alpha = 0.6) + 
            theme_classic() + 
            labs(title = "Single Cell DE (all patient data)", color = "Gene Regulation") + 
            theme(plot.title = element_text(hjust=0.5)) +
            scale_color_manual(values = c("coral1", "cornflowerblue", "gray")) +
            scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn)
          
            
          
          
#Heatmap:
contact_derm <- calc_agg_bulk_srt(contact_derm, aggregate_by = c("Lesion", "Semi_Detailed_Celltype"))

#Set DE data source:
sc_nofilter_DE_f <- fread("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE_mincell0/Semi_Detailed_Celltype_DE/Merged_DE_File/Semi_Detailed_CellType_DE.merged_filtered.txt")
sc_nofilter_DE_f <- sc_nofilter_DE_f[,-1]

DE.data <- sc_nofilter_DE_f

#Filter first by padj and fold change:
genes.for.heatmap.pval.fc.filtered <- DE.data %>% filter(padj < 0.01 & abs(logFC > 0.5)) %>% select(Gene) %>% unlist() %>% as.character() %>% unique()

#Next filter Mitochondrial and Ribosomal genes:
if (any(str_detect(genes.for.heatmap.pval.fc.filtered, pattern = "^RPS|^RPL|^MT-")) == TRUE) {
  genes.for.heatmap.pval.fc.filtered <- genes.for.heatmap.pval.fc.filtered[-which(str_detect(genes.for.heatmap.pval.fc.filtered, pattern = "^RPS|^RPL|^MT-"))]
}

  #Next filter by Cell counts and faction positive:
          raw_counts <- t(contact_derm@assays$RNA@counts)
          gene.info.df <- cbind(contact_derm@meta.data[,which(colnames(contact_derm@meta.data) %in% c("Lesion", "Semi_Detailed_Celltype"))], raw_counts)
          gene.info_long <- pivot_longer(gene.info.df, cols = colnames(gene.info.df)[-which(colnames(gene.info.df) %in% c("Lesion", "Semi_Detailed_Celltype"))], names_to = "Gene", values_to = "raw_counts")
          gene.info_long$gene_expression <- gene.info_long$raw_counts
          gene.info_long$gene_expression[gene.info_long$gene_expression >0] <- 1
          
          gene.info_long_filtered <- gene.info_long %>% filter(Gene %in% genes.for.heatmap.pval.fc.filtered)
                
          Cell_count.df <- gene.info.df %>% group_by(Lesion, Semi_Detailed_Celltype) %>% summarize(cell_count = n())
        
          gene.stats.df <- gene.info_long_filtered %>% group_by(Lesion, Semi_Detailed_Celltype, Gene) %>% summarize(cells.expressing = sum(gene_expression),
                                                                                                           mean.expression = mean(raw_counts))
          
          gene.info.df <- full_join(Cell_count.df,gene.stats.df)
          gene.info.df$cells.frac.expressing <- gene.info.df$cells.expressing/gene.info.df$cell_count
          gene.info.df$Lesion <- as.character(gene.info.df$Lesion)
          colnames(gene.info.df)[which(colnames(gene.info.df) == "Semi_Detailed_Celltype")] <- "CellType"
          
          gene.info.df <- DE.data %>% filter(Gene %in% genes.for.heatmap.pval.fc.filtered) %>% left_join(gene.info.df)
          
          gene.info.df$Lesion <- str_replace_all(gene.info.df$Lesion, pattern = "_", replacement = " ")
        
          gene.info.df <- gene.info.df %>% mutate(inclusion.factor = case_when(Lesion == Comparison ~ "include",
                                                                       Lesion == Reference ~ "include",
                                                                       Lesion != Comparison & Lesion != Reference ~ "no"))
          
          gene.info.df <- gene.info.df %>% filter(inclusion.factor == "include")
          
          # Max expected cell number estimation: 
          gene.info.df$max_cell_frac <- NA
          gene.info.df$max_cell_frac <- as.numeric(gene.info.df$max_cell_frac)
          
          for (g in 1:length(unique(gene.info.df$Gene))) {
             
            gene <- unique(gene.info.df$Gene)[g]
            
            gene.data <- gene.info.df %>% filter(Gene == gene)
            
                for (c in 1:length(unique(gene.data$CellType))) {
                  
                 cell <- unique(gene.data$CellType)[c]
                  
                 gene.cell.data <- gene.data %>% filter(CellType == cell)
                  
                 max.frac <- max(gene.cell.data$cells.frac.expressing) 
                 
                 gene.data[which(gene.data$CellType == cell), "max_cell_frac"] <- max.frac
                }
            
            gene.info.df$max_cell_frac[which(gene.info.df$Gene == gene)] <- gene.data$max_cell_frac
            
          }
          
          gene.info.df <- gene.info.df %>% mutate(max_expected_pos_cells = round(cell_count*max_cell_frac, digits = 1))
              
    #Filter by cell stats here:

          gene.filter.df <- data.frame(Gene = unique(gene.info.df$Gene),
                                       max_min_max_expected_pos_cells = NA)
          
          for (g in 1:length(unique(gene.info.df$Gene))) {
            
            gene <- unique(gene.info.df$Gene)[g]
            
            gene.data <- gene.info.df %>% filter(Gene == gene)
            
            gene.comparison.min.expected.cells <- gene.data %>% group_by(CellType, Comparison, Reference) %>% summarise(min_max_expected_pos_cells = min(max_expected_pos_cells))
            
            max_min_max_expected_pos_cells = as.numeric(max(gene.comparison.min.expected.cells$min_max_expected_pos_cells))
          
            gene.filter.df[which(gene.filter.df$Gene == gene), "max_min_max_expected_pos_cells"] <- max_min_max_expected_pos_cells
              
          }
          
          plot(density(gene.filter.df$max_min_max_expected_pos_cells))
          
          
          #max_min_max_expected_pos_cells > 20 looks like a good threshold:
                    #plot_violin_srt(contact_derm, gene = "CBR3", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype")
                    #plot_violin_srt(contact_derm, gene = "PIK3R5", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype")
                    #plot_violin_srt(contact_derm, gene = "KLRG1", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype")
          
          expected_pos_cell_threshold_cutoff = 20
          
          genes.to.filter <- gene.filter.df %>% filter(max_min_max_expected_pos_cells < expected_pos_cell_threshold_cutoff) %>% select(Gene) %>% unlist() %>% as.character()
          
          
    # Finalize genes for DE heatmap:
    if (length(genes.to.filter) >0 ) {
      DE.heatmap.genes <- genes.for.heatmap.pval.fc.filtered[-which(genes.for.heatmap.pval.fc.filtered %in% genes.to.filter)]
    }
    
    if (length(genes.to.filter) == 0 ) {
      DE.heatmap.genes <- genes.for.heatmap.pval.fc.filtered
    }

          
#Scaled by Cell Type:
        all.filtered.DE.genes <- sc_nofilter_DE_f %>% filter(Color_factor != "Not Significant") %>% dplyr::select(Gene) %>% unlist() %>% as.character() %>% unique()

plot_heatmap_srt(contact_derm,
                 genes = all.filtered.DE.genes,
                 type = "bulk",
                 facet_by = "Semi_Detailed_Celltype",
                 scale_group = "Semi_Detailed_Celltype",
                 scale_by = "row",
                 text_angle = 90,
                 text_sizes = c(20, 10, 7, 10, 5, 5, 5),
                 title = "DE genes",
                 color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                 log_scale_base = NULL,
                 ceiling = F,
                 floor = F,
                 gene_names = F)


#Create separate plots for each Celltype:
celltypes <- sc_nofilter_DE_f %>% select(CellType) %>% unlist() %>% as.character() %>% unique()

sc_nofilter_DE_f <- as.data.frame(sc_nofilter_DE_f)

    #Get cell-specific DE gene vectors:
      for (c in 1:length(celltypes)) {
        
        cell <- celltypes[]
       
        cell.genes <- sc_nofilter_DE_f %>% filter(CellType == cell) %>% filter(Color_factor != "Not Significant") %>% dplyr::select(Gene) %>% unlist() %>% as.character() %>% unique()
          
        assign(paste(cell, "_heatmap_genes", sep = ""), cell.genes)
        
        rm(cell.genes)
        
      }

    
    Idents(contact_derm) <- "Semi_Detailed_Celltype"
    
    #Subset for each celltype:
        for (c in 1:length(celltypes)) {
          
          cell <- celltypes[c]
          
          cell.seurat.subset <- subset(contact_derm, idents = cell)
          cell.seurat.subset@assays$RNA@meta.features <- cell.seurat.subset@assays$RNA@meta.features[,which(str_detect(colnames(cell.seurat.subset@assays$RNA@meta.features), pattern = cell))]
          
          cell.seurat.subset$Lesion <- factor(cell.seurat.subset$Lesion)
        
          assign(paste(cell, "_data_subset", sep = ""), cell.seurat.subset)
      
        }

    #Plot for each celltype:
          
          data(geneList, package="DOSE")
          
          
          #CD4
          set.seed(100)

              CD4.heatmap <- plot_heatmap_srt(CD4_data_subset, 
                                           genes = CD4_heatmap_genes, 
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = "row", 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           cluster_type = "kmeans", 
                                           k = 6, 
                                           show_k = T,
                                           color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                           log_scale_base = NULL,
                                           ceiling = F,
                                           floor = F)
              #Gene order:
                CD4_genes <- c(2,3,   4,    6,   5,1)
                CD4_genes <- rev(CD4_genes)
                CD4_genes_reorder <- c()
                
                for (i in 1:length(CD4_genes)) {
                  int1 <- CD4_genes[i]
                  ind <- grep(paste0("^", int1, "$"), CD4.heatmap[[2]]$cluster)
                  reorder <- names(CD4.heatmap[[2]]$cluster)[ind]
                  CD4_genes_reorder <- c(CD4_genes_reorder, reorder)
                }
              
              #Heatmap Gene breaks:
                CD4_cluster_breaks <- c(4,6,5)
                CD4_cluster_breaks <- rev(CD4_cluster_breaks)
                
                CD4_gene_breaks <- c()
                
                for (i in 1:length(CD4_cluster_breaks)) {
                  int1 <- CD4_cluster_breaks[i]
                  ind <- grep(paste0("^", int1, "$"), CD4.heatmap[[2]]$cluster)
                  cluster_genes <- names(CD4.heatmap[[2]]$cluster)[ind]
                  CD4_gene_breaks <- c(CD4_gene_breaks, cluster_genes[length(cluster_genes)])
                }
                
              #Gene Ontology enrichment within gene groups:
                genes = CD4_genes_reorder
                gene_breaks = CD4_gene_breaks
                heatmap_data =CD4.heatmap
                
                      GO_terms <- list(c())
                      Gene_cluster_df <-data.frame(genes = genes,
                                                    gene_facet_group = NA)
                
                      for (s in 1:length(gene_breaks)) {
                        
                        if (s == 1) {
                          facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                        }
                        
                        if (s > 1) {
                          facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                          
                        }
                        
                        Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                        
                      }
                      
                      last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                      
                      Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                      
                      Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                              levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                      
                      
                      for (i in 1:length(gene_breaks)) {
      
                        gene <- gene_breaks[i]
                        
                        facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                        
                        cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                        
                        
                              cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                              cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                              
                              cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                              names(cluster_gene_vector) <- cluster_gene_data$genes
            
                              universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                              universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                              
                              GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                    pvalueCutoff  = 0.01,
                                                    qvalueCutoff  = 0.05,
                                                    readable      = TRUE)
                              
                              GO_result <- GO_result@result
                              
                              GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                              GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                              
                              GO_result <- GO_result %>% arrange(p.adjust)
                              
                              if (nrow(GO_result) > 3) {
                                cluster_GO_terms <- GO_result$Description[1:3]
                              }
                              
                              if (nrow(GO_result) <= 3) {
                                cluster_GO_terms <- GO_result$Description
                              }
                              
                        GO_terms[[i]] <- cluster_GO_terms
                        
                }
              
                      CD4_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                
                      
                                   plot_heatmap_srt(CD4_data_subset, 
                                                       genes = CD4_genes_reorder, 
                                                       gene_breaks = CD4_gene_breaks,
                                                       gene_labels = NULL,
                                                       gene_label_side = "left",
                                                       gene_labels_size = 3,
                                                       gene_labels_nudge = 0.3,
                                                       facet_row_descriptions = CD4_GO_terms,
                                                       facet_row_descriptions_side = "right",
                                                       facet_row_descriptions_size = 3,
                                                       facet_row_descriptions_nudge = 0.1,
                                                       type = "bulk", 
                                                       facet_by = "Semi_Detailed_Celltype",
                                                       scale_group = "Semi_Detailed_Celltype",
                                                       cluster_by = FALSE, 
                                                       gene_names = F, 
                                                       pdf_format = "tile", 
                                                       scale_by = "row",
                                                       text_angle = 90,
                                                       color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                       log_scale_base = NULL,
                                                       ceiling = 1.5,
                                                       floor = -1.5,
                                                       panel_spacing = 1,
                                                       plot_axis = TRUE)
                                  
                                   
          #CD8
                   set.seed(100)
                                   
                    CD8.heatmap <- plot_heatmap_srt(CD8_data_subset, 
                                                    genes = CD8_heatmap_genes, 
                                                    type = "bulk", 
                                                    facet_by = "Semi_Detailed_Celltype",
                                                    scale_group = "Semi_Detailed_Celltype",
                                                    cluster_by = "row", 
                                                    pdf_format = "tile", 
                                                    scale_by = "row",
                                                    text_angle = 90,
                                                    cluster_type = "kmeans", 
                                                    k = 6, 
                                                    show_k = T,
                                                    color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                    log_scale_base = NULL,
                                                    ceiling = F,
                                                    floor = F)
                # Gene Order:
                    CD8_genes <- c(2,   1,6,   3,   5,4)
                    
                    CD8_genes <- rev(CD8_genes)
                    CD8_genes_reorder <- c()
                    
                    for (i in 1:length(CD8_genes)) {
                      int1 <- CD8_genes[i]
                      ind <- grep(paste0("^", int1, "$"), CD8.heatmap[[2]]$cluster)
                      reorder <- names(CD8.heatmap[[2]]$cluster)[ind]
                      CD8_genes_reorder <- c(CD8_genes_reorder, reorder)
                    }
                
                # Heatmap Gene Breaks:    
                    CD8_cluster_breaks <- c(1,3,5)
                    
                    
                    CD8_cluster_breaks <- rev(CD8_cluster_breaks)
                
                    CD8_gene_breaks <- c()
                    
                    for (i in 1:length(CD8_cluster_breaks)) {
                      int1 <- CD8_cluster_breaks[i]
                      ind <- grep(paste0("^", int1, "$"), CD8.heatmap[[2]]$cluster)
                      cluster_genes <- names(CD8.heatmap[[2]]$cluster)[ind]
                      CD8_gene_breaks <- c(CD8_gene_breaks, cluster_genes[length(cluster_genes)])
                    }
                    
                # Gene Ontology enrichment within gene groups:
                    genes = CD8_genes_reorder
                    gene_breaks = CD8_gene_breaks
                    heatmap_data =CD8.heatmap
                    
                              GO_terms <- list(c())
                              Gene_cluster_df <-data.frame(genes = genes,
                                                           gene_facet_group = NA)
                              
                              for (s in 1:length(gene_breaks)) {
                                
                                if (s == 1) {
                                  facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                                }
                                
                                if (s > 1) {
                                  facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                                  
                                }
                                
                                Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                                
                    }
                              
                              last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                              
                              Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                              
                              Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                         levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                              
                              
                              for (i in 1:length(gene_breaks)) {
                                
                                gene <- gene_breaks[i]
                                
                                facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                                
                                cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                                
                                
                                cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                                cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                                
                                cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                                names(cluster_gene_vector) <- cluster_gene_data$genes
                                
                                universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                                universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                                
                                GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.01,
                                                      qvalueCutoff  = 0.05,
                                                      readable      = TRUE)
                                
                                GO_result <- GO_result@result
                                
                                GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                                GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                                
                                GO_result <- GO_result %>% arrange(p.adjust)
                                
                                if (nrow(GO_result) > 3) {
                                  cluster_GO_terms <- GO_result$Description[1:3]
                                }
                                
                                if (nrow(GO_result) <= 3) {
                                  cluster_GO_terms <- GO_result$Description
                                }
                                
                                GO_terms[[i]] <- cluster_GO_terms
                                
                    }
                    
                    CD8_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                    
                                        plot_heatmap_srt(CD8_data_subset, 
                                                         genes = CD8_genes_reorder, 
                                                         gene_breaks = CD8_gene_breaks,
                                                         gene_labels = NULL,
                                                         gene_label_side = "left",
                                                         gene_labels_size = 3,
                                                         gene_labels_nudge = 0.3,
                                                         facet_row_descriptions = CD8_GO_terms,
                                                         facet_row_descriptions_side = "right",
                                                         facet_row_descriptions_size = 3,
                                                         facet_row_descriptions_nudge = 0.1,
                                                         type = "bulk", 
                                                         facet_by = "Semi_Detailed_Celltype",
                                                         scale_group = "Semi_Detailed_Celltype",
                                                         cluster_by = FALSE, 
                                                         gene_names = F, 
                                                         pdf_format = "tile", 
                                                         scale_by = "row",
                                                         text_angle = 90,
                                                         color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                         log_scale_base = NULL,
                                                         ceiling = 1.5,
                                                         floor = -1.5,
                                                         panel_spacing = 1,
                                                         plot_axis = TRUE)
                                    
            
          #KRT
                  set.seed(100)
                                        
                  KRT.heatmap <- plot_heatmap_srt(KRT_data_subset, 
                                                  genes = KRT_heatmap_genes, 
                                                  type = "bulk", 
                                                  facet_by = "Semi_Detailed_Celltype",
                                                  scale_group = "Semi_Detailed_Celltype",
                                                  cluster_by = "row", 
                                                  pdf_format = "tile", 
                                                  scale_by = "row",
                                                  text_angle = 90,
                                                  cluster_type = "kmeans", 
                                                  k = 6, 
                                                  show_k = T,
                                                  color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                  log_scale_base = NULL,
                                                  ceiling = F,
                                                  floor = F)
                  
              # Gene Order:
                  KRT_genes <- c(2,   5,4,    1,     3,6)   

                  KRT_genes <- rev(KRT_genes)
                  KRT_genes_reorder <- c()
                
                  for (i in 1:length(KRT_genes)) {
                    int1 <- KRT_genes[i]
                    ind <- grep(paste0("^", int1, "$"), KRT.heatmap[[2]]$cluster)
                    reorder <- names(KRT.heatmap[[2]]$cluster)[ind]
                    KRT_genes_reorder <- c(KRT_genes_reorder, reorder)
                  }
                  
              # Heatmap Gene Breaks:    
                  KRT_cluster_breaks <- c(5,  1,  3)
                  
                  KRT_cluster_breaks <- rev(KRT_cluster_breaks)
                  
                  KRT_gene_breaks <- c()
                  
                  for (i in 1:length(KRT_cluster_breaks)) {
                    int1 <- KRT_cluster_breaks[i]
                    ind <- grep(paste0("^", int1, "$"), KRT.heatmap[[2]]$cluster)
                    cluster_genes <- names(KRT.heatmap[[2]]$cluster)[ind]
                    KRT_gene_breaks <- c(KRT_gene_breaks, cluster_genes[length(cluster_genes)])
                    
                  }
              
              # Gene Ontology enrichment within gene groups:
                  genes = KRT_genes_reorder
                  gene_breaks = KRT_gene_breaks
                  heatmap_data = KRT.heatmap
                  
                            GO_terms <- list(c())
                            Gene_cluster_df <-data.frame(genes = genes,
                                                         gene_facet_group = NA)
                            
                            for (s in 1:length(gene_breaks)) {
                              
                              if (s == 1) {
                                facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                              }
                              
                              if (s > 1) {
                                facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                                
                              }
                              
                              Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                              
                            }
                            
                            last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                            
                            Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                            
                            Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                       levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                            
                            
                            for (i in 1:length(gene_breaks)) {
                              
                              gene <- gene_breaks[i]
                              
                              facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                              
                              cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                              
                              
                              cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                              cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                              
                              cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                              names(cluster_gene_vector) <- cluster_gene_data$genes
                              
                              universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                              universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                              
                              GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                    pvalueCutoff  = 0.01,
                                                    qvalueCutoff  = 0.05,
                                                    readable      = TRUE)
                              
                              GO_result <- GO_result@result
                              
                              GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                              GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                              
                              GO_result <- GO_result %>% arrange(p.adjust)
                              
                              if (nrow(GO_result) > 3) {
                                cluster_GO_terms <- GO_result$Description[1:3]
                              }
                              
                              if (nrow(GO_result) <= 3) {
                                cluster_GO_terms <- GO_result$Description
                              }
                              
                              GO_terms[[i]] <- cluster_GO_terms
                              
                            }
                            
                    KRT_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                            
                            plot_heatmap_srt(KRT_data_subset, 
                                             genes = KRT_genes_reorder, 
                                             gene_breaks = KRT_gene_breaks,
                                             gene_labels = NULL,
                                             gene_label_side = "left",
                                             gene_labels_size = 3,
                                             gene_labels_nudge = 0.3,
                                             facet_row_descriptions = KRT_GO_terms,
                                             facet_row_descriptions_side = "right",
                                             facet_row_descriptions_size = 3,
                                             facet_row_descriptions_nudge = 0.1,
                                             type = "bulk", 
                                             facet_by = "Semi_Detailed_Celltype",
                                             scale_group = "Semi_Detailed_Celltype",
                                             cluster_by = FALSE, 
                                             gene_names = F, 
                                             pdf_format = "tile", 
                                             scale_by = "row",
                                             text_angle = 90,
                                             color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                             log_scale_base = NULL,
                                             ceiling = 1.5,
                                             floor = -1.5,
                                             panel_spacing = 1,
                                             plot_axis = TRUE)
                            
                                          
          #MAC1
                          set.seed(100)
                            
                          Mac1.heatmap <- plot_heatmap_srt(Mac1_data_subset, 
                                                          genes = Mac1_heatmap_genes, 
                                                          type = "bulk", 
                                                          facet_by = "Semi_Detailed_Celltype",
                                                          scale_group = "Semi_Detailed_Celltype",
                                                          cluster_by = "row", 
                                                          pdf_format = "tile", 
                                                          scale_by = "row",
                                                          text_angle = 90,
                                                          cluster_type = "kmeans", 
                                                          k = 6, 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                          log_scale_base = NULL,
                                                          ceiling = F,
                                                          floor = F)
                  # Gene order:
                          Mac1_genes <- c(2,   4,5,    1,   6,3)

                          
                          Mac1_genes <- rev(Mac1_genes)
                          Mac1_genes_reorder <- c()
                          
                          for (i in 1:length(Mac1_genes)) {
                            int1 <- Mac1_genes[i]
                            ind <- grep(paste0("^", int1, "$"), Mac1.heatmap[[2]]$cluster)
                            reorder <- names(Mac1.heatmap[[2]]$cluster)[ind]
                            Mac1_genes_reorder <- c(Mac1_genes_reorder, reorder)
                          }
                          
                  # Heatmap Gene Breaks:    
                          Mac1_cluster_breaks <- c(4, 1, 6)

                          
                          Mac1_cluster_breaks <- rev(Mac1_cluster_breaks)
                          
                          Mac1_gene_breaks <- c()
                          
                          for (i in 1:length(Mac1_cluster_breaks)) {
                            int1 <- Mac1_cluster_breaks[i]
                            ind <- grep(paste0("^", int1, "$"), Mac1.heatmap[[2]]$cluster)
                            cluster_genes <- names(Mac1.heatmap[[2]]$cluster)[ind]
                            Mac1_gene_breaks <- c(Mac1_gene_breaks, cluster_genes[length(cluster_genes)])
                            
                          }
                          
                  # Gene Ontology enrichment within gene groups:
                          genes = Mac1_genes_reorder
                          gene_breaks = Mac1_gene_breaks
                          heatmap_data = Mac1.heatmap
                          
                          GO_terms <- list(c())
                          Gene_cluster_df <-data.frame(genes = genes,
                                                       gene_facet_group = NA)
                          
                          for (s in 1:length(gene_breaks)) {
                            
                            if (s == 1) {
                              facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                            }
                            
                            if (s > 1) {
                              facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                              
                            }
                            
                            Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                            
                          }
                          
                          last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                          
                          Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                          
                          Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                     levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                          
                          
                          for (i in 1:length(gene_breaks)) {
                            
                            gene <- gene_breaks[i]
                            
                            facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                            
                            cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                            
                            
                            cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                            cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                            
                            cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                            names(cluster_gene_vector) <- cluster_gene_data$genes
                            
                            universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                            universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                            
                            GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                  pvalueCutoff  = 0.01,
                                                  qvalueCutoff  = 0.05,
                                                  readable      = TRUE)
                            
                            GO_result <- GO_result@result
                            
                            GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                            GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                            
                            GO_result <- GO_result %>% arrange(p.adjust)
                            
                            if (nrow(GO_result) > 3) {
                              cluster_GO_terms <- GO_result$Description[1:3]
                            }
                            
                            if (nrow(GO_result) <= 3) {
                              cluster_GO_terms <- GO_result$Description
                            }
                            
                            GO_terms[[i]] <- cluster_GO_terms
                            
                          }
                          
                          Mac1_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                          
                          plot_heatmap_srt(Mac1_data_subset, 
                                           genes = Mac1_genes_reorder, 
                                           gene_breaks = Mac1_gene_breaks,
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           facet_row_descriptions = Mac1_GO_terms,
                                           facet_row_descriptions_side = "right",
                                           facet_row_descriptions_size = 3,
                                           facet_row_descriptions_nudge = 0.1,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           gene_names = F, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                           log_scale_base = NULL,
                                           ceiling = 1.5,
                                           floor = -1.5,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          
                          
              #MAC2
                  set.seed(100)
                          
                  Mac2.heatmap <- plot_heatmap_srt(Mac2_data_subset, 
                                                   genes = Mac2_heatmap_genes, 
                                                   type = "bulk", 
                                                   facet_by = "Semi_Detailed_Celltype",
                                                   scale_group = "Semi_Detailed_Celltype",
                                                   cluster_by = "row", 
                                                   pdf_format = "tile", 
                                                   scale_by = "row",
                                                   text_angle = 90,
                                                   cluster_type = "kmeans", 
                                                   k = 6, 
                                                   show_k = T,
                                                   color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                   log_scale_base = NULL,
                                                   ceiling = F,
                                                   floor = F)
                  
                  # Gene order:
                      Mac2_genes <- c(2,     1,     3,5,   6,4)
                      #Mac2_genes <- c(3,     2,     4,   1,5,6)
                  
                      Mac2_genes <- rev(Mac2_genes)
                      Mac2_genes_reorder <- c()
                      
                      for (i in 1:length(Mac2_genes)) {
                        int1 <- Mac2_genes[i]
                        ind <- grep(paste0("^", int1, "$"), Mac2.heatmap[[2]]$cluster)
                        reorder <- names(Mac2.heatmap[[2]]$cluster)[ind]
                        Mac2_genes_reorder <- c(Mac2_genes_reorder, reorder)
                      }
                  
                  # Heatmap Gene Breaks:    
                      Mac2_cluster_breaks <- c(1,3,6)
                      #Mac2_cluster_breaks <- c(2,4,1)
                      
                      Mac2_cluster_breaks <- rev(Mac2_cluster_breaks)
                      
                      Mac2_gene_breaks <- c()
                      
                      for (i in 1:length(Mac2_cluster_breaks)) {
                        int1 <- Mac2_cluster_breaks[i]
                        ind <- grep(paste0("^", int1, "$"), Mac2.heatmap[[2]]$cluster)
                        cluster_genes <- names(Mac2.heatmap[[2]]$cluster)[ind]
                        Mac2_gene_breaks <- c(Mac2_gene_breaks, cluster_genes[length(cluster_genes)])
                        
                      }
                  
                  # Gene Ontology enrichment within gene groups:
                      genes = Mac2_genes_reorder
                      gene_breaks = Mac2_gene_breaks
                      heatmap_data = Mac2.heatmap
                      
                      GO_terms <- list(c())
                      Gene_cluster_df <-data.frame(genes = genes,
                                                   gene_facet_group = NA)
                      
                      for (s in 1:length(gene_breaks)) {
                        
                        if (s == 1) {
                          facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                        }
                        
                        if (s > 1) {
                          facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                          
                        }
                        
                        Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                        
                      }
                      
                      last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                      
                      Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                      
                      Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                 levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                      
                      for (i in 1:length(gene_breaks)) {
                        
                        gene <- gene_breaks[i]
                        
                        facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                        
                        cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                        
                        
                        cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                        cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                        
                        cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                        names(cluster_gene_vector) <- cluster_gene_data$genes
                        
                        universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                        universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                        
                        GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                              pvalueCutoff  = 0.01,
                                              qvalueCutoff  = 0.05,
                                              readable      = TRUE)
                        
                        GO_result <- GO_result@result
                        
                        GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                        GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                        
                        GO_result <- GO_result %>% arrange(p.adjust)
                        
                        if (nrow(GO_result) > 3) {
                          cluster_GO_terms <- GO_result$Description[1:3]
                        }
                        
                        if (nrow(GO_result) <= 3) {
                          cluster_GO_terms <- GO_result$Description
                        }
                        
                        GO_terms[[i]] <- cluster_GO_terms
                        
                      }
                      
                      Mac2_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                      
                                    plot_heatmap_srt(Mac2_data_subset, 
                                                     genes = Mac2_genes_reorder, 
                                                     gene_breaks = Mac2_gene_breaks,
                                                     gene_labels = NULL,
                                                     gene_label_side = "left",
                                                     gene_labels_size = 3,
                                                     gene_labels_nudge = 0.3,
                                                     facet_row_descriptions = Mac2_GO_terms,
                                                     facet_row_descriptions_side = "right",
                                                     facet_row_descriptions_size = 3,
                                                     facet_row_descriptions_nudge = 0.1,
                                                     type = "bulk", 
                                                     facet_by = "Semi_Detailed_Celltype",
                                                     scale_group = "Semi_Detailed_Celltype",
                                                     cluster_by = FALSE, 
                                                     gene_names = F, 
                                                     pdf_format = "tile", 
                                                     scale_by = "row",
                                                     text_angle = 90,
                                                     color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                     log_scale_base = NULL,
                                                     ceiling = 1.5,
                                                     floor = -1.5,
                                                     panel_spacing = 1,
                                                     plot_axis = TRUE)
                          
            #LC
                    set.seed(100)
                                    
                    LC.heatmap <- plot_heatmap_srt(LC_data_subset, 
                                                     genes = LC_heatmap_genes, 
                                                     type = "bulk", 
                                                     facet_by = "Semi_Detailed_Celltype",
                                                     scale_group = "Semi_Detailed_Celltype",
                                                     cluster_by = "row", 
                                                     pdf_format = "tile", 
                                                     scale_by = "row",
                                                     text_angle = 90,
                                                     cluster_type = "kmeans", 
                                                     k = 6, 
                                                     show_k = T,
                                                     color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                     log_scale_base = NULL,
                                                     ceiling = F,
                                                     floor = F)
                    
                  # Gene order:
                        LC_genes <- c(6,3,4,   1,   5,2)
                        LC_genes <- rev(LC_genes)
                        LC_genes_reorder <- c()
                        
                        for (i in 1:length(LC_genes)) {
                          int1 <- LC_genes[i]
                          ind <- grep(paste0("^", int1, "$"), LC.heatmap[[2]]$cluster)
                          reorder <- names(LC.heatmap[[2]]$cluster)[ind]
                          LC_genes_reorder <- c(LC_genes_reorder, reorder)
                        }
                    
                  # Heatmap Gene Breaks:
                    
                      LC_cluster_breaks <- c(1, 5)
                      LC_cluster_breaks <- rev(LC_cluster_breaks)
                      
                      LC_gene_breaks <- c()
                      
                      for (i in 1:length(LC_cluster_breaks)) {
                        int1 <- LC_cluster_breaks[i]
                        ind <- grep(paste0("^", int1, "$"), LC.heatmap[[2]]$cluster)
                        cluster_genes <- names(LC.heatmap[[2]]$cluster)[ind]
                        LC_gene_breaks <- c(LC_gene_breaks, cluster_genes[length(cluster_genes)])
                      }
                      
                  # Gene Ontology enrichment within gene groups:
                    genes = LC_genes_reorder
                    gene_breaks = LC_gene_breaks
                    heatmap_data = LC.heatmap
                    
                            GO_terms <- list(c())
                            Gene_cluster_df <-data.frame(genes = genes,
                                                         gene_facet_group = NA)
                            
                            for (s in 1:length(gene_breaks)) {
                              
                              if (s == 1) {
                                facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                              }
                              
                              if (s > 1) {
                                facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                                
                              }
                              
                              Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                              
                            }
                            
                            last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                            
                            Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                            
                            Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                       levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                            
                            
                            for (i in 1:length(gene_breaks)) {
                              
                              gene <- gene_breaks[i]
                              
                              facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                              
                              cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                              
                              
                              cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                              cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                              
                              cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                              names(cluster_gene_vector) <- cluster_gene_data$genes
                              
                              universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                              universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                              
                              GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                    pvalueCutoff  = 0.01,
                                                    qvalueCutoff  = 0.05,
                                                    readable      = TRUE)
                              
                              GO_result <- GO_result@result
                              
                              GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                              GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                              
                              GO_result <- GO_result %>% arrange(p.adjust)
                              
                              if (nrow(GO_result) > 3) {
                                cluster_GO_terms <- GO_result$Description[1:3]
                              }
                              
                              if (nrow(GO_result) <= 3) {
                                cluster_GO_terms <- GO_result$Description
                              }
                              
                              GO_terms[[i]] <- cluster_GO_terms
                              
                            }
                            
                    LC_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                                  
                                  plot_heatmap_srt(LC_data_subset, 
                                                   genes = LC_genes_reorder, 
                                                   gene_breaks = LC_gene_breaks,
                                                   gene_labels = NULL,
                                                   gene_label_side = "left",
                                                   gene_labels_size = 3,
                                                   gene_labels_nudge = 0.3,
                                                   facet_row_descriptions = LC_GO_terms,
                                                   facet_row_descriptions_side = "right",
                                                   facet_row_descriptions_size = 3,
                                                   facet_row_descriptions_nudge = 0.1,
                                                   type = "bulk", 
                                                   facet_by = "Semi_Detailed_Celltype",
                                                   scale_group = "Semi_Detailed_Celltype",
                                                   cluster_by = FALSE, 
                                                   gene_names = F, 
                                                   pdf_format = "tile", 
                                                   scale_by = "row",
                                                   text_angle = 90,
                                                   color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                   log_scale_base = NULL,
                                                   ceiling = 1.5,
                                                   floor = -1.5,
                                                   panel_spacing = 1,
                                                   plot_axis = TRUE)
                                            
                                    
          #MEL
              set.seed(100)
                                  
              MEL.heatmap <- plot_heatmap_srt(MEL_data_subset, 
                                             genes = MEL_heatmap_genes, 
                                             type = "bulk", 
                                             facet_by = "Semi_Detailed_Celltype",
                                             scale_group = "Semi_Detailed_Celltype",
                                             cluster_by = "row", 
                                             pdf_format = "tile", 
                                             scale_by = "row",
                                             text_angle = 90,
                                             cluster_type = "kmeans", 
                                             k = 6, 
                                             show_k = T,
                                             color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                             log_scale_base = NULL,
                                             ceiling = F,
                                             floor = F)
              
              # Gene order:
                    MEL_genes <- c(3,2,     5,     6,  4,1)
                    MEL_genes <- rev(MEL_genes)
                    MEL_genes_reorder <- c()
                    
                    for (i in 1:length(MEL_genes)) {
                      int1 <- MEL_genes[i]
                      ind <- grep(paste0("^", int1, "$"), MEL.heatmap[[2]]$cluster)
                      reorder <- names(MEL.heatmap[[2]]$cluster)[ind]
                      MEL_genes_reorder <- c(MEL_genes_reorder, reorder)
                    }
              
              # Heatmap Gene Breaks:
                    MEL_cluster_breaks <- c(5,6,4)
                    MEL_cluster_breaks <- rev(MEL_cluster_breaks)
                    
                    MEL_gene_breaks <- c()
                    
                    for (i in 1:length(MEL_cluster_breaks)) {
                      int1 <- MEL_cluster_breaks[i]
                      ind <- grep(paste0("^", int1, "$"), MEL.heatmap[[2]]$cluster)
                      cluster_genes <- names(MEL.heatmap[[2]]$cluster)[ind]
                      MEL_gene_breaks <- c(MEL_gene_breaks, cluster_genes[length(cluster_genes)])
                      
                    }
              
              # Gene Ontology enrichment within gene groups:
                    genes = MEL_genes_reorder
                    gene_breaks = MEL_gene_breaks
                    heatmap_data = MEL.heatmap
                    
                              GO_terms <- list(c())
                              Gene_cluster_df <-data.frame(genes = genes,
                                                           gene_facet_group = NA)
                              
                              for (s in 1:length(gene_breaks)) {
                                
                                if (s == 1) {
                                  facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                                }
                                
                                if (s > 1) {
                                  facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                                  
                                }
                                
                                Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                                
                    }
                              
                              last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                              
                              Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                              
                              Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                         levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                              
                              
                              for (i in 1:length(gene_breaks)) {
                                
                                gene <- gene_breaks[i]
                                
                                facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                                
                                cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                                
                                
                                cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                                cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                                
                                cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                                names(cluster_gene_vector) <- cluster_gene_data$genes
                                
                                universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                                universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                                
                                GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.01,
                                                      qvalueCutoff  = 0.05,
                                                      readable      = TRUE)
                                
                                GO_result <- GO_result@result
                                
                                GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                                GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                                
                                GO_result <- GO_result %>% arrange(p.adjust)
                                
                                if (nrow(GO_result) > 3) {
                                  cluster_GO_terms <- GO_result$Description[1:3]
                                }
                                
                                if (nrow(GO_result) <= 3) {
                                  cluster_GO_terms <- GO_result$Description
                                }
                                
                                GO_terms[[i]] <- cluster_GO_terms
                                
                    }
                              
                    MEL_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                              
                              plot_heatmap_srt(MEL_data_subset, 
                                               genes = MEL_genes_reorder, 
                                               gene_breaks = MEL_gene_breaks,
                                               gene_labels = NULL,
                                               gene_label_side = "left",
                                               gene_labels_size = 3,
                                               gene_labels_nudge = 0.3,
                                               facet_row_descriptions = MEL_GO_terms,
                                               facet_row_descriptions_side = "right",
                                               facet_row_descriptions_size = 3,
                                               facet_row_descriptions_nudge = 0.1,
                                               type = "bulk", 
                                               facet_by = "Semi_Detailed_Celltype",
                                               scale_group = "Semi_Detailed_Celltype",
                                               cluster_by = FALSE, 
                                               gene_names = F, 
                                               pdf_format = "tile", 
                                               scale_by = "row",
                                               text_angle = 90,
                                               color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                               log_scale_base = NULL,
                                               ceiling = 1.5,
                                               floor = -1.5,
                                               panel_spacing = 1,
                                               plot_axis = TRUE)
                    
                              
        #Treg
            set.seed(100)
                              
            Treg.heatmap <- plot_heatmap_srt(Treg_data_subset, 
                                            genes = Treg_heatmap_genes, 
                                            type = "bulk", 
                                            facet_by = "Semi_Detailed_Celltype",
                                            scale_group = "Semi_Detailed_Celltype",
                                            cluster_by = "row", 
                                            pdf_format = "tile", 
                                            scale_by = "row",
                                            text_angle = 90,
                                            cluster_type = "kmeans", 
                                            k = 6, 
                                            show_k = T,
                                            color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                            log_scale_base = NULL,
                                            ceiling = F,
                                            floor = F)
            
            # Gene Order:
                  Treg_genes <- c(3,2,1,     6,    5,     4)
                  Treg_genes <- rev(Treg_genes)
                  Treg_genes_reorder <- c()
                  
                  for (i in 1:length(Treg_genes)) {
                    int1 <- Treg_genes[i]
                    ind <- grep(paste0("^", int1, "$"), Treg.heatmap[[2]]$cluster)
                    reorder <- names(Treg.heatmap[[2]]$cluster)[ind]
                    Treg_genes_reorder <- c(Treg_genes_reorder, reorder)
                  }
            
            # Heatmap Gene Breaks:
                  Treg_cluster_breaks <- c(6,5,4)
                  Treg_cluster_breaks <- rev(Treg_cluster_breaks)
                  
                  Treg_gene_breaks <- c()
                  
                  for (i in 1:length(Treg_cluster_breaks)) {
                    int1 <- Treg_cluster_breaks[i]
                    ind <- grep(paste0("^", int1, "$"), Treg.heatmap[[2]]$cluster)
                    cluster_genes <- names(Treg.heatmap[[2]]$cluster)[ind]
                    Treg_gene_breaks <- c(Treg_gene_breaks, cluster_genes[length(cluster_genes)])  
                  }
              
            # Gene Ontology enrichment within gene groups:
                  genes = Treg_genes_reorder
                  gene_breaks = Treg_gene_breaks
                  heatmap_data = Treg.heatmap
                  
                              GO_terms <- list(c())
                              Gene_cluster_df <-data.frame(genes = genes,
                                                           gene_facet_group = NA)
                              
                              for (s in 1:length(gene_breaks)) {
                                
                                if (s == 1) {
                                  facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                                }
                                
                                if (s > 1) {
                                  facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                                  
                                }
                                
                                Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                                
                              }
                              
                              last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                              
                              Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                              
                              Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                         levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                              
                              
                              for (i in 1:length(gene_breaks)) {
                                
                                gene <- gene_breaks[i]
                                
                                facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                                
                                cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                                
                                
                                cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                                cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                                
                                cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                                names(cluster_gene_vector) <- cluster_gene_data$genes
                                
                                universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                                universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                                
                                GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                      pvalueCutoff  = 0.01,
                                                      qvalueCutoff  = 0.05,
                                                      readable      = TRUE)
                                
                                GO_result <- GO_result@result
                                
                                GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                                GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                                
                                GO_result <- GO_result %>% arrange(p.adjust)
                                
                                if (nrow(GO_result) > 3) {
                                  cluster_GO_terms <- GO_result$Description[1:3]
                                }
                                
                                if (nrow(GO_result) <= 3) {
                                  cluster_GO_terms <- GO_result$Description
                                }
                                
                                GO_terms[[i]] <- cluster_GO_terms
                                
                              }
                              
                      Treg_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                              
                              plot_heatmap_srt(Treg_data_subset, 
                                               genes = Treg_genes_reorder, 
                                               gene_breaks = Treg_gene_breaks,
                                               gene_labels = NULL,
                                               gene_label_side = "left",
                                               gene_labels_size = 3,
                                               gene_labels_nudge = 0.3,
                                               facet_row_descriptions = Treg_GO_terms,
                                               facet_row_descriptions_side = "right",
                                               facet_row_descriptions_size = 3,
                                               facet_row_descriptions_nudge = 0.1,
                                               type = "bulk", 
                                               facet_by = "Semi_Detailed_Celltype",
                                               scale_group = "Semi_Detailed_Celltype",
                                               cluster_by = FALSE, 
                                               gene_names = F, 
                                               pdf_format = "tile", 
                                               scale_by = "row",
                                               text_angle = 90,
                                               color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                               log_scale_base = NULL,
                                               ceiling = 1.5,
                                               floor = -1.5,
                                               panel_spacing = 1,
                                               plot_axis = TRUE)
                  
            
        
                                              
          #NK
            set.seed(100)
                              
            NK.heatmap <- plot_heatmap_srt(NK_data_subset, 
                                             genes = NK_heatmap_genes, 
                                             type = "bulk", 
                                             facet_by = "Semi_Detailed_Celltype",
                                             scale_group = "Semi_Detailed_Celltype",
                                             cluster_by = "row", 
                                             pdf_format = "tile", 
                                             scale_by = "row",
                                             text_angle = 90,
                                             cluster_type = "kmeans", 
                                             k = 6, 
                                             show_k = T,
                                             color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                             log_scale_base = NULL,
                                             ceiling = F,
                                             floor = F)
            
            # Gene Order:
                  NK_genes <- c(2,6,   1,    5,   3,4)
                  NK_genes <- rev(NK_genes)
                  NK_genes_reorder <- c()
                  
                  for (i in 1:length(NK_genes)) {
                    int1 <- NK_genes[i]
                    ind <- grep(paste0("^", int1, "$"), NK.heatmap[[2]]$cluster)
                    reorder <- names(NK.heatmap[[2]]$cluster)[ind]
                    NK_genes_reorder <- c(NK_genes_reorder, reorder)
                  }
                  
            # Heatmap Gene Breaks:
                      
                NK_cluster_breaks <- c(1,5,3)
                NK_cluster_breaks <- rev(NK_cluster_breaks)
                
                NK_gene_breaks <- c()
                
                for (i in 1:length(NK_cluster_breaks)) {
                  int1 <- NK_cluster_breaks[i]
                  ind <- grep(paste0("^", int1, "$"), NK.heatmap[[2]]$cluster)
                  cluster_genes <- names(NK.heatmap[[2]]$cluster)[ind]
                  NK_gene_breaks <- c(NK_gene_breaks, cluster_genes[length(cluster_genes)])
                }
                                    
            # Gene Ontology enrichment within gene groups:
                genes = NK_genes_reorder
                gene_breaks = NK_gene_breaks
                heatmap_data = NK.heatmap
                
                        GO_terms <- list(c())
                        Gene_cluster_df <-data.frame(genes = genes,
                                                     gene_facet_group = NA)
                        
                        for (s in 1:length(gene_breaks)) {
                          
                          if (s == 1) {
                            facet_section_genes <- genes[1:which(genes == gene_breaks[s])]
                          }
                          
                          if (s > 1) {
                            facet_section_genes <- genes[c(which(genes == gene_breaks[s-1])+1):which(genes == gene_breaks[s])]
                            
                          }
                          
                          Gene_cluster_df[which(Gene_cluster_df$genes %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                          
                        }
                        
                        last_facet_section_genes <- genes[c(which(genes == gene_breaks[length(gene_breaks)])+1):length(genes)]
                        
                        Gene_cluster_df[which(Gene_cluster_df$genes %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
                        
                        Gene_cluster_df$gene_facet_group <- factor(Gene_cluster_df$gene_facet_group,
                                                                   levels = rev(unique(Gene_cluster_df$gene_facet_group)))
                        
                        
                        for (i in 1:length(gene_breaks)) {
                          
                          gene <- gene_breaks[i]
                          
                          facet_section <- Gene_cluster_df %>% filter(genes == gene) %>% dplyr::select(gene_facet_group) %>% unlist()
                          
                          cluster_genes <- Gene_cluster_df %>% filter(gene_facet_group == facet_section) %>% dplyr::select(genes) %>% unlist() %>% as.character()
                          
                          
                          cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% cluster_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
                          cluster_gene_data <- cluster_gene_data %>% arrange(desc(max_expression))
                          
                          cluster_gene_vector <- as.numeric(cluster_gene_data$max_expression)
                          names(cluster_gene_vector) <- cluster_gene_data$genes
                          
                          universe_id <- bitr(names(cluster_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                          universe_id$max_expression <- cluster_gene_vector[match(universe_id$SYMBOL, names(cluster_gene_vector))]
                          
                          GO_result <- enrichGO(gene = universe_id$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.01,
                                                qvalueCutoff  = 0.05,
                                                readable      = TRUE)
                          
                          GO_result <- GO_result@result
                          
                          GO_result$set_length <- sapply(GO_result$GeneRatio, FUN = function(x) as.numeric(str_split(x, pattern = "/")[[1]])[2])
                          GO_result$Gene_frac <- GO_result$Count/GO_result$set_length
                          
                          GO_result <- GO_result %>% arrange(p.adjust)
                          
                          if (nrow(GO_result) > 3) {
                            cluster_GO_terms <- GO_result$Description[1:3]
                          }
                          
                          if (nrow(GO_result) <= 3) {
                            cluster_GO_terms <- GO_result$Description
                          }
                          
                          GO_terms[[i]] <- cluster_GO_terms
                          
                        }
                        
                  NK_GO_terms <- c("", rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
                        
                            plot_heatmap_srt(NK_data_subset, 
                                             genes = NK_genes_reorder, 
                                             gene_breaks = NK_gene_breaks,
                                             gene_labels = NULL,
                                             gene_label_side = "left",
                                             gene_labels_size = 3,
                                             gene_labels_nudge = 0.3,
                                             facet_row_descriptions = NK_GO_terms,
                                             facet_row_descriptions_side = "right",
                                             facet_row_descriptions_size = 3,
                                             facet_row_descriptions_nudge = 0.1,
                                             type = "bulk", 
                                             facet_by = "Semi_Detailed_Celltype",
                                             scale_group = "Semi_Detailed_Celltype",
                                             cluster_by = FALSE, 
                                             gene_names = F, 
                                             pdf_format = "tile", 
                                             scale_by = "row",
                                             text_angle = 90,
                                             color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                             log_scale_base = NULL,
                                             ceiling = 1.5,
                                             floor = -1.5,
                                             panel_spacing = 1,
                                             plot_axis = TRUE)
                                
                                
  ## Redo all celltypes merged with only genes that show reasonable pattern above:
      filtered_DE_genes <- c()
                                            
          for (c in 1:length(celltypes)) {
          
            gene_breaks <- get(paste(celltypes[c], "_gene_breaks", sep = ""))
            
            full_DE_genes <- get(paste(celltypes[c], "_genes_reorder", sep = ""))

            last_facet_section_genes <- full_DE_genes[c(which(full_DE_genes == gene_breaks[length(gene_breaks)])+1):length(full_DE_genes)]
            
            filtered_cell_DE_genes <- full_DE_genes[-which(full_DE_genes %in% last_facet_section_genes)]
            
            filtered_DE_genes <- unique(c(filtered_DE_genes, filtered_cell_DE_genes))
            
          }
      
      
      plot_heatmap_srt(contact_derm,
                       genes = filtered_DE_genes,
                       type = "bulk",
                       facet_by = "Semi_Detailed_Celltype",
                       scale_group = "Semi_Detailed_Celltype",
                       scale_by = "row",
                       text_angle = 90,
                       text_sizes = c(20, 10, 7, 10, 5, 5, 5),
                       title = "Filtered DE genes",
                       color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                       log_scale_base = NULL,
                       ceiling = F,
                       floor = F,
                       gene_names = F)
      
      
# Gene Pattern Analysis:

      for (c in 1:length(celltypes)) {
        
        gene_breaks <- get(paste(celltypes[c], "_gene_breaks", sep = ""))
        
        full_DE_genes <- get(paste(celltypes[c], "_genes_reorder", sep = ""))
         
        gene_cell_pattern_df <- data.frame(gene = full_DE_genes,
                                             celltype = celltypes[c],
                                             gene_facet_group = NA)  
        
        
        
              for (s in 1:length(gene_breaks)) {
                
                if (s == 1) {
                  facet_section_genes <- full_DE_genes[1:which(full_DE_genes == gene_breaks[s])]
                }
                
                if (s > 1) {
                  facet_section_genes <- full_DE_genes[c(which(full_DE_genes == gene_breaks[s-1])+1):which(full_DE_genes == gene_breaks[s])]
                  
                }
                
                gene_cell_pattern_df[which(gene_cell_pattern_df$gene %in% facet_section_genes), "gene_facet_group"] <- as.numeric(s)
                
              }
        
              last_facet_section_genes <- full_DE_genes[c(which(full_DE_genes == gene_breaks[length(gene_breaks)])+1):length(full_DE_genes)]
              gene_cell_pattern_df[which(gene_cell_pattern_df$gene %in% last_facet_section_genes), "gene_facet_group"] <- c(length(gene_breaks)+1)
              
        
              if (c == 1) {
                gene_DE_pattern_df <- gene_cell_pattern_df
              }
              if (c > 1) {
                gene_DE_pattern_df <- bind_rows(gene_DE_pattern_df, gene_cell_pattern_df)
              }
        
      }
      
      gene_DE_pattern_df <- gene_DE_pattern_df %>% mutate(lesion_group = case_when(gene_facet_group == 1 ~ "Allergy",
                                                                                   gene_facet_group == 2 ~ "Irritant",
                                                                                   gene_facet_group == 3 ~ "NL",
                                                                                   gene_facet_group == 4 ~ "Acetone"))
      
      gene_DE_pattern_df <- gene_DE_pattern_df %>% group_by(gene) %>% summarize(count = n(),
                                                                  celltypes = str_c(sort(unique(celltype)), collapse = ", "),
                                                                  celltype_count = length(unique(celltype)),
                                                                  lesion_groups = str_c(sort(unique(lesion_group)), collapse = ", "),
                                                                  lesion_count = length(unique(lesion_group)))
      

      ggplot(gene_DE_pattern_df, aes(x = lesion_count)) + geom_bar() + theme_classic() + scale_x_discrete(limits = seq(1:4)) + scale_y_continuous(limits = c(0, 500))
      ggplot(gene_DE_pattern_df, aes(x = celltype_count)) + geom_bar() + theme_classic() + scale_x_discrete(limits = seq(1:9)) + scale_y_continuous(limits = c(0, 400))
      
      ggplot(gene_DE_pattern_df, aes(x = celltype_count, y = lesion_count, color = lesion_groups)) + geom_jitter() + theme_classic() + scale_x_discrete(limits = seq(1:9)) + scale_y_continuous(limits = c(0, 5)) + labs(title = "DE Gene Specificity Analysis") + theme(plot.title = element_text(hjust = 0.5))
      

# GO Enrichment Pattern Analysis:
      
      for (c in 1:length(celltypes)) {
        
        GO_terms <- get(paste(celltypes[c], "_GO_terms", sep = ""))
          
              if (any(GO_terms == "")) {
                GO_terms <- GO_terms[-which(GO_terms == "")]
              }
      
              GO_terms <- str_split(GO_terms, pattern = "\n")
              GO_term_df <- as.data.frame(GO_terms)
              colnames(GO_term_df) <- seq(1:ncol(GO_term_df))
              GO_term_df <- pivot_longer(GO_term_df, cols = c(1:ncol(GO_term_df)), names_to = "gene_facet_group", values_to = "GO_term")
              
          GO_term_df$celltype = celltypes[c]
              
              if (c == 1) {
                GO_term_pattern_df <- GO_term_df
              }
              if (c > 1) {
                GO_term_pattern_df <- bind_rows(GO_term_pattern_df, GO_term_df)
              }
            
      }
      
      GO_term_pattern_df <- GO_term_pattern_df %>% mutate(lesion_group = case_when(gene_facet_group == 1 ~ "NL",
                                                                                   gene_facet_group == 2 ~ "Irritant",
                                                                                   gene_facet_group == 3 ~ "Allergy"))

      GO_term_pattern_df <- GO_term_pattern_df %>% group_by(GO_term) %>% summarize(count = n(),
                                                                                celltypes = str_c(sort(unique(celltype)), collapse = ", "),
                                                                                celltype_count = length(unique(celltype)),
                                                                                lesion_groups = str_c(sort(unique(lesion_group)), collapse = ", "),
                                                                                lesion_count = length(unique(lesion_group)))
      
      
      ggplot(GO_term_pattern_df, aes(x = lesion_count)) + geom_bar() + theme_classic() + scale_x_discrete(limits = seq(1:4)) + scale_y_continuous(limits = c(0, 500))
      ggplot(GO_term_pattern_df, aes(x = celltype_count)) + geom_bar() + theme_classic() + scale_x_discrete(limits = seq(1:9)) + scale_y_continuous(limits = c(0, 400))
      
      ggplot(GO_term_pattern_df, aes(x = celltype_count, y = lesion_count, color = lesion_groups)) + geom_jitter(width = 0.25, height = 0.25) + theme_classic() + scale_x_discrete(limits = seq(1:9)) + scale_y_continuous(limits = c(0, 4)) + labs(title = "GO Term Specificity Analysis") + theme(plot.title = element_text(hjust = 0.5))
      
      
# Gene LogFC comparison:
      DE.data.wider <- DE.data.wider %>% arrange()
      
      DE.data.wider <- full_join(DE.data.wider, gene_DE_pattern_df, by = c("Gene" = "gene"))
      
      DE.data.wider$lesion_count[is.na(DE.data.wider$lesion_count)] <- 0
      DE.data.wider$lesion_count <- factor(DE.data.wider$lesion_count,
                                           levels = rev(c(0,1,2,3,4)))
      
      ggplot(DE.data.wider, aes(x = FoldChange_sc, y = `FoldChange_sc patient batched; patch patients only`, color = lesion_count)) + 
        facet_grid(CellType ~ Comparison, scales = "free") +
        geom_point(alpha = 0.6) + 
        theme_classic() + 
        labs(title = "Gene LogFC comparison using patient batch DE covariate vs. non-batched all data") + 
        theme(plot.title = element_text(hjust=0.5)) +
        scale_x_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn) +
        scale_y_continuous(breaks = c(-10,-5, -2, 1, 2,5, 10,20), labels=c(-10,-5, -2, 1, 2,5, 10,20), trans = tn) +
        geom_vline(aes(xintercept = 2^0.5), color = "#BD3535", linetype = "dashed") +
        geom_vline(aes(xintercept = -2^0.5), color = "#BD3535", linetype = "dashed") +
        geom_hline(aes(yintercept = 2^0.5), color = "#BD3535", linetype = "dashed") +
        geom_hline(aes(yintercept = -2^0.5), color = "#BD3535", linetype = "dashed") +
        scale_color_manual(values = rev(c("#747474", "#127546", "#184EA5", "#FFA122", "#E54127")))
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
##########################################                            
                # Old Code:
##########################################                            
                                    
                                
                                
                                
                                
                                
                                
                                
                                
                                


                      #With Gene Order:
                      set.seed(100)
                      de_genes <- plot_heatmap_srt(contact_derm, 
                                                    genes = DE.heatmap.genes, 
                                                    type = "bulk", 
                                                    facet_by = "Semi_Detailed_Celltype",
                                                    scale_group = "Semi_Detailed_Celltype",
                                                    cluster_by = "row", 
                                                    pdf_format = "tile", 
                                                    scale_by = "row",
                                                    text_angle = 90,
                                                    cluster_type = "kmeans", 
                                                    k = 10, 
                                                    show_k = T,
                                                    color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                                    log_scale_base = NULL,
                                                    ceiling = F,
                                                    floor = F)

                      
                      
                      int_de_genes <- c()
                      int_de_genes <- rev(int_de_genes)
                      de_genes_reorder <- c()
                      
                      for (i in 1:length(int_de_genes)) {
                        int1 <- int_de_genes[i]
                        ind <- grep(paste0("^", int1, "$"), de_genes[[2]]$cluster)
                        reorder <- names(de_genes[[2]]$cluster)[ind]
                        de_genes_reorder <- c(de_genes_reorder, reorder)
                      }
                      
                      plot_heatmap_srt(contact_derm, 
                                       genes = de_genes_reorder, 
                                       type = "bulk", 
                                       facet_by = "Semi_Detailed_Celltype",
                                       scale_group = "Semi_Detailed_Celltype",
                                       cluster_by = FALSE, 
                                       gene_names = F, 
                                       pdf_format = "tile", 
                                       scale_by = "row",
                                       text_angle = 90,
                                       color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                       log_scale_base = NULL,
                                       ceiling = 1.5,
                                       floor = -1.5)


                    #Scaled individually by celltype:
                      plot_heatmap_srt(contact_derm, 
                                       genes = de_genes_reorder, 
                                       type = "bulk",
                                       facet_by = "Semi_Detailed_Celltype",
                                       scale_group = "Semi_Detailed_Celltype",
                                       cluster_by = FALSE, 
                                       pdf_format = "tile", 
                                       scale_by = "row", 
                                       text_angle = 90,
                                       text_sizes = c(20, 10, 7, 10, 5, 5, 5),
                                       gene_names = F, 
                                       col_names = TRUE, 
                                       title = "DE Genes",
                                       color_pal = colorRampPalette(c("#97A0CD","#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                                       log_scale_base = NULL,
                                       ceiling = 2,
                                       floor = -2)

plot_violin_srt(contact_derm, gene = "CXCL13", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype")




genes.for.heatmap.pval.fc.filtered[1]
genes.to.filter




#Barplot of Significant Gene Counts:
Up_genecount <- DE.data %>% group_by(Comparison, CellType) %>% 
  filter(Color_factor == "Significantly Up") %>% 
  summarise(Significantly_Up = n(),
            Mean.FoldChange = mean(logfc2fc(logFC))) %>% 
  gather(key = "direction", value = "gene_count", -c(Comparison, CellType, Mean.FoldChange))

Down_genecount <- DE.data %>% group_by(Comparison, CellType) %>% 
  filter(Color_factor == "Significantly Down") %>% 
  summarise(Significantly_Down = n(),
            Mean.FoldChange = mean(logfc2fc(logFC))) %>% 
  gather(key = "direction", value = "gene_count", -c(Comparison, CellType, Mean.FoldChange))

significant_genecount <- bind_rows(Up_genecount, Down_genecount)
significant_genecount$direction <- factor(significant_genecount$direction, levels = c("Significantly_Up", "Significantly_Down"))



ggplot(significant_genecount, aes(x = CellType, y = gene_count, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE)) +
  facet_grid(~Comparison) +
  theme_classic() + 
  labs(title = "Lesional vs. Nonlesional DE Gene Counts for padj < 0.01", x = "CellType", y = "DE Gene Count", fill = "Gene Regulation") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_fill_manual(values = c("coral1", "cornflowerblue"))    

ggplot(significant_genecount, aes(x = CellType, y = abs(Mean.FoldChange), fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE)) +
  facet_grid(~Comparison) +
  theme_classic() + 
  labs(title = "Effect size of DE Gene Counts for padj < 0.01", x = "CellType", y = "Abs(Mean Fold Change)", fill = "Gene Regulation") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_fill_manual(values = c("coral1", "cornflowerblue")) +
  scale_y_continuous(breaks = c(0, 2,4,6,8,10,12), labels=c(0, 2,4,6,8,10,12))





#Dotplots for top most significant genes for each Patient (BY Semi_Detailed_CellType):    

#Get Aggregate bulk gene data:
contact_derm <- calc_agg_bulk_srt(contact_derm, aggregate_by = c("Lesion", "Semi_Detailed_Celltype", "Patient"))

bulk_aggregate_terms = c("Lesion", "Semi_Detailed_Celltype", "Patient")


gene.data <- contact_derm@assays$RNA@meta.features[, grep("bulk", colnames(contact_derm@assays$RNA@meta.features))]

#remove "dummy cells" that have all zeros:

gene.data <- gene.data[,-which(str_detect(colnames(gene.data), pattern = "num__cells__1__"))]

gene.data_long <- tidyr::gather(gene.data, key = "group", "Expression", 
                                1:ncol(gene.data), factor_key = "TRUE")
gene.data_long$Gene <- rep(factor(rownames(gene.data), levels = rownames(gene.data)), 
                           ncol(gene.data))

for (i in 1:length(bulk_aggregate_terms)) {
  facs <- as.character(unique(contact_derm@meta.data[, bulk_aggregate_terms[i]]))
  
  gene.data_long$new_col <- NA
  
  for (j in 1:length(facs)) {
    gene.data_long$new_col[which(str_detect(string = gene.data_long$group, pattern = facs[j]))] <- str_extract(string = gene.data_long$group, pattern = facs[j])[!is.na(str_extract(string = gene.data_long$group, pattern = facs[j]))]
  }
  
  gene.data_long <- gene.data_long %>% dplyr::rename(!!bulk_aggregate_terms[i] := new_col)
  
}

gene.data_long$Gene <- as.character(gene.data_long$Gene)          

#Filter Mitochondrial and ribosomal genes out of list:
  filtered.gene.data <- gene.data_long
  filtered.gene.data <- gene.data_long[-which(str_detect(gene.data_long$Gene, pattern = "^RPL|^RPS|^MT")),]
  
filtered.DE.data <- DE.data[-which(str_detect(DE.data$Gene, pattern = "^RPL|^RPS|^MT")),]

# Choose Genes by Top DE filter:
top.results <- filtered.DE.data %>% filter(FoldChange_Direction == "Up") %>% group_by(Comparison, CellType) %>% top_n(n = 6, wt = -log(padj)) %>% ungroup
top.results <- top.results[,which(colnames(top.results) %in% c("Gene", "CellType", "Comparison"))]
top.results <- top.results %>% gather(key = "DE.perspective", value = "Lesion",-c(Gene, CellType))
top.results$Lesion <- str_replace_all(top.results$Lesion, pattern = " ", replacement = "_")



#Get Data only for Comparisons and Reference for each CellType:
comparison.data <- gene.data_long[0,]

for (i in 1:nrow(top.results)) {
  
  lesion <- as.character(top.results[i,"Lesion"])
  lesion <- str_replace_all(lesion, pattern = " ", replacement = "_")
  CellType <- as.character(unlist(top.results[i, "CellType"]))
  gene <- as.character(unlist(top.results[i, "Gene"]))
  #Edit for Celltype Factor
  data <- gene.data_long %>% filter(Lesion == lesion, Semi_Detailed_Celltype == CellType, Gene == gene)
  
  data$Comparison.lesion <- lesion
  
  comparison.data <- full_join(comparison.data, data)                  
  
}
#Edit for Celltype Factor
comparison.data <- comparison.data %>% dplyr::rename(CellType = Semi_Detailed_Celltype)

top.data <- full_join(top.results, comparison.data)

reference.data <- gene.data_long[0,]

for (i in 1:nrow(top.results)) {
  
  lesion <- Reference
  CellType <- as.character(unlist(top.results[i, "CellType"]))
  gene <- as.character(unlist(top.results[i, "Gene"]))
  comparison.lesion <- as.character(top.results[i,"Lesion"])
  
  
  #Edit for Celltype Factor
  data <- gene.data_long %>% filter(Lesion == lesion, Semi_Detailed_Celltype == CellType, Gene == gene)
  
  data$Comparison.lesion <- comparison.lesion
  reference.data <- full_join(reference.data, data)                  
  
}


#Edit for Celltype Factor
reference.data <- reference.data %>% dplyr::rename(CellType = Semi_Detailed_Celltype)

reference.data <- reference.data[,which(colnames(reference.data) %in% c("Gene", "Expression", "CellType", "Patient", "Comparison.lesion"))]
reference.data$DE.perspective <- "Reference"

top.data <- full_join(top.data, reference.data)

top.data$DE.perspective <- factor(top.data$DE.perspective, levels = c("Reference", "Comparison"))
top.data$Comparison.lesion <- factor(top.data$Comparison.lesion, levels = c("Irritant", "Day2_Allergy", "Day4_Allergy"))
top.data$CellType <- factor(top.data$CellType, levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))
top.data$Comparison.lesion <- factor(str_replace_all(top.data$Comparison.lesion, pattern = "_", replacement = " "),
                                     levels = str_replace_all(levels(top.data$Comparison.lesion), pattern = "_", replacement = " "))



ggplot(top.data, aes(x = Gene, y = Expression, color = DE.perspective)) +
  facet_wrap(Comparison.lesion ~ CellType, scales = "free") +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_classic() + 
  labs(title = "Lesional vs NonLesional aggregate bulk patient values for Top Genes by padj", x = "Genes", y = "Aggregate bulk CPM", color = "Lesion") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_color_manual(values = c("gray30", "coral1")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  theme(axis.text.x = element_text(size = 10, angle = 90))
  


#### Heatmaps:
contact_derm@meta.data$Semi_Detailed_Celltype <- factor(contact_derm@meta.data$Semi_Detailed_Celltype,
                                                        levels = c("CD8","CD4", "Treg", "NK", "LC", "Mac1", "Mac2", "Bcell", "KRT", "MEL"))


contact_derm <- calc_agg_bulk_srt(contact_derm, aggregate_by = c("Lesion", "Semi_Detailed_Celltype"), group_by = "Lesion")

# Choose Genes by Top DE filter:
Significant.genes <- filtered.DE.data %>% filter(FoldChange_Direction == "Up" & padj < (10^-15) & Comparison %in% c("Day2 Allergy", "Day4 Allergy"))
Significant.genes <- unique(Significant.genes$Gene)
    
  #Filter further by means in other control groups:
    aggregate.counts <- contact_derm@assays$RNA@meta.features[,which(str_detect(string = colnames(contact_derm@assays$RNA@meta.features), pattern = "bulk"))]
    
    rowsums <- rowSums(aggregate.counts)
    aggregate.counts <- aggregate.counts[which(rowsums > 0),]
    
    aggregate.counts2 <- t(apply(aggregate.counts, 1, scale))
    colnames(aggregate.counts2) <- colnames(aggregate.counts)
    aggregate.counts <- aggregate.counts2
    colnames(aggregate.counts) <- sub("__num_.*", "", colnames(aggregate.counts))
    
    aggregate.counts_lng <- tidyr::gather(as.data.frame(aggregate.counts), key = "group", "Expression", 
                                  1:ncol(aggregate.counts), factor_key = "TRUE")
    aggregate.counts_lng$genes <- rep(factor(rownames(aggregate.counts), levels = rownames(aggregate.counts)), 
                              ncol(aggregate.counts))
    
    aggregate.counts_lng$Lesion <- sapply(aggregate.counts_lng$group, FUN = function(x) unlist(str_split(x, pattern = "__")[[1]])[1])
    aggregate.counts_lng$Celltype <- sapply(aggregate.counts_lng$group, FUN = function(x) unlist(str_split(x, pattern = "__")[[1]])[2])
    

    gene.filter.df <- data.frame(gene = unique(aggregate.counts_lng$genes),
                                 max.celltype.expres.mean = NA)
    
    gene.filter.df <- gene.filter.df %>% filter(gene %in% Significant.genes)
    
    for (g in 1:length(gene.filter.df$gene)) {
      
      gene <- gene.filter.df$gene[g]
      
      gene.data <- filter(aggregate.counts_lng, genes == gene & Lesion %in% c("Nonlesional", "Irritant", "Acetone_Vehicle"))
      gene.filter.df[g,"max.celltype.expres.mean"] <- as.numeric(max(gene.data$Expression))
      
      }
    
    genes.for.heatmap <- gene.filter.df %>% filter(max.celltype.expres.mean <2) %>% dplyr::select(gene) %>% unlist() %>% as.character()
    
    plot_heatmap_srt(contact_derm,
                 genes = genes.for.heatmap,
                 type = "bulk",
                 facet_by = "Semi_Detailed_Celltype",
                 calc_agg_bulk_group.term = "Lesion",
                 text_angle = 90,
                 text_sizes = c(20, 10, 7, 10, 5, 5, 5))




        #With Gene Order:
        set.seed(100)
        DE.gene.plot <- plot_heatmap_srt(contact_derm, 
                                          genes = Significant.genes, 
                                          type = "bulk",
                                          facet_by = "Semi_Detailed_Celltype",
                                          calc_agg_bulk_group.term = "Lesion",
                                          cluster_by = "row", 
                                          pdf_format = "tile", 
                                          scale_by = "row",
                                          cluster_type = "kmeans", 
                                          k = 15, 
                                          show_k = T, 
                                          text_angle = 90,
                                          text_sizes = c(20, 10, 7, 10, 5, 5, 5))
        
        int_cytokines <- c()
        int_cytokines <- rev(int_cytokines)
        cytokines_reorder <- c()
        
        for (i in 1:length(int_cytokines)) {
          int1 <- int_cytokines[i]
          ind <- grep(paste0("^", int1, "$"), cytokines[[2]]$cluster)
          reorder <- names(cytokines[[2]]$cluster)[ind]
          cytokines_reorder <- c(cytokines_reorder, reorder)
        }
        
        plot_heatmap_srt(contact_derm, 
                         genes = cytokines_reorder, 
                         type = "bulk", 
                         cluster_by = FALSE, 
                         pdf_format = "tile", 
                         scale_by = "row", 
                         gene_names = T, 
                         group_names = TRUE, 
                         title = "Cytokines")
