#######################################################
# DE Analysis: Concatenation of DE results:
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
library(formattable)

#Load data:

setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load data:
load(file = "Data/contact_derm_filtered3.Rdata")

              #Set Celltypes as factors:
                      contact_derm@meta.data$Basic_Celltype <- factor(contact_derm@meta.data$Basic_Celltype,
                                                                      levels = c("KRT", "MEL", "APC", "LYMPH"))
                      
                      contact_derm@meta.data$Semi_Detailed_Celltype <- factor(contact_derm@meta.data$Semi_Detailed_Celltype, 
                                                                              levels = c("KRT", "MEL", "LC","Myeloid","pDC", "CD4", "Treg", "CD8", "NK"))
                      
                      contact_derm@meta.data$Detailed_Celltype <- factor(contact_derm@meta.data$Detailed_Celltype, 
                                                                         levels = c("KRT-b1", "KRT-b2", "KRT-sp", "KRT-g", "KRT-77","KRT-mucl","KRT-wr", "Mel1", "Mel2", "LC", "Myeloid_1", "Myeloid_ccl22", "M2_mac", "cDC1","pDC", 
                                                                                    "CD4_conv1", "CD4_hsp", "CD4_cd161", "CD4_tcm", "CD4_crem", "Treg", "CD8_CD4_active","CD8", "NK_ctl", "NK_areg"))

                      #DE Subset:
                      DE.input = contact_derm
                      
                      DE.input@meta.data <- DE.input@meta.data %>% mutate(DE_inclusion_factor = case_when(Lesion == "Nonlesional" ~ "yes",
                                                                                                          Lesion == "Acetone_Vehicle" & Lesion_visual_score < 1 ~ "yes",
                                                                                                          Lesion == "Irritant" ~ "yes",
                                                                                                          Lesion == "Day2_Allergy" ~ "yes",
                                                                                                          Lesion == "Day4_Allergy" ~ "yes",
                                                                                                          TRUE ~ "no"))
                      
                      
                      
                      Idents(DE.input) <- "DE_inclusion_factor"
                      DE.input <- subset(DE.input, idents = "yes")
                      Idents(DE.input) <- "Semi_Detailed_Celltype"
                      
                      plot_violin_bar_srt(DE.input, gene = "CXCL14", color_by = "Lesion", size = 0.5,facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors, plot_patient_pos_frac = FALSE)
                      plot_tsne_gene_srt(DE.input, gene = "CXCR3")

# First merge DE files to single dataframe:

            #Single Lesion comparisons:
                  file.path = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_single_lesions/"
                  
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
                    
                    
                    assign(paste(CellType,"_",Comparison,"_vs_", Reference,".data", sep = ""), DE.data)
                    
                    DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,"_vs_", Reference,".data", sep = ""))
                    
                  }
                  
                  DE.data <- purrr::reduce(mget(DE.data.list), bind_rows)

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
                      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & abs(FoldChange) >= 2^0.5 ~ "Significantly Up", 
                                                                             FoldChange_Direction == "Down" & padj <= 0.01 & abs(FoldChange) >= 2^0.5 ~ "Significantly Down",
                                                                             padj > 0.01 ~ "Not Significant",
                                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
                      
                      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
                      
                      
                      DE.data$Comparison <- factor(DE.data$Comparison, 
                                                   levels = c("Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
                      
                      DE.data$CellType <- factor(DE.data$CellType, 
                                                 levels = levels(contact_derm@meta.data$Semi_Detailed_Celltype))
                      
                  #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt", sep = ""), sep = "\t")
                  
          #Lesion combo comparisons:
                  
                  file.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_combo_lesions/"
                  
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
                    
                    #Clarify what the "Controls" were:
                    if (Reference == "Control") {
                      Reference = "Nonlesional/Acetone/Irritant"
                    }
                    if (Reference == "Inflammatory") {
                      Reference = "Irritant/Day2_Allergy/Day4_Allergy"
                    }
                    
                    #Add labels:
                    DE.data$CellType <- CellType
                    DE.data$Comparison <- Comparison
                    DE.data$Reference <- Reference
                    
                    
                    assign(paste(CellType,"_",Comparison,".data", sep = ""), DE.data)
                    
                    DE.data.list <- c(DE.data.list, paste(CellType,"_",Comparison,".data", sep = ""))
                    
                  }
                  
                  DE.data <- purrr::reduce(mget(DE.data.list), bind_rows)
                  
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
                      DE.data <- DE.data %>% mutate(Color_factor = case_when(FoldChange_Direction == "Up" & padj <= 0.01 & abs(FoldChange) >= 2^0.5 ~ "Significantly Up", 
                                                                             FoldChange_Direction == "Down" & padj <= 0.01 & abs(FoldChange) >= 2^0.5 ~ "Significantly Down",
                                                                             padj > 0.01 ~ "Not Significant",
                                                                             abs(FoldChange) < 2^0.5 ~ "Not Significant"))
                      
                      DE.data$Color_factor <- factor(DE.data$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
                      
                      DE.data$Comparison <- factor(DE.data$Comparison, 
                                                   levels = c("Nonlesional","Day2_Allergy", "Day4_Allergy"))
                      
                      DE.data$CellType <- factor(DE.data$CellType, 
                                                 levels = levels(contact_derm@meta.data$Semi_Detailed_Celltype))
                  
                  #write.table(DE.data, paste(file.path,"Merged_DE_File/Semi_Detailed_CellType_DEcombo.merged.txt", sep = ""), sep = "\t")
                  
    single.lesion.DE <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_single_lesions/Merged_DE_File/Semi_Detailed_CellType_DE.merged.txt")
    combo.lesion.DE <- read.table("/pi/john.harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_combo_lesions/Merged_DE_File/Semi_Detailed_CellType_DEcombo.merged.txt")
    
    contact_derm_DE <- bind_rows(single.lesion.DE, combo.lesion.DE)    
         
    #write.table(contact_derm_DE, "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE.txt", sep = "\t")


#Filter DE results to exclude comparisons with very few cell numbers and unconvincingly low expression values:
    contact_derm_DE <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE.txt")
    
    single.lesion.DE <- contact_derm_DE %>% filter(!Reference %in% c("Nonlesional/Acetone/Irritant", "Irritant/Day2_Allergy/Day4_Allergy"))             
    combo.lesion.DE <- contact_derm_DE %>% filter(Reference %in% c("Nonlesional/Acetone/Irritant"))            
    combo.lesion.DE2 <- contact_derm_DE %>% filter(Reference %in% c("Irritant/Day2_Allergy/Day4_Allergy"))            
    
    
    contact_derm@meta.data$Lesion <- as.character(contact_derm@meta.data$Lesion)
    contact_derm@meta.data <- contact_derm@meta.data %>% mutate(Lesion_combos = case_when(Lesion %in% c("Nonlesional", "Irritant", "Acetone_Vehicle") ~ "Nonlesional/Acetone/Irritant",
                                                                                  Lesion == "Day2_Allergy" ~ "Day2_Allergy",
                                                                                  Lesion == "Day4_Allergy" ~ "Day4_Allergy"))
    
    contact_derm@meta.data <- contact_derm@meta.data %>% mutate(Lesion_combos2 = case_when(Lesion %in% c("Irritant", "Day2_Allergy", "Day4_Allergy") ~ "Irritant/Day2_Allergy/Day4_Allergy",
                                                                                   Lesion == "Nonlesional" ~ "Nonlesional"))
    
    contact_derm@meta.data$Lesion <- factor(contact_derm@meta.data$Lesion,
                                               levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
    
    contact_derm@meta.data$Lesion_combos <- factor(contact_derm@meta.data$Lesion_combos,
                                               levels = c("Nonlesional/Acetone/Irritant", "Day2_Allergy", "Day4_Allergy"))
    
    contact_derm@meta.data$Lesion_combos2 <- factor(contact_derm@meta.data$Lesion_combos2,
                                                   levels = c("Irritant/Day2_Allergy/Day4_Allergy", "Nonlesional"))
    
            #DE Subset:
            DE.input = contact_derm
            
            DE.input@meta.data <- DE.input@meta.data %>% mutate(DE_inclusion_factor = case_when(Lesion == "Nonlesional" ~ "yes",
                                                                                                Lesion == "Acetone_Vehicle" & Lesion_visual_score < 1 ~ "yes",
                                                                                                Lesion == "Irritant" ~ "yes",
                                                                                                Lesion == "Day2_Allergy" ~ "yes",
                                                                                                Lesion == "Day4_Allergy" ~ "yes",
                                                                                                TRUE ~ "no"))
            
            
            
            Idents(DE.input) <- "DE_inclusion_factor"
            DE.input <- subset(DE.input, idents = "yes")
            Idents(DE.input) <- "Semi_Detailed_Celltype"
            
    
    single.lesion.DE_filtered <- DE_filter_srt(DE_data = single.lesion.DE, single_cell_data = DE.input, sc_metadata_cols = c("Lesion", "Semi_Detailed_Celltype"))
    combo.lesion.DE_filtered <- DE_filter_srt(DE_data = combo.lesion.DE, single_cell_data = DE.input, sc_metadata_cols = c("Lesion_combos", "Semi_Detailed_Celltype"))
    combo.lesion.DE_filtered2 <- DE_filter_srt(DE_data = combo.lesion.DE2, single_cell_data = DE.input, sc_metadata_cols = c("Lesion_combos2", "Semi_Detailed_Celltype"))
    
    contact_derm_DE_filtered <- bind_rows(single.lesion.DE_filtered, combo.lesion.DE_filtered, combo.lesion.DE_filtered2)    
    #write.table(contact_derm_DE_filtered, "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE_filtered.txt", sep = "\t")
    
    
###############################################################
# DE Analysis: Visualization 
###############################################################
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
library(scales)
library(clusterProfiler)
library(aliases2entrez)
library(writexl)
library(formattable)
library(sparkline)
library(org.Hs.eg.db)
    
#Load data:

setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load data:
contact_derm_DE_filtered <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE_filtered.txt")
    
load(file = "Data/contact_derm_filtered3.Rdata")

      #Set Celltypes as factors:
          contact_derm@meta.data$Basic_Celltype <- factor(contact_derm@meta.data$Basic_Celltype,
                                                          levels = c("KRT", "MEL", "APC", "LYMPH"))
          
          contact_derm@meta.data$Semi_Detailed_Celltype <- factor(contact_derm@meta.data$Semi_Detailed_Celltype, 
                                                                  levels = c("KRT", "MEL", "LC","Myeloid","pDC", "CD4", "Treg", "CD8", "NK"))
          
          contact_derm@meta.data$Detailed_Celltype <- factor(contact_derm@meta.data$Detailed_Celltype, 
                                                             levels = c("KRT-b1", "KRT-b2", "KRT-sp", "KRT-g", "KRT-krt77","KRT-mucl1","KRT-wr", "MEL-1", "MEL-2", "LC", "Myeloid-1", "Myeloid-ccl22", "MAC-fcn1", "cDC1","pDC", 
                                                                        "CD4-conv1", "CD4-hspa6", "CD4-klrb1", "CD4-il7r", "CD4-crem", "Treg", "CD8-CD4-active","CD8", "NK-ctl", "NK-areg"))
          
          #DE Subset:
          DE.input = contact_derm
          
          DE.input@meta.data <- DE.input@meta.data %>% mutate(DE_inclusion_factor = case_when(Lesion == "Nonlesional" ~ "yes",
                                                                                              Lesion == "Acetone_Vehicle" & Lesion_visual_score < 1 ~ "yes",
                                                                                              Lesion == "Irritant" ~ "yes",
                                                                                              Lesion == "Day2_Allergy" ~ "yes",
                                                                                              Lesion == "Day4_Allergy" ~ "yes",
                                                                                              TRUE ~ "no"))
          
          
          
          Idents(DE.input) <- "DE_inclusion_factor"
          DE.input <- subset(DE.input, idents = "yes")
          Idents(DE.input) <- "Semi_Detailed_Celltype"
          

    #Mitochondrial genes didn't get filtered properly:
    if (any(str_detect(contact_derm_DE_filtered$Gene, pattern = "^MT.|^MT-"))) {
      contact_derm_DE_filtered <- contact_derm_DE_filtered[-which(str_detect(contact_derm_DE_filtered$Gene, pattern = "^MT.|^MT-")),]
    }
    
#Should we add a heterogeniety filter to only retain genes that significantly change in a certain number of patients?
          
#Shrinkage MA plot:
  
  ggplot(contact_derm_DE_filtered, aes(x = logCPM, y = logFC, color = Color_factor))+
    geom_point(size = 0.1)+
    theme_classic() +
    facet_wrap(CellType~.) +
    scale_color_manual(values = c("black", "blue","red"))
  
  plotMA <- last_plot()

#Mean Variance Plot:
  for (c in 1:length(unique(DE.input@meta.data$Semi_Detailed_Celltype))) {
      celltype = as.character(unique(DE.input@meta.data$Semi_Detailed_Celltype))[c]
      celltype.ind <- which(DE.input@meta.data$Semi_Detailed_Celltype == celltype)
      
      celltype.raw.counts <- DE.input@assays$RNA@counts[,celltype.ind]
      
      celltype.vars <- apply(celltype.raw.counts, MARGIN = 1, FUN = var)
      
      if (c == 1){
        Celltype.vars.df <- data.frame(CellType = celltype,
                                       Gene = names(celltype.vars),
                                       Var = as.numeric(celltype.vars))
      }
      if (c > 1){
        new.vars <- data.frame(CellType = celltype,
                                       Gene = names(celltype.vars),
                                       Var = as.numeric(celltype.vars))
        
        Celltype.vars.df <- bind_rows(Celltype.vars.df, new.vars)
        rm(new.vars)
      }
      
    }
    
    contact_derm_DE_filtered <- left_join(contact_derm_DE_filtered, Celltype.vars.df)
    
    
    contact_derm_DE_filtered %>% filter(Var != 0) %>%
    ggplot(., aes(x = logCPM, y = log10(Var), color = Color_factor))+
      geom_point(size = 0.5)+
      theme_classic() +
      facet_wrap(CellType~., scales = "free") +
      scale_color_manual(values = c("black", "blue","red"))
    
    mean_var_plot <- last_plot()

    #setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps")
    
    #ggsave(plotMA, width = 20, height = 20,units = "cm",filename = "plotMA.png", limitsize = FALSE) 
    #ggsave(mean_var_plot, width = 20, height = 20,units = "cm",filename = "mean_var_plot.png", limitsize = FALSE) 



                                                 
                                                                                   
   #Prep for Heatmaps:
         data(geneList, package="DOSE")
         
         DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Semi_Detailed_Celltype", "Lesion"))
         
          Idents(DE.input) <- "Semi_Detailed_Celltype"
          
          DE.input@meta.data <- DE.input@meta.data[,-which(colnames(DE.input@meta.data) == "Lesion_type")]
          
          #new function:
           plot_heatmap_srt_with_GO_analysis <- function(seurat.data, k.means.heatmap, genes_reorder, cluster_breaks, include_GSEA = FALSE) {
            
            
            gene_breaks <- c()
            
            if (length(cluster_breaks) > 0) {
              
              
              for (i in 1:length(cluster_breaks)) {
                int1 <- cluster_breaks[i]
                ind <- grep(paste0("^", int1, "$"), k.means.heatmap[[2]]$cluster)
                cluster_genes <- names(k.means.heatmap[[2]]$cluster)[ind]
                gene_breaks <- c(gene_breaks, cluster_genes[length(cluster_genes)])
              }
              
              #Gene Ontology enrichment within gene groups:
              genes = genes_reorder
              gene_breaks = gene_breaks
              heatmap_data =k.means.heatmap
              
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
              
              #Add GO term analysis for last facet section:
              
              cluster_gene_data <- heatmap_data[[1]]$data %>% filter(genes %in% last_facet_section_genes) %>% group_by(genes) %>% summarize(max_expression = max(Expression))
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
              
              GO_terms[[length(GO_terms)+1]] <- cluster_GO_terms
              
              
              GO_terms <- c( rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
            }
            
            if (length(cluster_breaks) == 0) {
              
              heatmap_data =k.means.heatmap
              
              GO_terms <- list(c())
              
              cluster_gene_data <- heatmap_data[[1]]$data %>% group_by(genes) %>% summarize(max_expression = max(Expression))
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
              
              GO_terms[[1]] <- cluster_GO_terms
              
              
              GO_terms <- c( rev(sapply(GO_terms, FUN = function(x) str_c(x, collapse = "\n"))))
              
              gene_breaks = FALSE
            } 
            
            if (include_GSEA == TRUE) {
              
                    plot_heatmap_srt(seurat.data, 
                                   genes = genes_reorder, 
                                   gene_breaks = gene_breaks,
                                   gene_labels = NULL,
                                   gene_label_side = "left",
                                   gene_labels_size = 3,
                                   gene_labels_nudge = 0.3,
                                   facet_row_descriptions = GO_terms,
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
                                   color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#F7FAEE","#F7FAEE","#FF8FAE","#EC1F1F"))(256),
                                   log_scale_base = NULL,
                                   ceiling = 1.5,
                                   floor = -1.5,
                                   panel_spacing = 1,
                                   plot_axis = TRUE)
            }
            if (include_GSEA == FALSE) {
                  plot_heatmap_srt(seurat.data, 
                                   genes = genes_reorder, 
                                   gene_breaks = gene_breaks,
                                   gene_labels = NULL,
                                   gene_label_side = "left",
                                   gene_labels_size = 3,
                                   gene_labels_nudge = 0.3,
                                   facet_row_descriptions = NULL,
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
                                   color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#F7FAEE","#F7FAEE","#FF8FAE","#EC1F1F"))(256),
                                   log_scale_base = NULL,
                                   ceiling = 1.5,
                                   floor = -1.5,
                                   panel_spacing = 1,
                                   plot_axis = TRUE)
            }
            
          }
           find_gene_breaks <- function(k.means.heatmap, genes_reorder, cluster_breaks) {
             
             gene_breaks = c()
             if (length(cluster_breaks) > 0) {
               
               
               for (i in 1:length(cluster_breaks)) {
                 int1 <- cluster_breaks[i]
                 ind <- grep(paste0("^", int1, "$"), k.means.heatmap[[2]]$cluster)
                 cluster_genes <- names(k.means.heatmap[[2]]$cluster)[ind]
                 gene_breaks <- c(gene_breaks, cluster_genes[length(cluster_genes)])
               }
               
             }
             return(gene_breaks)
           }
           
           
        #Subset for each celltype:
        for (c in 1:length(unique(DE.input@meta.data$Semi_Detailed_Celltype))) {
          
          cell <- unique(as.character(DE.input@meta.data$Semi_Detailed_Celltype))[c]
          
          cell.seurat.subset <- subset(DE.input, idents = cell)
          
          cell.seurat.subset@assays$RNA@meta.features <- cell.seurat.subset@assays$RNA@meta.features[,which(str_detect(colnames(cell.seurat.subset@assays$RNA@meta.features), pattern = cell))]
          
          cell.seurat.subset$Lesion <- factor(cell.seurat.subset$Lesion)
          
          assign(paste(cell, "_data_subset", sep = ""), cell.seurat.subset)
          
        }
    
          
           #CD4:
              #CD4_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "CD4" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
              CD4_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "CD4" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
              

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
                                                          k = round(length(CD4_heatmap_genes)^0.5, digits =0),
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = F,
                                                          floor = F)
                          
                          #Gene order:      
                          #CD4_genes <- c(3,    6,    2,5,1)
                          CD4_genes <- c(3,6,    2,5,    1) 
                          
                          #groups.to.not.plot <- c(4)
                          groups.to.not.plot <- c(4)

                          CD4_genes <- rev(CD4_genes)
                          CD4_genes_reorder <- c()
                          
                          for (i in 1:length(CD4_genes)) {
                            int1 <- CD4_genes[i]
                            ind <- grep(paste0("^", int1, "$"), CD4.heatmap[[2]]$cluster)
                            reorder <- names(CD4.heatmap[[2]]$cluster)[ind]
                            CD4_genes_reorder <- c(CD4_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #CD4_cluster_breaks <- c(6,2)
                          CD4_cluster_breaks <- c(2,1)
                          
                          CD4_cluster_breaks <- rev(CD4_cluster_breaks)
                        
                          CD4_gene_breaks <- find_gene_breaks(CD4.heatmap, CD4_genes_reorder,CD4_cluster_breaks)
                        
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), CD4.heatmap[[2]]$cluster)
                            reorder <- names(CD4.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }

                          
                              #Manual gene editing:
                                    #CD4_genes_reorder <- CD4_genes_reorder[-which(CD4_genes_reorder %in% c("HSPA6", "HSPA1B", "SELL"))]
                                    #CD4_genes_reorder <- c(c("HSPA6", "HSPA1B", "SELL"), CD4_genes_reorder)
                                    
                                    #CD4_genes_reorder <- CD4_genes_reorder[-which(CD4_genes_reorder %in% c("LMO4"))]
                                    #CD4_genes_reorder <- c(CD4_genes_reorder[1:which(CD4_genes_reorder == "LINC01871")],"LMO4", CD4_genes_reorder[(which(CD4_genes_reorder == "LINC01871")+1):length(CD4_genes_reorder)])
                                    #CD4_gene_breaks[which(CD4_gene_breaks == "LINC01871")] <- "LMO4"
                                    
                                    #genes.to.swap <- CD4_genes_reorder[which(CD4_genes_reorder == "JAML"):which(CD4_genes_reorder == "LMO4")]
                                    #genes.to.swap <- rev(genes.to.swap)
                                    #CD4_genes_reorder <- CD4_genes_reorder[-(which(CD4_genes_reorder == "JAML"):which(CD4_genes_reorder == "LMO4"))]
                                    #CD4_genes_reorder <- c(CD4_genes_reorder[1:which(CD4_genes_reorder == "SLC34A2")],genes.to.swap, CD4_genes_reorder[(which(CD4_genes_reorder == "SLC34A2")+1):length(CD4_genes_reorder)])
                                    #CD4_gene_breaks[which(CD4_gene_breaks == "LMO4")] <- "JAML"
                                    
                                    #genes.to.swap <- CD4_genes_reorder[which(CD4_genes_reorder == "GZMK"):which(CD4_genes_reorder == "CCR7")]
                                    #CD4_genes_reorder <- CD4_genes_reorder[-(which(CD4_genes_reorder == "GZMK"):which(CD4_genes_reorder == "CCR7"))]
                                    #CD4_genes_reorder <- c(CD4_genes_reorder[1:which(CD4_genes_reorder == "SLC34A2")],genes.to.swap, CD4_genes_reorder[(which(CD4_genes_reorder == "SLC34A2")+1):length(CD4_genes_reorder)])
                                    #CD4_gene_breaks[which(CD4_gene_breaks == "SLC34A2")] <- "CCR7"
                                    
                                    #CD4_genes_reorder <- CD4_genes_reorder[-which(CD4_genes_reorder %in% c("BAG3", "DNAJB1"))]
                                    #CD4_genes_reorder <- c(CD4_genes_reorder[1:which(CD4_genes_reorder == "GIMAP7")],"DNAJB1","BAG3", CD4_genes_reorder[(which(CD4_genes_reorder == "GIMAP7")+1):length(CD4_genes_reorder)])

                                    #CD4_genes_reorder <- CD4_genes_reorder[-which(CD4_genes_reorder %in% c("HNRNPA1P48"))]
                                    #CD4_genes_reorder <- c(CD4_genes_reorder[1:which(CD4_genes_reorder == "SLC34A2")],"HNRNPA1P48", CD4_genes_reorder[(which(CD4_genes_reorder == "SLC34A2")+1):length(CD4_genes_reorder)])
                               
                                    #CD4_genes_reorder <- CD4_genes_reorder[-which(CD4_genes_reorder %in% c("HSPH1"))]
                                    #CD4_genes_reorder <- c(CD4_genes_reorder[1:which(CD4_genes_reorder == "GIMAP7")],"HSPH1", CD4_genes_reorder[(which(CD4_genes_reorder == "GIMAP7")+1):length(CD4_genes_reorder)])
                                    
                                    #CD4_gene_breaks <- c("BAG3", CD4_gene_breaks)
                                    #CD4_gene_breaks <- CD4_gene_breaks[-which(CD4_gene_breaks == "JAML")]
                                  
                                  CD4_heatmap <- plot_heatmap_srt(CD4_data_subset, 
                                                   genes = CD4_genes_reorder, 
                                                   gene_breaks = CD4_gene_breaks,
                                                   gene_names = TRUE, 
                                                   gene_labels = NULL,
                                                   gene_label_side = "left",
                                                   gene_labels_size = 3,
                                                   gene_labels_nudge = 0.3,
                                                   type = "bulk", 
                                                   facet_by = "Semi_Detailed_Celltype",
                                                   scale_group = "Semi_Detailed_Celltype",
                                                   cluster_by = FALSE, 
                                                   pdf_format = "tile", 
                                                   scale_by = "row",
                                                   text_angle = 90,
                                                   color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                   log_scale_base = 10,
                                                   ceiling = F,
                                                   floor = F,
                                                   panel_spacing = 1,
                                                   plot_axis = TRUE)

                                  CD4.plot.all.labels <- last_plot()
                                  
                                  plot_heatmap_srt(CD4_data_subset, 
                                                  genes = c(CD4_genes_reorder, genes.to.not.plot),
                                                  gene_breaks = c(CD4_gene_breaks,CD4_genes_reorder[length(CD4_genes_reorder)]),
                                                  gene_names = TRUE, 
                                                  gene_labels = NULL,
                                                  gene_label_side = "left",
                                                  gene_labels_size = 3,
                                                  gene_labels_nudge = 0.3,
                                                  type = "bulk", 
                                                  facet_by = "Semi_Detailed_Celltype",
                                                  scale_group = "Semi_Detailed_Celltype",
                                                  cluster_by = FALSE, 
                                                  pdf_format = "tile", 
                                                  scale_by = "row",
                                                  text_angle = 90,
                                                  color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                  log_scale_base = 10,
                                                  ceiling = F,
                                                  floor = F,
                                                  panel_spacing = 1,
                                                  plot_axis = TRUE)
                                  
                                  CD4.plot.all.labels.full <- last_plot()
                                  
                                  plot_heatmap_srt(CD4_data_subset, 
                                                   genes = CD4_genes_reorder, 
                                                   gene_breaks = CD4_gene_breaks,
                                                   row_enrichGO = TRUE,
                                                   gene_names = FALSE, 
                                                   gene_label_side = "left",
                                                   gene_labels_size = 3,
                                                   gene_labels_nudge = 0.5,
                                                   type = "bulk", 
                                                   cluster_by = FALSE, 
                                                   pdf_format = "tile", 
                                                   scale_by = "row",
                                                   text_angle = 90,
                                                   color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                   log_scale_base = 10,
                                                   ceiling = F,
                                                   floor = F,
                                                   panel_spacing = 1,
                                                   plot_axis = FALSE)

                                  CD4.plot <- last_plot()

                  #CD8:
                          #CD8_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "CD8" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          CD8_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "CD8" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()

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
                                                          k = round(length(CD8_heatmap_genes)^0.5, digits =0), 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = F,
                                                          floor = F)
                          
                          #Gene order:
                          #CD8_genes <- c(4,3,1,    6)
                          CD8_genes <- c(4,2,     14,13,6,11,     7,12)
                          CD8_genes <- rev(CD8_genes)
                          CD8_genes_reorder <- c()
                          
                          #groups.to.not.plot <- c(2,5)
                          
                          groups.to.not.plot <- c(1,3,5,8,   9,10)
                          
                          for (i in 1:length(CD8_genes)) {
                            int1 <- CD8_genes[i]
                            ind <- grep(paste0("^", int1, "$"), CD8.heatmap[[2]]$cluster)
                            reorder <- names(CD8.heatmap[[2]]$cluster)[ind]
                            CD8_genes_reorder <- c(CD8_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #CD8_cluster_breaks <- c(6)
                          CD8_cluster_breaks <- c(14,7)
                          CD8_cluster_breaks <- rev(CD8_cluster_breaks)
                          
                          CD8_gene_breaks <- find_gene_breaks(CD8.heatmap, CD8_genes_reorder,CD8_cluster_breaks)
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), CD8.heatmap[[2]]$cluster)
                            reorder <- names(CD8.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }
                                        #Manual gene editing:
                                        #CD8_genes_reorder <- CD8_genes_reorder[-(which(CD8_genes_reorder == "GAPDH"):which(CD8_genes_reorder == "RCBTB2"))]
                                        #CD8_genes_reorder <- CD8_genes_reorder[-(which(CD8_genes_reorder == "FAM153CP"):which(CD8_genes_reorder == "GPR15"))]
                                        
                                        #CD8_genes_reorder <- CD8_genes_reorder[-which(CD8_genes_reorder %in% c("LTB","GZMA", "MYL6", "ARHGDIB","YARS1","MRPL12", "PPIF","FURIN"))]
                                        #CD8_genes_reorder <- c(c("LTB","GZMA", "MYL6", "ARHGDIB","YARS1","MRPL12", "PPIF","FURIN"), CD8_genes_reorder)
                                        
                                        #genes.to.swap <- CD8_genes_reorder[which(CD8_genes_reorder == "MMP25"):which(CD8_genes_reorder == "LINC01641")]
                                        #CD8_genes_reorder <- CD8_genes_reorder[-(which(CD8_genes_reorder == "MMP25"):which(CD8_genes_reorder == "LINC01641"))]
                                        #CD8_genes_reorder <- c(CD8_genes_reorder[1:which(CD8_genes_reorder == "TMBIM6")],genes.to.swap, CD8_genes_reorder[(which(CD8_genes_reorder == "TMBIM6")+1):length(CD8_genes_reorder)])
                                        #CD8_gene_breaks[which(CD8_gene_breaks == "TMBIM6")] <- "LINC01641"
                                        
                                        #genes.to.swap <- CD8_genes_reorder[which(CD8_genes_reorder == "EIF4G1"):which(CD8_genes_reorder == "TRA2B")]
                                        #genes.to.swap <- c(genes.to.swap, CD8_genes_reorder[which(CD8_genes_reorder == "CD69"):which(CD8_genes_reorder == "SLFN11")])
                                        #genes.to.swap <- c(genes.to.swap, CD8_genes_reorder[which(CD8_genes_reorder == "FSCN1"):which(CD8_genes_reorder == "MS4A6E")])
                                        
                                        #CD8_genes_reorder <- CD8_genes_reorder[-(which(CD8_genes_reorder %in% genes.to.swap))]
                                        #CD8_genes_reorder <- c(CD8_genes_reorder, genes.to.swap)
                                        
                                        #genes.to.swap <- CD8_genes_reorder[which(CD8_genes_reorder == "HSP90AB1"):which(CD8_genes_reorder == "TIGIT")]
                                        #CD8_genes_reorder <- CD8_genes_reorder[-(which(CD8_genes_reorder %in% genes.to.swap))]
                                        #CD8_genes_reorder <- c(CD8_genes_reorder[1:which(CD8_genes_reorder == "WDR5")],genes.to.swap, CD8_genes_reorder[(which(CD8_genes_reorder == "WDR5")+1):length(CD8_genes_reorder)])
                          
                          
                                  plot_heatmap_srt(CD8_data_subset, 
                                                   genes = CD8_genes_reorder, 
                                                   gene_breaks = CD8_gene_breaks,
                                                   gene_names = TRUE, 
                                                   gene_labels = NULL,
                                                   gene_label_side = "left",
                                                   gene_labels_size = 3,
                                                   gene_labels_nudge = 0.3,
                                                   type = "bulk", 
                                                   facet_by = "Semi_Detailed_Celltype",
                                                   scale_group = "Semi_Detailed_Celltype",
                                                   cluster_by = FALSE, 
                                                   pdf_format = "tile", 
                                                   scale_by = "row",
                                                   text_angle = 90,
                                                   color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                   log_scale_base = 10,
                                                   ceiling = F,
                                                   floor = F,
                                                   panel_spacing = 1,
                                                   plot_axis = TRUE)
                                  
                                  CD8.plot.all.labels <- last_plot()
                                  
                                  plot_heatmap_srt(CD8_data_subset, 
                                                   genes = c(CD8_genes_reorder, genes.to.not.plot),
                                                   gene_breaks = c(CD8_gene_breaks, CD8_genes_reorder[length(CD8_genes_reorder)]),
                                                   gene_names = TRUE, 
                                                   gene_labels = NULL,
                                                   gene_label_side = "left",
                                                   gene_labels_size = 3,
                                                   gene_labels_nudge = 0.3,
                                                   type = "bulk", 
                                                   facet_by = "Semi_Detailed_Celltype",
                                                   scale_group = "Semi_Detailed_Celltype",
                                                   cluster_by = FALSE, 
                                                   pdf_format = "tile", 
                                                   scale_by = "row",
                                                   text_angle = 90,
                                                   color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                   log_scale_base = 10,
                                                   ceiling = F,
                                                   floor = F,
                                                   panel_spacing = 1,
                                                   plot_axis = TRUE)
                                  
                                  CD8.plot.all.labels.full <- last_plot()
                                  
                                  plot_heatmap_srt(CD8_data_subset, 
                                                   genes = CD8_genes_reorder, 
                                                   gene_breaks = CD8_gene_breaks,
                                                   gene_names = FALSE, 
                                                   row_enrichGO = TRUE,
                                                   gene_label_side = "left",
                                                   gene_labels_size = 3,
                                                   gene_labels_nudge = 0.5,
                                                   type = "bulk", 
                                                   cluster_by = FALSE, 
                                                   pdf_format = "tile", 
                                                   scale_by = "row",
                                                   text_angle = 90,
                                                   color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                   log_scale_base = 10,
                                                   ceiling = F,
                                                   floor = F,
                                                   panel_spacing = 1,
                                                   plot_axis = FALSE)
                                  
                                  CD8.plot <- last_plot()
                                  

                #KRT:
                          #KRT_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "KRT" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          KRT_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "KRT" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
                          
                                
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
                                                          k = round(length(KRT_heatmap_genes)^0.5, digits =0), 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = F,
                                                          floor = F)
                          
                          #Gene order:
                          #KRT_genes <- c(1,    2,6,5,   4)
                          KRT_genes <- c(1,    2,6,   5,4)
                          KRT_genes <- rev(KRT_genes)
                          KRT_genes_reorder <- c()
                          
                          groups.to.not.plot <- c(3)
                          
                          for (i in 1:length(KRT_genes)) {
                            int1 <- KRT_genes[i]
                            ind <- grep(paste0("^", int1, "$"), KRT.heatmap[[2]]$cluster)
                            reorder <- names(KRT.heatmap[[2]]$cluster)[ind]
                            KRT_genes_reorder <- c(KRT_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #KRT_cluster_breaks <- c(2, 4)
                          KRT_cluster_breaks <- c(2, 5)
                          
                          KRT_cluster_breaks <- rev(KRT_cluster_breaks)
                          
                          KRT_gene_breaks <- find_gene_breaks(KRT.heatmap, KRT_genes_reorder,KRT_cluster_breaks)
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), KRT.heatmap[[2]]$cluster)
                            reorder <- names(KRT.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }                
                                                    
                          plot_heatmap_srt(KRT_data_subset, 
                                           genes = KRT_genes_reorder, 
                                           gene_breaks = KRT_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          KRT.plot.all.labels <- last_plot()
                          
                          
                          plot_heatmap_srt(KRT_data_subset, 
                                           genes = c(KRT_genes_reorder, genes.to.not.plot),
                                           gene_breaks = c(KRT_gene_breaks,KRT_genes_reorder[length(KRT_genes_reorder)]),
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          KRT.plot.all.labels.full <- last_plot()
                          
                          plot_heatmap_srt(KRT_data_subset, 
                                           genes = KRT_genes_reorder, 
                                           gene_breaks = KRT_gene_breaks,
                                           gene_names = FALSE, 
                                           row_enrichGO = TRUE,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.5,
                                           type = "bulk", 
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = FALSE)
                          
                          KRT.plot <- last_plot()
                          
              #MEL:
                          #MEL_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "MEL" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          MEL_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "MEL" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
                          
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
                                                          k = round(length(MEL_heatmap_genes)^0.5, digits =0), 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = F,
                                                          floor = F)
                          
                          #Gene order:
                          #MEL_genes <- c(5,7,     1)
                          MEL_genes <- c(7,10,6,12,     2,        1)
                          
                          MEL_genes <- rev(MEL_genes)
                          MEL_genes_reorder <- c()
                          
                          #groups.to.not.plot <- c(2,3,4,6,8)
                          groups.to.not.plot <- c(3,4,5,8,9,11,13,14,15,16,17,18,19)
                          
                          for (i in 1:length(MEL_genes)) {
                            int1 <- MEL_genes[i]
                            ind <- grep(paste0("^", int1, "$"), MEL.heatmap[[2]]$cluster)
                            reorder <- names(MEL.heatmap[[2]]$cluster)[ind]
                            MEL_genes_reorder <- c(MEL_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #MEL_cluster_breaks <- c(1)
                          MEL_cluster_breaks <- c(2,1)
                          
                          MEL_cluster_breaks <- rev(MEL_cluster_breaks)
                          
                          MEL_gene_breaks <- find_gene_breaks(MEL.heatmap, MEL_genes_reorder,MEL_cluster_breaks)
                          
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), MEL.heatmap[[2]]$cluster)
                            reorder <- names(MEL.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }         
                          
                          
                          plot_heatmap_srt(MEL_data_subset, 
                                           genes = MEL_genes_reorder, 
                                           gene_breaks = MEL_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          MEL.plot.all.labels <- last_plot()
                          
                          
                          plot_heatmap_srt(MEL_data_subset, 
                                           genes = c(MEL_genes_reorder, genes.to.not.plot),
                                           gene_breaks = c(MEL_gene_breaks,MEL_genes_reorder[length(MEL_genes_reorder)]),
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          MEL.plot.all.labels.full <- last_plot()
                          
                          plot_heatmap_srt(MEL_data_subset, 
                                           genes = MEL_genes_reorder, 
                                           gene_breaks = MEL_gene_breaks,
                                           gene_names = FALSE, 
                                           row_enrichGO = TRUE,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.5,
                                           type = "bulk", 
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = FALSE)
                          
                          MEL.plot <- last_plot()
                          
              #Myeloid:
                          #Myeloid_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "Myeloid" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          Myeloid_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "Myeloid" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
                          
                          set.seed(100)
                          Myeloid.heatmap <- plot_heatmap_srt(Myeloid_data_subset, 
                                                          genes = Myeloid_heatmap_genes, 
                                                          type = "bulk", 
                                                          facet_by = "Semi_Detailed_Celltype",
                                                          scale_group = "Semi_Detailed_Celltype",
                                                          cluster_by = "row", 
                                                          pdf_format = "tile", 
                                                          scale_by = "row",
                                                          text_angle = 90,
                                                          cluster_type = "kmeans", 
                                                          k = round(length(Myeloid_heatmap_genes)^0.5, digits =0), 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = F,
                                                          floor = F)
                          
                          #Gene order:
                          #Myeloid_genes <- c(1,6,    7,3,5,4,    8 )
                          Myeloid_genes <- c(10,         1,2,3,8,6,      7,11,14)
                          
                          Myeloid_genes <- rev(Myeloid_genes)
                          Myeloid_genes_reorder <- c()
                          
                          groups.to.not.plot <- c(4,5,9,13,12)
                          
                          for (i in 1:length(Myeloid_genes)) {
                            int1 <- Myeloid_genes[i]
                            ind <- grep(paste0("^", int1, "$"), Myeloid.heatmap[[2]]$cluster)
                            reorder <- names(Myeloid.heatmap[[2]]$cluster)[ind]
                            Myeloid_genes_reorder <- c(Myeloid_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #Myeloid_cluster_breaks <- c(7,8)
                          Myeloid_cluster_breaks <- c(1,7)
                          
                          Myeloid_cluster_breaks <- rev(Myeloid_cluster_breaks)
                          
                          
                          Myeloid_gene_breaks <- find_gene_breaks(Myeloid.heatmap, Myeloid_genes_reorder,Myeloid_cluster_breaks)
                          
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), Myeloid.heatmap[[2]]$cluster)
                            reorder <- names(Myeloid.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }         
                          
                                  #Manual gene editing:
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "STMN1"):which(Myeloid_genes_reorder == "H4C3")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "DSP")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "DSP")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "IRF1"):which(Myeloid_genes_reorder == "NR4A2")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "GLIPR2")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "GLIPR2")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "LYZ"):which(Myeloid_genes_reorder == "RASSF5")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "HSP90AA1")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "HSP90AA1")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "HSPA1A"):which(Myeloid_genes_reorder == "HSP90AA1")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "LAP3")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "LAP3")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "EEF1A1"):which(Myeloid_genes_reorder == "CTSZ")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "LAP3")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "LAP3")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "ASAP1"):which(Myeloid_genes_reorder == "ASGR2")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "GLIPR2")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "GLIPR2")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "ALOX15"):which(Myeloid_genes_reorder == "LGALS1")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "JAML")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "JAML")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "MMP12"):which(Myeloid_genes_reorder == "MIF")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "TNFRSF1B")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "TNFRSF1B")+1):length(Myeloid_genes_reorder)])
                                  #Myeloid_gene_breaks[which(Myeloid_gene_breaks =="TNFRSF1B")] <- "MIF"
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "EEF1A1"):which(Myeloid_genes_reorder == "CDKN1A")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "HSP90AA1")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "HSP90AA1")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "IFITM2"):which(Myeloid_genes_reorder == "ALOX15")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(genes.to.swap, Myeloid_genes_reorder)
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "LGALS2"):which(Myeloid_genes_reorder == "IER5")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(genes.to.swap, Myeloid_genes_reorder)
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "G0S2"):which(Myeloid_genes_reorder == "C1QB")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "GLIPR2")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "GLIPR2")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "ANXA11"):which(Myeloid_genes_reorder == "ARHGDIB")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "JAML")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "JAML")+1):length(Myeloid_genes_reorder)])
                                  
                                  #genes.to.swap <- Myeloid_genes_reorder[which(Myeloid_genes_reorder == "CALHM6"):which(Myeloid_genes_reorder == "GLIPR2")]
                                  #Myeloid_genes_reorder <- Myeloid_genes_reorder[-(which(Myeloid_genes_reorder %in% genes.to.swap))]
                                  #Myeloid_genes_reorder <- c(Myeloid_genes_reorder[1:which(Myeloid_genes_reorder == "ALOX15")],genes.to.swap, Myeloid_genes_reorder[(which(Myeloid_genes_reorder == "ALOX15")+1):length(Myeloid_genes_reorder)])
                                  #Myeloid_gene_breaks[which(Myeloid_gene_breaks =="GLIPR2")] <- "HSP90AA1"
                                  

                                  
                          plot_heatmap_srt(Myeloid_data_subset, 
                                           genes = Myeloid_genes_reorder, 
                                           gene_breaks = Myeloid_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          
                          Myeloid.plot.all.labels <- last_plot()
                          
                          plot_heatmap_srt(Myeloid_data_subset, 
                                           genes = c(Myeloid_genes_reorder, genes.to.not.plot),
                                           gene_breaks = c(Myeloid_gene_breaks,Myeloid_genes_reorder[length(Myeloid_genes_reorder)]),
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          Myeloid.plot.all.labels.full <- last_plot()

                          plot_heatmap_srt(Myeloid_data_subset, 
                                           genes = Myeloid_genes_reorder, 
                                           gene_breaks = Myeloid_gene_breaks,
                                           gene_names = FALSE, 
                                           row_enrichGO = TRUE,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.5,
                                           type = "bulk", 
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = FALSE)
                          
                          Myeloid.plot <- last_plot()
                          
              #LC
                          #LC_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "LC" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          LC_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "LC" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
                          
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
                                                              k = round(length(LC_heatmap_genes)^0.5, digits =0), 
                                                              show_k = T,
                                                              color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                              log_scale_base = 10,
                                                              ceiling = F,
                                                              floor = F)
                          
                          #Gene order:
                          #LC_genes <- c(2,    6,3,1,4)
                          LC_genes <- c(7,       6,8,4,2)
                          
                          LC_genes <- rev(LC_genes)
                          LC_genes_reorder <- c()
                          
                          #groups.to.not.plot <- c(5)
                          groups.to.not.plot <- c(1,3,5)
                          
                          
                          for (i in 1:length(LC_genes)) {
                            int1 <- LC_genes[i]
                            ind <- grep(paste0("^", int1, "$"), LC.heatmap[[2]]$cluster)
                            reorder <- names(LC.heatmap[[2]]$cluster)[ind]
                            LC_genes_reorder <- c(LC_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #LC_cluster_breaks <- c(6)
                          LC_cluster_breaks <- c(6)
                          
                          LC_cluster_breaks <- rev(LC_cluster_breaks)
                          
                          
                          LC_gene_breaks <- find_gene_breaks(LC.heatmap, LC_genes_reorder,LC_cluster_breaks)
                         
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), LC.heatmap[[2]]$cluster)
                            reorder <- names(LC.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }         
                          
                                              #Manuel editing:
                          
                                              #genes.to.swap <- LC_genes_reorder[(which(LC_genes_reorder == "TTC28.AS1"):which(LC_genes_reorder == "GINS2"))]
                                              #LC_genes_reorder <- LC_genes_reorder[-which(LC_genes_reorder %in% genes.to.swap)]
                                              #LC_genes_reorder <- c(LC_genes_reorder[1:(which(LC_genes_reorder == "SRGN"))], genes.to.swap, LC_genes_reorder[(which(LC_genes_reorder == "SRGN")+1):length(LC_genes_reorder)])
                                              #LC_gene_breaks[which(LC_gene_breaks == "SRGN")] <- "GINS2"
                                              
                                              #genes.to.swap <- LC_genes_reorder[(which(LC_genes_reorder == "SH3BGRL3"):which(LC_genes_reorder == "C15orf48"))]
                                              #LC_genes_reorder <- LC_genes_reorder[-which(LC_genes_reorder %in% genes.to.swap)]
                                              #LC_genes_reorder <- c(genes.to.swap, LC_genes_reorder)
                                              
                          
                                      #New Manuel editing:
                                          #ACTB and PFN1 moved to common response portion nex tto CDKN1A:
                                            genes.to.swap <- LC_genes_reorder[(which(LC_genes_reorder == "PFN1"):which(LC_genes_reorder == "ACTB"))]
                                            LC_genes_reorder <- LC_genes_reorder[-which(LC_genes_reorder %in% genes.to.swap)]
                                            LC_genes_reorder <- c(LC_genes_reorder[1:(which(LC_genes_reorder == "CDKN1A"))], genes.to.swap, LC_genes_reorder[(which(LC_genes_reorder == "CDKN1A")+1):length(LC_genes_reorder)])
                          

                          
                          plot_heatmap_srt(LC_data_subset, 
                                           genes = LC_genes_reorder, 
                                           gene_breaks = LC_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          LC.plot.all.labels <- last_plot()
                          
                          plot_violin_bar_srt(DE.input, gene = "CXCR5", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE)
                          
                          
                          plot_heatmap_srt(LC_data_subset, 
                                           genes = c(LC_genes_reorder, genes.to.not.plot),
                                           gene_breaks = c(LC_gene_breaks,LC_genes_reorder[length(LC_genes_reorder)]),
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          LC.plot.all.labels.full <- last_plot()
                          
                          plot_heatmap_srt(LC_data_subset, 
                                           genes = LC_genes_reorder, 
                                           gene_breaks = LC_gene_breaks,
                                           gene_names = FALSE, 
                                           row_enrichGO = TRUE,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.5,
                                           type = "bulk", 
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = FALSE)
                          
                          LC.plot <- last_plot()
                          
                          
              #NK  
                          #NK_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "NK" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          NK_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "NK" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
                          
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
                                                         k = round(length(NK_heatmap_genes)^0.5, digits =0), 
                                                         show_k = T,
                                                         color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                         log_scale_base = 10,
                                                         ceiling = F,
                                                         floor = F)
                          

                          #Gene order:
                          #NK_genes <- c(1,   7,5,8,2,    4,  3)
                          NK_genes <- c(9,12,2,       1,     8,3,4)
                          
                          NK_genes <- rev(NK_genes)
                          NK_genes_reorder <- c()
                          
                          #groups.to.not.plot <- c(6)
                          groups.to.not.plot <- c(5,6,7,10,11)
                          
                          for (i in 1:length(NK_genes)) {
                            int1 <- NK_genes[i]
                            ind <- grep(paste0("^", int1, "$"), NK.heatmap[[2]]$cluster)
                            reorder <- names(NK.heatmap[[2]]$cluster)[ind]
                            NK_genes_reorder <- c(NK_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #NK_cluster_breaks <- c(4)
                          NK_cluster_breaks <- c(1,8)
                          
                          NK_cluster_breaks <- rev(NK_cluster_breaks)
                          
                          
                          NK_gene_breaks <- find_gene_breaks(NK.heatmap, NK_genes_reorder,NK_cluster_breaks)
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), NK.heatmap[[2]]$cluster)
                            reorder <- names(NK.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }         
                          
                          #Manuel editing:

                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "PFN1"):which(NK_genes_reorder == "RSF1"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder, genes.to.swap)
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "CREM"):which(NK_genes_reorder == "TAGAP"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder[1:(which(NK_genes_reorder == "TYROBP"))], genes.to.swap, NK_genes_reorder[(which(NK_genes_reorder == "TYROBP")+1):length(NK_genes_reorder)])
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "GEM"):which(NK_genes_reorder == "CAMK2D"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder, genes.to.swap)
                                    #NK_gene_breaks <- c(NK_gene_breaks,"RSF1")
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "KRT14"):which(NK_genes_reorder == "VSTM2L"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder, genes.to.swap)
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "ITM2A"):which(NK_genes_reorder == "CDKN1B"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder, genes.to.swap)
                                    #NK_gene_breaks[which(LC_gene_breaks == "SRGN")] <- "GINS2"
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "ITM2A"):which(NK_genes_reorder == "CDKN1B"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder[1:(which(NK_genes_reorder == "RSF1"))], genes.to.swap, NK_genes_reorder[(which(NK_genes_reorder == "RSF1")+1):length(NK_genes_reorder)])
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "LGALS3"):which(NK_genes_reorder == "C11orf58"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder[1:(which(NK_genes_reorder == "RSF1"))], genes.to.swap, NK_genes_reorder[(which(NK_genes_reorder == "RSF1")+1):length(NK_genes_reorder)])
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "RGS10"):which(NK_genes_reorder == "SAMSN1"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder[1:(which(NK_genes_reorder == "RSF1"))], genes.to.swap, NK_genes_reorder[(which(NK_genes_reorder == "RSF1")+1):length(NK_genes_reorder)])
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "NSA2"):which(NK_genes_reorder == "GLUL"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "AREG"):which(NK_genes_reorder == "CAT"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder[1:(which(NK_genes_reorder == "RSF1"))], genes.to.swap, NK_genes_reorder[(which(NK_genes_reorder == "RSF1")+1):length(NK_genes_reorder)])
                                    
                                    #genes.to.swap <- NK_genes_reorder[(which(NK_genes_reorder == "TXNIP"):which(NK_genes_reorder == "TAGAP"))]
                                    #NK_genes_reorder <- NK_genes_reorder[-which(NK_genes_reorder %in% genes.to.swap)]
                                    #NK_genes_reorder <- c(NK_genes_reorder[1:(which(NK_genes_reorder == "RSF1"))], genes.to.swap, NK_genes_reorder[(which(NK_genes_reorder == "RSF1")+1):length(NK_genes_reorder)])
                                    
                                    #NK_gene_breaks[which(LC_gene_breaks == "NHP2")] <- "MANF"
                          
                      
                          plot_heatmap_srt(NK_data_subset, 
                                           genes = NK_genes_reorder, 
                                           gene_breaks = NK_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          NK.plot.all.labels <- last_plot()
                          
                          
                          plot_heatmap_srt(NK_data_subset, 
                                           genes = c(NK_genes_reorder, genes.to.not.plot),
                                           gene_breaks = c(NK_gene_breaks,NK_genes_reorder[length(NK_genes_reorder)]),
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          NK.plot.all.labels.full <- last_plot()
                          
                          plot_heatmap_srt(NK_data_subset, 
                                           genes = NK_genes_reorder, 
                                           gene_breaks = NK_gene_breaks,
                                           gene_names = FALSE, 
                                           row_enrichGO = TRUE,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.5,
                                           type = "bulk", 
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = FALSE)
                          
                          NK.plot <- last_plot()
                          
                          
              #pDC
                          #pDC_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "pDC" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          pDC_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "pDC" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
                          
                          set.seed(100)
                          pDC.heatmap <- plot_heatmap_srt(pDC_data_subset, 
                                                         genes = pDC_heatmap_genes, 
                                                         type = "bulk", 
                                                         facet_by = "Semi_Detailed_Celltype",
                                                         scale_group = "Semi_Detailed_Celltype",
                                                         cluster_by = "row", 
                                                         pdf_format = "tile", 
                                                         scale_by = "row",
                                                         text_angle = 90,
                                                         cluster_type = "kmeans", 
                                                         k = round(length(pDC_heatmap_genes)^0.5, digits =0), 
                                                         show_k = T,
                                                         color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                         log_scale_base = 10,
                                                         ceiling = F,
                                                         floor = F)
                          
                          #Gene order:
                          #pDC_genes <- c(4,  5,     6,     3,   2,1)
                          pDC_genes <- c(4,23,          27,19,6,20,21,22,8,9,14,16,17,          3,5,7,11,12,28,10)
                          
                          pDC_genes <- rev(pDC_genes)
                          pDC_genes_reorder <- c()
                          
                          #groups.to.not.plot <- c()
                          
                          groups.to.not.plot <- c(1,2,13,15,18,24,25,26)
                          
                          for (i in 1:length(pDC_genes)) {
                            int1 <- pDC_genes[i]
                            ind <- grep(paste0("^", int1, "$"), pDC.heatmap[[2]]$cluster)
                            reorder <- names(pDC.heatmap[[2]]$cluster)[ind]
                            pDC_genes_reorder <- c(pDC_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          #pDC_cluster_breaks <- c(5, 2)
                          pDC_cluster_breaks <- c(27,3)
                          
                          pDC_cluster_breaks <- rev(pDC_cluster_breaks)
                          
                          pDC_gene_breaks <- find_gene_breaks(pDC.heatmap, pDC_genes_reorder,pDC_cluster_breaks)
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), pDC.heatmap[[2]]$cluster)
                            reorder <- names(pDC.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }         
                             
                                 #Manuel editing:
        
                                        #genes.to.swap <- pDC_genes_reorder[(which(pDC_genes_reorder == "HLA.DMA"):which(pDC_genes_reorder == "SNRPB"))]
                                        #pDC_genes_reorder <- pDC_genes_reorder[-which(pDC_genes_reorder %in% genes.to.swap)]
                                        #pDC_genes_reorder <- c(pDC_genes_reorder[1:(which(pDC_genes_reorder == "FOXP1"))], genes.to.swap, pDC_genes_reorder[(which(pDC_genes_reorder == "FOXP1")+1):length(pDC_genes_reorder)])
        
                                        #genes.to.swap <- pDC_genes_reorder[(which(pDC_genes_reorder == "UQCRH"):which(pDC_genes_reorder == "GSN"))]
                                        #pDC_genes_reorder <- pDC_genes_reorder[-which(pDC_genes_reorder %in% genes.to.swap)]
                                        #pDC_genes_reorder <- c(pDC_genes_reorder[1:(which(pDC_genes_reorder == "SNRPB"))], genes.to.swap, pDC_genes_reorder[(which(pDC_genes_reorder == "SNRPB")+1):length(pDC_genes_reorder)])
        
                                        #genes.to.swap <- pDC_genes_reorder[(which(pDC_genes_reorder == "INO80E"):which(pDC_genes_reorder == "NAP1L1"))]
                                        #pDC_genes_reorder <- pDC_genes_reorder[-which(pDC_genes_reorder %in% genes.to.swap)]
                                        

                          plot_heatmap_srt(pDC_data_subset, 
                                           genes = pDC_genes_reorder, 
                                           gene_breaks = pDC_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          pDC.plot.all.labels <- last_plot()
                          
                          plot_heatmap_srt(pDC_data_subset, 
                                           genes = c(pDC_genes_reorder, genes.to.not.plot),
                                           gene_breaks = c(pDC_gene_breaks,pDC_genes_reorder[length(pDC_genes_reorder)]),
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                        
                          pDC.plot.all.labels.full <- last_plot()
                          
                          
                          plot_heatmap_srt(pDC_data_subset, 
                                           genes = pDC_genes_reorder, 
                                           gene_breaks = pDC_gene_breaks,
                                           gene_names = FALSE, 
                                           row_enrichGO = TRUE,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.5,
                                           type = "bulk", 
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = FALSE)
                          
                          pDC.plot <- last_plot()
                          
            #Treg
                          #Treg_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "Treg" & logFC > 1) %>% dplyr::pull(Gene) %>% unique()
                          Treg_heatmap_genes <- contact_derm_DE_filtered %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(CellType == "Treg" & logFC > 1 & Color_factor != "Not Significant") %>% dplyr::pull(Gene) %>% unique()
                          
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
                                                            k = round(length(Treg_heatmap_genes)^0.5, digits =0), 
                                                            show_k = T,
                                                            color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                            log_scale_base = 10,
                                                            ceiling = F,
                                                            floor = F)
                          
                          #Gene order:
                          #Treg_genes <- c(3,    5,4,    2)
                          Treg_genes <- c(11,        3,1,     12)
                          
                          Treg_genes <- rev(Treg_genes)
                          Treg_genes_reorder <- c()
                          
                          #groups.to.not.plot <- c(1,6)
                          groups.to.not.plot <- c(2,4,5,6,7,8,9,10)
                          
                          for (i in 1:length(Treg_genes)) {
                            int1 <- Treg_genes[i]
                            ind <- grep(paste0("^", int1, "$"), Treg.heatmap[[2]]$cluster)
                            reorder <- names(Treg.heatmap[[2]]$cluster)[ind]
                            Treg_genes_reorder <- c(Treg_genes_reorder, reorder)
                          }

                          #Heatmap Gene breaks:
                          #Treg_cluster_breaks <- c(5,2)
                          Treg_cluster_breaks <- c(3,12)
                          
                          Treg_cluster_breaks <- rev(Treg_cluster_breaks)
                          
                          Treg_gene_breaks <- find_gene_breaks(Treg.heatmap, Treg_genes_reorder,Treg_cluster_breaks)
                          
                          
                          #Get full list of DE genes that I filtered out, to include in supplemental:
                          genes.to.not.plot <- c()
                          for (i in 1:length(groups.to.not.plot)) {
                            int1 <- groups.to.not.plot[i]
                            ind <- grep(paste0("^", int1, "$"), Treg.heatmap[[2]]$cluster)
                            reorder <- names(Treg.heatmap[[2]]$cluster)[ind]
                            genes.to.not.plot <- c(genes.to.not.plot, reorder)
                          }         
                                  
                                      #Manual Editing:
                                            #Treg_genes_reorder <- Treg_genes_reorder[-which(Treg_genes_reorder %in% c("CXCL14", "KRT1"))]
                                            #Treg_genes_reorder <- c(Treg_genes_reorder, "CXCL14", "KRT1")
                                            
                                            #genes.to.swap <- Treg_genes_reorder[(which(Treg_genes_reorder == "TUBB"):which(Treg_genes_reorder == "RAC2"))]
                                            #Treg_genes_reorder <- Treg_genes_reorder[-which(Treg_genes_reorder %in% genes.to.swap)]
                                            #Treg_genes_reorder <- c(genes.to.swap, Treg_genes_reorder)
                                
                                      #New Manual Editing:
                                          #Filter EEF1G through TUBB which look like irritant drop-out genes:
                                            genes.to.swap <- Treg_genes_reorder[(which(Treg_genes_reorder == "EEF1G"):which(Treg_genes_reorder == "TUBB"))]
                                            Treg_genes_reorder <- Treg_genes_reorder[-which(Treg_genes_reorder %in% genes.to.swap)]
                                            genes.to.not.plot <- c(genes.to.not.plot,genes.to.swap)
                                            
                          plot_heatmap_srt(Treg_data_subset, 
                                           genes = Treg_genes_reorder, 
                                           gene_breaks = Treg_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          Treg.plot.all.labels <- last_plot()
                          
                          plot_heatmap_srt(Treg_data_subset, 
                                           genes = c(Treg_genes_reorder, genes.to.not.plot),
                                           gene_breaks = c(Treg_gene_breaks,Treg_genes_reorder[length(Treg_genes_reorder)]),
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          Treg.plot.all.labels.full <- last_plot()
                          
                          plot_heatmap_srt(Treg_data_subset, 
                                           genes = Treg_genes_reorder, 
                                           gene_breaks = Treg_gene_breaks,
                                           gene_names = FALSE, 
                                           row_enrichGO = TRUE,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.5,
                                           type = "bulk", 
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = FALSE)
                        
                          Treg.plot <- last_plot()
                          
                          #Save heatmaps:
                          
                          setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps")
                          
                          #ggsave(CD4.plot, width = 600, height = 5*length(CD4_genes_reorder),units = "mm",filename = "CD4_DE.pdf", limitsize = FALSE) 
                          #ggsave(CD8.plot, width = 600, height = 5*length(CD8_genes_reorder),units = "mm",filename = "CD8_DE.pdf", limitsize = FALSE) 
                          #ggsave(KRT.plot, width = 600, height = 5*length(KRT_genes_reorder),units = "mm", filename = "KRT_DE.pdf", limitsize = FALSE) 
                          #ggsave(MEL.plot, width = 600, height = 5*length(MEL_genes_reorder),units = "mm",filename = "MEL_DE.pdf", limitsize = FALSE) 
                          #ggsave(Myeloid.plot, width = 600, height = 5*length(Myeloid_genes_reorder),units = "mm",filename = "Myeloid_DE.pdf", limitsize = FALSE)        
                          #ggsave(LC.plot, width = 600, height = 5*length(LC_genes_reorder),units = "mm",filename = "LC_DE.pdf", limitsize = FALSE)    
                          #ggsave(NK.plot, width = 600, height = 5*length(NK_genes_reorder),units = "mm",filename = "NK_DE.pdf", limitsize = FALSE)   
                          #ggsave(pDC.plot, width = 600, height = 5*length(pDC_genes_reorder),units = "mm",filename = "pDC_DE.pdf", limitsize = FALSE)  
                          #ggsave(Treg.plot, width = 600, height = 5*length(Treg_genes_reorder),units = "mm", filename = "Treg_DE.pdf", limitsize = FALSE)  
                          
                          
                          #ggsave(CD4.plot.all.labels, width = 200, height = 5*length(CD4_genes_reorder),units = "mm",filename = "CD4_DE.all.labels.pdf", limitsize = FALSE) 
                          #ggsave(CD8.plot.all.labels, width = 200, height = 5*length(CD8_genes_reorder),units = "mm",filename = "CD8_DE.all.labels.pdf", limitsize = FALSE) 
                          #ggsave(KRT.plot.all.labels, width = 200, height = 5*length(KRT_genes_reorder),units = "mm", filename = "KRT_DE.all.labels.pdf", limitsize = FALSE) 
                          #ggsave(MEL.plot.all.labels, width = 200, height = 5*length(MEL_genes_reorder),units = "mm",filename = "MEL_DE.all.labels.pdf", limitsize = FALSE) 
                          #ggsave(Myeloid.plot.all.labels, width = 200, height = 5*length(Myeloid_genes_reorder),units = "mm",filename = "Myeloid_DE.all.labels.pdf", limitsize = FALSE)        
                          #ggsave(LC.plot.all.labels, width = 200, height = 5*length(LC_genes_reorder),units = "mm",filename = "LC_DE.all.labels.pdf", limitsize = FALSE)    
                          #ggsave(NK.plot.all.labels, width = 200, height = 5*length(NK_genes_reorder),units = "mm",filename = "NK_DE.all.labels.pdf", limitsize = FALSE)   
                          #ggsave(pDC.plot.all.labels, width = 200, height = 5*length(pDC_genes_reorder),units = "mm",filename = "pDC_DE.all.labels.pdf", limitsize = FALSE)  
                          #ggsave(Treg.plot.all.labels, width = 200, height = 5*length(Treg_genes_reorder),units = "mm", filename = "Treg_DE.all.labels.pdf", limitsize = FALSE)  
                          
                          #ggsave(CD4.plot.all.labels.full, width = 200, height = 5*length(CD4_genes_reorder),units = "mm",filename = "CD4_DE.all.labels.full.pdf", limitsize = FALSE) 
                          #ggsave(CD8.plot.all.labels.full, width = 200, height = 5*length(CD8_genes_reorder),units = "mm",filename = "CD8_DE.all.labels.full.pdf", limitsize = FALSE) 
                          #ggsave(KRT.plot.all.labels.full, width = 200, height = 5*length(KRT_genes_reorder),units = "mm", filename = "KRT_DE.all.labels.full.pdf", limitsize = FALSE) 
                          #ggsave(MEL.plot.all.labels.full, width = 200, height = 5*length(MEL_genes_reorder),units = "mm",filename = "MEL_DE.all.labels.full.pdf", limitsize = FALSE) 
                          #ggsave(Myeloid.plot.all.labels.full, width = 200, height = 5*length(Myeloid_genes_reorder),units = "mm",filename = "Myeloid_DE.all.labels.full.pdf", limitsize = FALSE)        
                          #ggsave(LC.plot.all.labels.full, width = 200, height = 5*length(LC_genes_reorder),units = "mm",filename = "LC_DE.all.labels.full.pdf", limitsize = FALSE)    
                          #ggsave(NK.plot.all.labels.full, width = 200, height = 5*length(NK_genes_reorder),units = "mm",filename = "NK_DE.all.labels.full.pdf", limitsize = FALSE)   
                          #ggsave(pDC.plot.all.labels.full, width = 200, height = 5*length(pDC_genes_reorder),units = "mm",filename = "pDC_DE.all.labels.full.pdf", limitsize = FALSE)  
                          #ggsave(Treg.plot.all.labels.full, width = 200, height = 5*length(Treg_genes_reorder),units = "mm", filename = "Treg_DE.all.labels.full.pdf", limitsize = FALSE)  
                          
                          
                          
      #DE Gene Data table:
                    
                  DE.gene.stats.df <- data.frame(matrix(ncol = 8, nrow =0))
                  DE.gene.stats.df[,1] <- as.character(DE.gene.stats.df[,1])
                  DE.gene.stats.df[,2:8] <- lapply(DE.gene.stats.df[,2:8], FUN = as.numeric)
                  colnames(DE.gene.stats.df) <- c("Cell_type", "Total_cells", "Total_genes", "Total_DE_genes", 
                                                  "Allergy_specific_fraction","Nonspecific_upregulated_fraction", "Homeostatic_downregulated_fraction", "Unclear_trend_fraction")
                  
                  for (c in 1:length(unique(contact_derm_DE_filtered$CellType))) {
                    
                    cell = unique(contact_derm_DE_filtered$CellType)[c]
                    
                    Cell.count <- DE.input@meta.data %>% filter(Semi_Detailed_Celltype == cell) %>% nrow() 
                    total.genes <- get(paste(cell, "_data_subset", sep = ""))[which(rowSums(get(paste(cell, "_data_subset", sep = ""))) > 0), ]%>% nrow()
                    De.genes <- contact_derm_DE_filtered %>% filter(CellType == cell) %>% filter(Comparison != "Acetone_Vehicle" & Reference != "Acetone") %>% filter(logFC > 1 & padj <= 0.01) %>% dplyr::pull(Gene) %>% unique()
                    Total.de.genes.count <- as.numeric(length(De.genes))
                    
                    cell.heatmap.genes <- get(paste(cell, "_genes_reorder", sep = ""))

                    cell.heatmap.breaks <- get(paste(cell, "_gene_breaks", sep = ""))
                    
                    if (length(cell.heatmap.breaks) == 2) {
                      
                      Homeostatic.downregulated.de.genes <- length(cell.heatmap.genes[c(c(which(cell.heatmap.genes == cell.heatmap.breaks[2])+1):length(cell.heatmap.genes))])
                      nonspecific.upregulated.de.genes <- length(cell.heatmap.genes[c(c(which(cell.heatmap.genes == cell.heatmap.breaks[1])+1):which(cell.heatmap.genes == cell.heatmap.breaks[2]))])
                      allergy.specific.upregulated.de.genes <- length(cell.heatmap.genes[1:which(cell.heatmap.genes == cell.heatmap.breaks[1])])
                      unclear.trend.de.genes <- Total.de.genes.count - Homeostatic.downregulated.de.genes - nonspecific.upregulated.de.genes - allergy.specific.upregulated.de.genes
                      
                    }
                    
                    if (length(cell.heatmap.breaks) == 1) {
                      if (cell != "LC") {
                        stop ("I wrote this just for LC gene pattern")
                      }
                      
                      Homeostatic.downregulated.de.genes <- length(cell.heatmap.genes[c(c(which(cell.heatmap.genes == cell.heatmap.breaks)+1):length(cell.heatmap.genes))])
                      nonspecific.upregulated.de.genes <- length(cell.heatmap.genes[c(1:which(cell.heatmap.genes == cell.heatmap.breaks))])
                      allergy.specific.upregulated.de.genes <- 0
                      unclear.trend.de.genes <- Total.de.genes.count - Homeostatic.downregulated.de.genes - nonspecific.upregulated.de.genes - allergy.specific.upregulated.de.genes
                      
                    }
                    
                    row.data <- data.frame(Cell_type = cell,
                                           Total_cells = Cell.count, 
                                           Total_genes = total.genes,
                                           Total_DE_genes = Total.de.genes.count ,
                                           Allergy_specific_fraction = round(allergy.specific.upregulated.de.genes/Total.de.genes.count, digits = 2),
                                           Nonspecific_upregulated_fraction = round(nonspecific.upregulated.de.genes/Total.de.genes.count, digits = 2),
                                           Homeostatic_downregulated_fraction = round(Homeostatic.downregulated.de.genes/Total.de.genes.count,digits = 2),
                                           Unclear_trend_fraction = round(unclear.trend.de.genes/Total.de.genes.count, digits = 2))
                    
                    DE.gene.stats.df <- bind_rows(DE.gene.stats.df, row.data)
                    
                    Homeostatic.downregulated.de.genes <-0
                    nonspecific.upregulated.de.genes <-0
                    allergy.specific.upregulated.de.genes <-0
                    unclear.trend.de.genes <-0
                  }
                  
                  colnames(DE.gene.stats.df) <- c("Cell type", "Total cells", "Total genes", "Total DE genes", "Allergy-specific DE fraction", "Nonspecific upregulated DE fraction", "Homeostatic downregulated DE fraction", "Unclear trend DE fraction")
                  
                  DE.gene.stats.df <- DE.gene.stats.df %>% arrange(desc(`Allergy-specific DE fraction`), desc(`Total cells`))
                  
                  cell.numbers.df <- DE.input@meta.data %>% group_by_at(c("Lesion", "Semi_Detailed_Celltype")) %>% summarize(cell.count = n())
                    

                  library(data.table)
                  
                  Lesion.cell.numbers <- pivot_wider(cell.numbers.df, values_from = "cell.count", names_from = "Lesion")
                  colnames(Lesion.cell.numbers) <- c("Cell type","Nonlesional cells", "Irritant cells", "Acetone cells", "Day 2 Allergy cells", "Day 4 Allergy cells")
                  Lesion.cell.numbers[is.na(Lesion.cell.numbers)]<-0
                  
                  cell.numbers.df.nonzero <- bind_rows(cell.numbers.df, data.frame(Lesion = "Acetone_Vehicle", Semi_Detailed_Celltype = "pDC", cell.count = 0))
                  cell.numbers.df.nonzero$Lesion <- factor(cell.numbers.df.nonzero$Lesion, 
                                                           levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
                  cell.numbers.df.nonzero <- cell.numbers.df.nonzero %>% arrange(desc(Lesion))
                  cell.numbers.df.nonzero <- cell.numbers.df.nonzero %>% mutate(cell.count = cell.count+1)
                  cell.numbers.df.nonzero <- data.table(cell.numbers.df.nonzero)
                  colnames(cell.numbers.df.nonzero)[which(colnames(cell.numbers.df.nonzero) == "Semi_Detailed_Celltype")] <- "Cell type"
                  
                  cell.numbers.df <- data.table(cell.numbers.df)
                  colnames(cell.numbers.df)[which(colnames(cell.numbers.df) == "Semi_Detailed_Celltype")] <- "Cell type"
                  
                  sparkline_html <- cell.numbers.df.nonzero[, .(`Lesion Cell distributions` = spk_chr(log10(cell.count), type = 'bar', colorMap = c("#399A4B", "#ED9D31", "#9D9D9D", "#CC4141", "#E10E0E"))), by = "Cell type"]
                  
                  
                  bg = function(start, end, color, ...) {
                    paste("linear-gradient(90deg,transparent ",percent(start),",",
                          color, percent(start), ",", color, percent(end),
                          ", transparent", percent(end),")")
                  } 
                  
                  color_bar2 =  function (color = "lightgray", fun = "proportion", ...) {
                    fun <- match.fun(fun)
                    formatter("span", style = function(x) style(display = "inline-block",
                                                                `unicode-bidi` = "plaintext", 
                                                                "background" = bg(1-fun(as.numeric(x), ...), 1, color), "width"="100%" ))
                  }
                  
                  DE.sparkplot <- merge(DE.gene.stats.df, sparkline_html, by = "Cell type") %>% merge(.,Lesion.cell.numbers)  %>%
                    relocate(c(`Cell type`,`Total cells`,`Lesion Cell distributions`, `Nonlesional cells`, `Irritant cells`, `Acetone cells`, `Day 2 Allergy cells`, `Day 4 Allergy cells`)) %>%
                    formattable::formattable(align = c("l","r", "c", rep("r", times = 11)),
                                             list(`Allergy-specific DE fraction` = color_bar2("#CB708F"),
                                                  `Total cells` = color_bar2("#A7B9C4", fun = function (x) log10(x)/max(abs(log10(x)))),
                                                  `Total genes` = color_bar2("#A7B9C4"),
                                                  `Total DE genes` = color_bar2("#A7B9C4"))) %>% 
                    arrange(desc(`Allergy-specific DE fraction`)) %>%
                    as.htmlwidget() %>% 
                    spk_add_deps()
                  
                  DE.sparkplot
                 
                  #Note: the sparkline barplots don't seem accurate. For publication, i recreated a regular cell count bar plot myself and copy-pasted it in.
                  
                  #setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps")
                  
                  #htmlwidgets::saveWidget(DE.sparkplot, "DE.sparkplot.html", selfcontained = TRUE)
                  
                  #webshot::webshot(url = "DE.sparkplot.html", file = "DE.sparkplot.jpeg",vwidth = 2000, vheight = 500)
              
                  #Proper Cell count Barplots:
                  total.cell.counts <- DE.input@meta.data %>% group_by(Semi_Detailed_Celltype, Lesion) %>% summarize(cellcount = n())
                  
                  Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
                  
                  ggplot(total.cell.counts, aes(x = Lesion, y = log2(cellcount), fill = Lesion, color = Lesion)) +
                    geom_col()+
                    facet_grid(Semi_Detailed_Celltype~., scales="free_y")+
                    theme_classic()+
                    scale_color_manual(values = Lesion.colors)+
                    scale_fill_manual(values = Lesion.colors)
                  
                  semi_detailed_log2cellcount<- last_plot()
                  
                  ggplot(total.cell.counts, aes(x = Lesion, y = cellcount, fill = Lesion, color = Lesion)) +
                    geom_col()+
                    facet_grid(Semi_Detailed_Celltype~., scales="free_y")+
                    theme_classic()+
                    scale_color_manual(values = Lesion.colors)+
                    scale_fill_manual(values = Lesion.colors)
                  
                  semi_detailed_cellcount<- last_plot()
                  
                  #ggsave(semi_detailed_log2cellcount, width = 25, height = 25,units = "cm",filename = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps/semi_detailed_log2cellcount.pdf", limitsize = FALSE) 
                  
                  
  ## Allergy-specific DE Gene GO ontology bar plots:                  
          library(clusterProfiler)
          library(org.Hs.eg.db)
                
                  enrichgo.result.string <- c()
                  
                  for (c in 1:length(unique(contact_derm_DE_filtered$CellType))) {
                    
                        cell = unique(contact_derm_DE_filtered$CellType)[c]
                        
                        cell.heatmap.genes <- get(paste(cell, "_genes_reorder", sep = ""))
                        cell.heatmap.breaks <- get(paste(cell, "_gene_breaks", sep = ""))
                        
                        if (length(cell.heatmap.breaks) == 2) {
                          
                          allergy.specific.upregulated.de.genes <- cell.heatmap.genes[1:which(cell.heatmap.genes == cell.heatmap.breaks[1])]
                          
                        }
                        
                        if (length(cell.heatmap.breaks) == 1) {
                          if (cell != "LC") {
                            stop ("I wrote this just for LC gene pattern")
                          }
                          
                          allergy.specific.upregulated.de.genes <- c()
                          
                        }
                        
                        if ( length(allergy.specific.upregulated.de.genes) > 0) {
                          
                          universe_id <- bitr(allergy.specific.upregulated.de.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                          universe_id$ENTREZID <- as.numeric(universe_id$ENTREZID)
                          
                          ego <- enrichGO(gene          = universe_id$ENTREZID,
                                          universe      = names(geneList),
                                          OrgDb         = org.Hs.eg.db,
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.01,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)
                          
                          
                          enrich.go.result <- ego@result
                          enrich.go.result$Celltype <- cell
                          assign(paste(cell, "_enrichGO_result", sep = ""),enrich.go.result)
                          
                          enrichgo.result.string <- c(enrichgo.result.string, paste(cell, "_enrichGO_result", sep = ""))
                        }
                    

                  }
                  
                  enrich.go.result <- purrr::reduce(mget(enrichgo.result.string), bind_rows)                  
                
                  
                  enrich.go.result$GeneRatioNumeric <- sapply(enrich.go.result[,"GeneRatio"], FUN = function(x) {
                    gene.ratio <- unlist(str_split(x, pattern = "/"))
                    gene.ratio.numeric <- as.numeric(gene.ratio[1])/as.numeric(gene.ratio[2])
                    return(gene.ratio.numeric)
                  })
                  enrich.go.result$wt <- enrich.go.result$GeneRatioNumeric*enrich.go.result$Count
                  
                  # Consolidate for each cell type so that each combination of genes is only assigned to a single GO term:
                  enrich.go.result.filtered <- enrich.go.result[0,]
                  for (c in 1:length(unique(enrich.go.result$Celltype))) {
                    cell <- unique(enrich.go.result$Celltype)[c]
                    
                    cell.data <- enrich.go.result %>% filter(Celltype == cell)
                    
                    cell.data <- cell.data[-which(duplicated(cell.data$geneID)),]
                    
                    enrich.go.result.filtered <- bind_rows(enrich.go.result.filtered, cell.data)
                  }                
                  
                  
                  #Plot top results for each cell type:
                  n.to.plot = 5
                  
                          enrich.go.results.to.plot <- enrich.go.result.filtered %>% group_by(Celltype) %>% arrange(desc(GeneRatioNumeric)) %>% top_n(n = n.to.plot, wt = GeneRatioNumeric*Count)

                          enrich.go.results.to.plot2 <- enrich.go.results.to.plot[0,]
                          for (c in 1:length(unique(enrich.go.results.to.plot$Celltype))) {
                            cell <- unique(enrich.go.results.to.plot$Celltype)[c]
                            
                            rows <- filter(enrich.go.results.to.plot, Celltype == cell)
                            
                            enrich.go.results.to.plot2 <- bind_rows(enrich.go.results.to.plot2, rows[1:n.to.plot,])
                            
                          }
                          
                          enrich.go.results.to.plot2 <- enrich.go.results.to.plot2 %>% arrange(Celltype, desc(wt))
                  
                          for (c in 1:length(unique(enrich.go.results.to.plot2$Celltype))) {
                          
                            filtered.cell.GO.to.plot <- enrich.go.results.to.plot2 %>% filter(Celltype == unique(enrich.go.results.to.plot2$Celltype)[c]) 
                            filtered.cell.GO.to.plot$Description <- factor( filtered.cell.GO.to.plot$Description,
                                                                            levels = unique(filtered.cell.GO.to.plot$Description))
                            
                            plot <- ggplot(filtered.cell.GO.to.plot, aes(x = Description, y = Count, color = -log10(pvalue), fill = -log10(pvalue))) +
                                    geom_col() +
                                    theme_classic() +
                                    coord_flip() +
                                    scale_fill_gradientn(colors = colorRampPalette(c("#0D43F9","#F90D0D"))(256), limits=c(0, max(-log10(enrich.go.results.to.plot2$pvalue)))) +
                                    scale_color_gradientn(colors = colorRampPalette(c("#0D43F9", "#F90D0D"))(256),limits=c(0, max(-log10(enrich.go.results.to.plot2$pvalue)))) +
                                    labs(x = "Description", title = paste(unique(enrich.go.results.to.plot2$Celltype)[c])) +
                                    theme(plot.title = element_text(hjust = 0.5)) +
                                    scale_x_discrete(limits = rev(filtered.cell.GO.to.plot$Description))
                            
                            assign(paste(unique(enrich.go.results.to.plot2$Celltype)[c], "_enrichGO_barplot", sep = ""), plot)
                            
                            if (c == 1) {
                              plot.list <- list(get(paste(unique(enrich.go.results.to.plot2$Celltype)[c], "_enrichGO_barplot", sep = "")))
                            }
                            if (c > 1) {
                              plot.list[[c]] <- get(paste(unique(enrich.go.results.to.plot2$Celltype)[c], "_enrichGO_barplot", sep = ""))
                            }
                            
                            
                          }

                          enrichgo.finished.plot <- ggpubr::ggarrange(plotlist = plot.list, common.legend = TRUE)
                          
                          enrichgo.finished.plot
                          
                            #ggsave(enrichgo.finished.plot, width = 25, height = 15,units = "in",filename = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps/enrichgo/enrichgo.n5.plot.pdf", limitsize = FALSE)  
                            #ggsave(enrichgo.finished.plot, width = 25, height = 15,units = "in",filename = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps/enrichgo/enrichgo.n10.plot.pdf", limitsize = FALSE)  
           
  
## GSEA Plots:
        library(clusterProfiler)
        library(org.Hs.eg.db)
                          
                  ENTREZID_to_SYMBOL <- function(x) {
                            universe_id <- bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
                            
                            symbols <- universe_id$SYMBOL
                            
                            return(symbols)
                          }
                             
                                                      
                            enrichgo.result.string <- c()
                            
                            for (c in 1:length(unique(contact_derm_DE_filtered$CellType))) {
                              
                              cell = unique(contact_derm_DE_filtered$CellType)[c]
                              
                              cell.heatmap.genes <- get(paste(cell, "_genes_reorder", sep = ""))
                              cell.heatmap.breaks <- get(paste(cell, "_gene_breaks", sep = ""))
                              
                              if (length(cell.heatmap.breaks) == 2) {
                                
                                allergy.specific.upregulated.de.genes <- cell.heatmap.genes[1:which(cell.heatmap.genes == cell.heatmap.breaks[1])]
                                
                              }
                              
                              if (length(cell.heatmap.breaks) == 1) {
                                if (cell != "LC") {
                                  stop ("I wrote this just for LC gene pattern")
                                }
                                
                                allergy.specific.upregulated.de.genes <- c()
                                
                              }
                              
                              if ( length(allergy.specific.upregulated.de.genes) > 0) {
                                
                                universe_id <- bitr(allergy.specific.upregulated.de.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                                universe_id$ENTREZID <- as.numeric(universe_id$ENTREZID)
                                
                                ego <- enrichGO(gene          = universe_id$ENTREZID,
                                                universe      = names(geneList),
                                                OrgDb         = org.Hs.eg.db,
                                                ont           = "BP",
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.01,
                                                qvalueCutoff  = 0.05,
                                                readable      = TRUE)
                                
                                
                                enrich.go.result <- ego@result
                                enrich.go.result$Celltype <- cell
                                assign(paste(cell, "_enrichGO_result", sep = ""),enrich.go.result)
                                
                                enrichgo.result.string <- c(enrichgo.result.string, paste(cell, "_enrichGO_result", sep = ""))
                              }
                              
                              
                            }
                            
                            enrich.go.result <- purrr::reduce(mget(enrichgo.result.string), bind_rows)                  
                            
                            
                            enrich.go.result$GeneRatioNumeric <- sapply(enrich.go.result[,"GeneRatio"], FUN = function(x) {
                              gene.ratio <- unlist(str_split(x, pattern = "/"))
                              gene.ratio.numeric <- as.numeric(gene.ratio[1])/as.numeric(gene.ratio[2])
                              return(gene.ratio.numeric)
                            })
                            enrich.go.result$wt <- enrich.go.result$GeneRatioNumeric*enrich.go.result$Count
                            
                            # Consolidate for each cell type so that each combination of genes is only assigned to a single GO term:
                            enrich.go.result.filtered <- enrich.go.result[0,]
                            for (c in 1:length(unique(enrich.go.result$Celltype))) {
                              cell <- unique(enrich.go.result$Celltype)[c]
                              
                              cell.data <- enrich.go.result %>% filter(Celltype == cell)
                              
                              cell.data <- cell.data[-which(duplicated(cell.data$geneID)),]
                              
                              enrich.go.result.filtered <- bind_rows(enrich.go.result.filtered, cell.data)
                            }                
                            
                            
                            #Plot top results for each cell type:
                            n.to.plot = 2
                                    
                                  
                                  enrich.go.results.to.plot <- enrich.go.result.filtered %>% group_by(Celltype) %>% arrange(desc(GeneRatioNumeric)) %>% top_n(n = n.to.plot, wt = GeneRatioNumeric*Count)
                                  
                                  enrich.go.results.to.plot2 <- enrich.go.results.to.plot[0,]
                                  for (c in 1:length(unique(enrich.go.results.to.plot$Celltype))) {
                                    cell <- unique(enrich.go.results.to.plot$Celltype)[c]
                                    
                                    rows <- filter(enrich.go.results.to.plot, Celltype == cell)
                                    
                                    enrich.go.results.to.plot2 <- bind_rows(enrich.go.results.to.plot2, rows[1:n.to.plot,])
                                    
                                  }
                                  
                                  enrich.go.results.to.plot2 <- enrich.go.results.to.plot2 %>% arrange(Celltype, desc(wt))
                                  
                                  
                                  
                            enrich.go.results.to.plot2$Celltype <- factor(enrich.go.results.to.plot2$Celltype, 
                                                                            levels = levels(DE.input@meta.data$Semi_Detailed_Celltype))
                                  
                    #Custom GSEA using above GO terms and Allergy vs Irritant Fold Enrichment:
                       
                       #GO IDs and Gene lists:
                       GO.IDs <- unique(enrich.go.results.to.plot2$ID)
                       
                       GOID_to_entrezID_list <- as.list(org.Hs.egGO2ALLEGS)
                       
                       GOID_to_Gene <- data.frame(matrix(nrow = 0, ncol = 2))
                       colnames(GOID_to_Gene) <- c("GOID", "Genes")
                       class(GOID_to_Gene$GOID) <- "character"
                       class(GOID_to_Gene$Genes) <- "character"
                       
                       for (element in 1:length(GO.IDs)) {
                         
                         GOID <- GO.IDs[element]
                         
                         ENTREZIDs <- unique(as.numeric(unlist(GOID_to_entrezID_list[[GO.IDs[element]]])))
                         
                         Genes <- ENTREZID_to_SYMBOL(ENTREZIDs)
                         
                         element.data <- data.frame(GOID = GOID,
                                                    Genes = str_c(Genes, collapse = ", "))
                         
                         
                         GOID_to_Gene <- bind_rows(GOID_to_Gene, element.data)
                         
                       }
                       
                       GOID_to_Gene$Description <- NA
                       for (r in 1:nrow(GOID_to_Gene)) {
                         id <- GOID_to_Gene$GOID[r]
                         
                         description <- enrich.go.results.to.plot2 %>% filter(ID == id) %>% dplyr::pull(Description) %>% unique() %>% as.character()
                         description <- description[1]
                         
                         GOID_to_Gene$Description[r] <- description
                         
                       }
                       
                       #Convert to term2gene dataframe:

                       for (r in 1:nrow(GOID_to_Gene)) {
                         id <- GOID_to_Gene$GOID[r]
                         
                         symbols <- unlist(str_split(GOID_to_Gene$Genes[r], pattern = ", "))
                         symbols_id <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                         symbols_id$ENTREZID <- as.numeric(symbols_id$ENTREZID)
                         symbols_id$GOID <- id
                         
                         rows.to.add <- right_join(GOID_to_Gene, symbols_id)
                        
                         if (r == 1) {
                           term2gene <- rows.to.add
                         }
                         if (r > 1) {
                           term2gene <- bind_rows(term2gene,rows.to.add)
                         }
                         
                       }
                       
                       term2gene <- term2gene[,-which(colnames(term2gene) == "Genes")]
                       
                      
                       #Run Custom GSEA for each cell type:
                       for (c in 1:length(unique(contact_derm_DE_filtered$CellType))) {
                         
                         cell = unique(contact_derm_DE_filtered$CellType)[c]
                         
                         cell.irr.all.DE <- contact_derm_DE_filtered %>% filter(Reference == "Irritant" & Comparison == "Day2_Allergy" & CellType == cell) %>% arrange(desc(logFC))
                         
                         if ( nrow(cell.irr.all.DE) > 0) {
                           
                           #GSEA set using Allergy vs Irritant LogFC:
                               
                               temp.genelist <- cell.irr.all.DE$logFC
                               names(temp.genelist) <- cell.irr.all.DE$Gene
                               
                               universe_id <- bitr(names(temp.genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                               universe_id$ENTREZID <- as.numeric(universe_id$ENTREZID)
                               universe_id$logFC <- temp.genelist[match(universe_id$SYMBOL, names(temp.genelist))]
                               
                               gsea_set <- as.vector(universe_id$logFC)
                               names(gsea_set) <- universe_id$ENTREZID
                               
                               gsea_set <- sort(gsea_set, decreasing = TRUE)
                           
                          #Custom GSEA:
                          
                               Term2gene.df <- data.frame(term = term2gene$Description,
                                                          gene = term2gene$ENTREZID)
                               
                               temp.gsea <- GSEA(geneList = gsea_set,
                                                 TERM2GENE = Term2gene.df,
                                                 pvalueCutoff = 1,
                                                 by = 'fgsea',
                                                 nPerm = 1000,
                                                 minGSSize = 0,
                                                 maxGSSize = 25000)
                               
                               temp.gsea <- temp.gsea@result
                               temp.gsea$Celltype <- cell
                               
                               if (c ==1) {
                                 custom.gsea <- temp.gsea
                               }
                               if (c > 1) {
                                 custom.gsea <- bind_rows(custom.gsea,temp.gsea)
                                 
                               }
                         }
                      }    
                       
                             #Filter to High Enrichment GO Terms:
                             GO.term.max.NES <- custom.gsea %>% group_by(ID) %>% summarize(max(abs(NES)))
                             colnames(GO.term.max.NES)[2] <- "max_abs_NES"
                             
                             GO.terms.to.plot <- GO.term.max.NES %>% filter(max_abs_NES > 1.6) %>% dplyr::pull(ID)
                             
                             custom.gsea.to.plot <- custom.gsea %>% filter(ID %in% GO.terms.to.plot)

                             
                             #Cluster GO Terms:
                             gsea.results.spread <- custom.gsea.to.plot %>% dplyr::select(Description, Celltype, NES) %>%as.data.frame() %>% pivot_wider(names_from = Celltype, values_from = NES)
                             
                             gsea.results.spread <- column_to_rownames(gsea.results.spread, var = "Description")
                             
                             
                             gsea.results.spread.scaled <- scale(gsea.results.spread)
                             
                             
                             ord2 <- hclust( dist(gsea.results.spread.scaled, method = "euclidean"), method = "ward.D" )$order
                             
                             custom.gsea.to.plot$Description <- factor(custom.gsea.to.plot$Description,
                                                               levels = rownames(gsea.results.spread.scaled)[ord2])
                             
                             
                             #Cluster Cell Types:
                             ord3 <- hclust( dist(t(gsea.results.spread.scaled), method = "euclidean"), method = "ward.D" )$order
                             
                             custom.gsea.to.plot$Celltype <- factor(custom.gsea.to.plot$Celltype, 
                                                            levels = levels(DE.input@meta.data$Semi_Detailed_Celltype))
                       
                       
                       ggplot(custom.gsea.to.plot, aes(x = Celltype, y = Description, fill = NES))+
                                           geom_tile(color = "#676767") +
                                           theme_classic() +
                                           scale_fill_gradientn(colors = colorRampPalette(c("#0059FF","#FFFFFF","#FFFFFF","#FFFFFF","#F90D0D"))(256), limits= c((-1*max(GO.term.max.NES$max_abs_NES)), max(GO.term.max.NES$max_abs_NES))) +
                                           scale_color_gradientn(colors = colorRampPalette(c("#0059FF","#FFFFFF","#FFFFFF","#FFFFFF","#F90D0D"))(256),limits=c((-1*max(GO.term.max.NES$max_abs_NES)), max(GO.term.max.NES$max_abs_NES)))
                                         
                        GO.term.GSEA <- last_plot()             
                       
                        setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps")
                        
                        #ggsave(GO.term.GSEA, width = 30, height = 25, units = "cm",filename = "GO_aller_vs_irr_GSEA.pdf", limitsize = FALSE) 
                       
                       
                       
                       
                       
                       
                       
                       