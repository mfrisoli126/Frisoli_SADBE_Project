##################################################################
#Run CellphoneDB in Python before processing it here with nichenet:
##################################################################
library(Seurat)
library(tibble)
library(dplyr)
library(umap)
library(FNN)
library(igraph)
library(stringr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(purrr)
library(data.table)
library(readxl)
library(aliases2entrez)
library(clusterProfiler)
library(tidyverse)


setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

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
                
    # Aggregate count data:
    DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Lesion", "Semi_Detailed_Celltype"))
    
    #Visualize Sample subset:
    plot_tsne_metadata_srt(DE.input, color_by = "Semi_Detailed_Celltype", size = 0.25)
    
    #Aggregate Network:
    network.file.path = "Data/RL_network2/out/"
    
    deconvoluted <- fread(paste(network.file.path, "deconvoluted.txt", sep = ""))
    significant_means <- fread(paste(network.file.path, "significant_means.txt", sep = ""))
    pvals <- fread(paste(network.file.path, "pvalues.txt", sep = ""))
    means <- fread(paste(network.file.path, "means.txt", sep = ""))
    
    means <- as.data.frame(means)
    interactions <- colnames(means)[-c(1:11)]
    interactions.df <- str_split(interactions, "\\|") %>% as.data.frame() %>% t() %>% as.data.frame()
    rownames(interactions.df) <- NULL
    
    interactions.df <- cbind(interactions, interactions.df)
    interactions.df <- interactions.df %>% mutate(Lesion1 = sapply(strsplit(V1,"__"), getElement, 1),
                                                  Lesion2 = sapply(strsplit(V2,"__"), getElement, 1))
    
    interactions.df <- interactions.df %>% filter(Lesion1 == Lesion2)
    
    #Filter to interactions between celltypes of the same lesion:
    CCI_network <- means[,which(colnames(means) %in% c("interacting_pair","partner_a", "partner_b", "receptor_a", "receptor_b", "is_integrin", "annotation_strategy",interactions.df$interactions))]
    
    #Remove duplicated interactions containing duplicated data:
    CCI_network <- CCI_network[-which(duplicated(CCI_network$interacting_pair)),]
    
    #Correct mis-labeled receptors:
    CCI_network[which(str_detect(CCI_network$interacting_pair, pattern = "_CXCL16")),"receptor_b"] <- FALSE
    CCI_network[which(str_detect(CCI_network$interacting_pair, pattern = "CD6_")),"receptor_a"] <- TRUE
    CCI_network[which(str_detect(CCI_network$interacting_pair, pattern = "CADM1_")),"receptor_a"] <- TRUE
    CCI_network[which(str_detect(CCI_network$interacting_pair, pattern = "_CADM1")),"receptor_b"] <- TRUE
    CCI_network[which(str_detect(CCI_network$interacting_pair, pattern = "VCAM1_")),"receptor_a"] <- TRUE
    CCI_network[which(str_detect(CCI_network$interacting_pair, pattern = "PDCD1_")),"receptor_a"] <- FALSE
    
    #Filter to true ligand:receptor interactions, excluding ligand-ligand interactions and receptor:receptor interactions (which mostly seem to be adhesion related):
    receptor.receptor.interactions.to.filter <- CCI_network %>% filter(receptor_a == TRUE & receptor_b == TRUE) %>% dplyr::pull(interacting_pair) %>% unique()
    lig.lig.adhesion.interactions.to.filter <- CCI_network %>% filter(receptor_a == FALSE & receptor_b == FALSE) %>% dplyr::pull(interacting_pair) %>% unique()
    
    CCI_network <- CCI_network %>% filter(!interacting_pair %in% c(receptor.receptor.interactions.to.filter, lig.lig.adhesion.interactions.to.filter))
  
    rownames(CCI_network) <- NULL
    CCI_network <- column_to_rownames(CCI_network, "interacting_pair")
    
    #Filter genes with standard dev = 0:
    interaction.sd <- apply(CCI_network[,-c(1:6)], MARGIN = 1, FUN = sd) %>% as.data.frame()
    colnames(interaction.sd) <- "sd"
    interactions.to.filter <- interaction.sd %>% filter(sd == 0) %>% rownames()
    
    CCI_network <- CCI_network[-which(rownames(CCI_network) %in% interactions.to.filter),]
    colnames <- colnames(CCI_network)
    
    #Scale by row:
    CCI_network[,-c(1:6)] <- as.data.frame(t(apply(CCI_network[,-c(1:6)], MARGIN = 1, FUN = function(x) scale(x, center = FALSE))))
    colnames(CCI_network) <- colnames
    
    #Convert to same structure as SignallingSingleCell Network Dataframe:
    CCI_network <- rownames_to_column(CCI_network, var = "interaction")
    
    #Clarify breaks in complex interaction names:
    CCI_network$interaction[which(sapply(CCI_network$interaction, FUN = function(x) as.numeric(length(unlist(str_split(x, pattern = "_"))))) == 2)] <- str_replace(CCI_network$interaction[which(sapply(CCI_network$interaction, FUN = function(x) as.numeric(length(unlist(str_split(x, pattern = "_"))))) == 2)], pattern = "_", replacement = "_._")
    
    for (row in 1:nrow(CCI_network)) {
       
      if (str_detect(CCI_network$interaction[row], pattern = "_._") == TRUE) {
        next
      }
      
      if (str_detect(CCI_network$partner_a[row], pattern = "complex") == TRUE) {
        
        complex.name <- str_replace(CCI_network$partner_a[row], pattern = "complex:", replacement = "")
      
        CCI_network$interaction[row] <- str_replace(CCI_network$interaction[row], pattern = complex.name, replacement = paste(complex.name, "_.", sep = ""))
          
        rm(complex.name)
        next
      }
      
      if (str_detect(CCI_network$partner_b[row], pattern = "complex") == TRUE) {
        
        complex.name <- str_replace(CCI_network$partner_b[row], pattern = "complex:", replacement = "")
        
        CCI_network$interaction[row] <- str_replace(CCI_network$interaction[row], pattern = complex.name, replacement = paste("._",complex.name, sep = ""))
        
        rm(complex.name)
        next
      }
        
    }
    
    CCI_network <- CCI_network %>% mutate(partner_a = sapply(strsplit(interaction,"_._"), getElement, 1),
                                          partner_b = sapply(strsplit(interaction,"_._"), getElement, 2))
    
    #Convert to format similar to signalling single cell network:
    CCI_network <- pivot_longer(CCI_network, cols = -c(1:7), names_to = "sample", values_to = "Connection_product") 
    
    CCI_network <- CCI_network %>% mutate(partner_a_source = sapply(strsplit(sample,"\\|"), getElement, 1),
                                          partner_b_source = sapply(strsplit(sample,"\\|"), getElement, 2))

    CCI_network <- CCI_network %>% mutate(Ligand_Source = case_when(receptor_a == FALSE ~ partner_a_source,
                                                                    receptor_a == TRUE ~ partner_b_source),
                                          Ligand = case_when(receptor_a == FALSE ~ partner_a,
                                                             receptor_a == TRUE ~ partner_b),
                                          
                                          Receptor_source = case_when(receptor_a == FALSE ~ partner_b_source,
                                                                      receptor_a == TRUE ~ partner_a_source),
                                          Receptor = case_when(receptor_a == FALSE ~ partner_b,
                                                               receptor_a == TRUE ~ partner_a))
    
    CCI_network <-CCI_network[,which(colnames(CCI_network) %in% c("Ligand", "Receptor", "Ligand_Source", "Receptor_source", "is_integrin","annotation_strategy","Connection_product"))]
    CCI_network <- CCI_network %>% filter(annotation_strategy == "curated")
    
    CCI_network <- CCI_network[,c("Ligand", "Receptor", "Ligand_Source", "Receptor_source", "Connection_product")]
    
    CCI_network <- CCI_network %>% filter(Connection_product > 0)
    CCI_network <- CCI_network %>% mutate(log10_Connection_product = log10(Connection_product))


#########################################################
# Infer Receptor ---> Ligand connections using Nichenet:
#########################################################
library(nichenetr)        
    
  #Load input references for Nichenet:
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
    
    lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
    
      # Mike Edit: Only KEGG connections seem legit:
      lr_network = lr_network %>% mutate(bonafide = database %in% c("kegg"))
      lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
      lr_network = lr_network %>% filter(bonafide == TRUE)
      
    weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
    weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
    
  #Set Comparisons:
    comparisons <- as.character(sort(unique(contact_derm@meta.data$Lesion)))
    comparisons <- expand.grid(comparisons, comparisons, stringsAsFactors = FALSE)
    colnames(comparisons) <- c("reference", "comparison")
    comparisons <- comparisons[c(6,11,16,21),]

    #Add any multi-lesion combo comparisons you want:
    comparisons <- bind_rows(comparisons, data.frame(reference = str_c(string = c("Irritant", "Day2_Allergy", "Day4_Allergy"), collapse = "_._"),
                                                     comparison = "Nonlesional"))
    
    #Load DE Table:
    contact_derm_DE_filtered <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE_filtered.txt")
    
    DE.table.input <- contact_derm_DE_filtered
    
  #Set seurat object for Nichenet:
    seurat_obj <- DE.input
  
    
                                        #Example code for Regular nichenet:
                                                    ## receiver = individual celltype:
                                                    receiver = "Myeloid"
                                                    expressed_genes_receiver = get_expressed_genes(receiver, seurat_obj, pct = 0.0001, assay_oi = "RNA")
                                                    
                                                    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

                                                    ## senders = all celltypes in data:
                                                    sender_celltypes = as.character(unique(Idents(seurat_obj)))
                                                    list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.0001, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
                                                    expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
                                                    
                                                    #Find DE genes within Receiver Celltype:
                                                    seurat_obj_receiver= subset(seurat_obj, idents = receiver)
                                                    seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Lesion"]])
                                                    
                                                    unique(seurat_obj_receiver$Lesion)
                                                    condition_oi = c("Day2_Allergy")
                                                    condition_reference = c("Nonlesional")
                                                    
                                                    #Using raw counts is best?
                                                    DE_table_receiver = FindMarkers(object = seurat_obj_receiver, slot = "counts", ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.001, test.use = "bimod") %>% rownames_to_column("gene")
                                                    
                                                    geneset_oi = DE_table_receiver %>% filter(abs(avg_log2FC) >= 0.5) %>% arrange(desc(avg_log2FC)) %>% pull(gene)
                                                    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
                                                    
                                                    #Define Potential Ligands for DE genes using previous reference networks:
                                                    ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
                                                    receptors = lr_network %>% dplyr::pull(receptor) %>% unique()
                                                    
                                                    expressed_ligands = intersect(ligands,expressed_genes_sender)
                                                    expressed_receptors = intersect(receptors,expressed_genes_receiver)
                                                    
                                                    potential_ligands = lr_network %>% filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>% dplyr::pull(ligand) %>% unique()
                                                    potential_ligands
                                                    
                                                    #Use Nichenet to estimate Ligand activities based on cell transcription profiles:
                                                    ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
                                                    
                                                    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
                                                    
                                                    #Visualize top ranked ligands:
                                                    best_upstream_ligands = ligand_activities %>% top_n(40, pearson) %>% arrange(-pearson) %>% dplyr::pull(test_ligand) %>% unique()
                                                    
                                                    DotPlot(seurat_obj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
                                                    
                                                    #Active Target Gene inference using receiver cell data:
                                                    active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
                                                    
                                                    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.90)
                                                    
                                                    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
                                                    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
                                                    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
                                                    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
                                                    
                                                    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
                                                    p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
                                                    p_ligand_target_network
                                            
                                                    
                                                    #Add LogFC data between conditions to this data:
                                                    # this uses a new nichenetr function - reinstall nichenetr if necessary!
                                                    get_lfc_celltype <- function (celltype_oi, seurat_obj, slot = "counts", condition_colname, condition_oi, 
                                                                                  condition_reference, celltype_col = "celltype", expression_pct = 0.1) {
                                                      requireNamespace("Seurat")
                                                      requireNamespace("dplyr")
                                                      if (!is.null(celltype_col)) {
                                                        seurat_obj_celltype = SetIdent(seurat_obj, value = seurat_obj[[celltype_col]])
                                                        seuratObj_sender = subset(seurat_obj_celltype, idents = celltype_oi)
                                                      }
                                                      else {
                                                        seuratObj_sender = subset(seurat_obj, idents = celltype_oi)
                                                      }
                                                      seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[[condition_colname]])
                                                      
                                                      conditions.exist <- seuratObj_sender[[condition_colname]] %>% unlist() %>% unique() %>% as.character()
                                                      condition_oi.exist <- condition_oi[which(condition_oi %in% conditions.exist)]
                                                      condition_reference.exist <- condition_reference[which(condition_reference %in% conditions.exist)] 
                                                      
                                                      cells.oi.num <- seuratObj_sender@meta.data %>% filter(eval(sym(condition_colname)) %in% condition_oi.exist) %>% nrow() %>% as.numeric()
                                                      cells.ref.num <- seuratObj_sender@meta.data %>% filter(eval(sym(condition_colname)) %in% condition_reference.exist) %>% nrow() %>% as.numeric()
                                                      
                                                      if (cells.oi.num < 3 | cells.ref.num <3) {
                                                        DE_table_sender = data.frame(gene = rownames(seuratObj_sender),
                                                                                     avg_log2FC = NA)
                                                        colnames(DE_table_sender)[2] = celltype_oi
                                                        return(DE_table_sender)
                                                      }
                                                      
                                                      if (cells.oi.num > 3 & cells.ref.num > 3) {
                                                        DE_table_sender = FindMarkers(object = seuratObj_sender, slot = slot,
                                                                                      ident.1 = condition_oi.exist, ident.2 = condition_reference.exist, 
                                                                                      min.pct = expression_pct, logfc.threshold = 0.05) %>% 
                                                          rownames_to_column("gene")
                                                        SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_sender)
                                                        if (SeuratV4 == TRUE) {
                                                          DE_table_sender = DE_table_sender %>% as_tibble() %>% 
                                                            dplyr::select(-p_val) %>% dplyr::select(gene, avg_log2FC)
                                                        }
                                                        else {
                                                          DE_table_sender = DE_table_sender %>% as_tibble() %>% 
                                                            dplyr::select(-p_val) %>% dplyr::elect(gene, avg_logFC)
                                                        }
                                                        colnames(DE_table_sender) = c("gene", celltype_oi)
                                                        return(DE_table_sender)
                                                      }
                                                    }
                                                    
                                                    DE_table_all = Idents(seurat_obj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seurat_obj, condition_colname = "Lesion", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.0001, celltype_col = "Semi_Detailed_Celltype") %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
                                                    DE_table_all[is.na(DE_table_all)] = 0
                                                    
                                                    
                                                    
                                                    if (any(rowSums(DE_table_all[,-1]) == 0)) {
                                                      DE_table_all <- DE_table_all[-which(rowSums(DE_table_all[,-1]) == 0),]
                                                      
                                                    }
                                                    
                                                    # Combine ligand activities with DE information
                                                    ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
                                                    ligand_activities_de[is.na(ligand_activities_de)] = 0
                                                    
                                                    # make LFC heatmap
                                                    lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
                                                    rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
                                                    
                                                    order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
                                                    vis_ligand_lfc = lfc_matrix[order_ligands,] %>% as.matrix()
                                                    
                                                    colnames(vis_ligand_lfc) = DE_table_all[,-1] %>% colnames() %>% make.names()
                                                    
                                                    p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
                                                    p_ligand_lfc
                                                    
                                                    
                                                    #Summary combined Plots:
                                                            # ligand activity heatmap
                                                            ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
                                                            
                                                            rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
                                                            colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
                                                            
                                                            vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
                                                            p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
                                                            # ligand expression Seurat dotplot
                                                            order_ligands_adapted = order_ligands
                                                            #order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
                                                            #order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
                                                            
                                                            rotated_dotplot = DotPlot(seurat_obj %>% subset(Semi_Detailed_Celltype %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
                                                                              figures_without_legend = cowplot::plot_grid(
                                                                                p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                                                                rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
                                                                                p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                                                                p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                                                                align = "hv",
                                                                                nrow = 1,
                                                                                rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))
                                                                              
                                                          legends = cowplot::plot_grid(
                                                                                ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
                                                                                ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
                                                                                ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
                                                                                ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
                                                                                nrow = 1,
                                                                                align = "h", rel_widths = c(1.5, 1, 1, 1))
                                                            
                                                            combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
                                                            combined_plot
                                                            
    # Analyze all celltypes and conditions with Nichenet:
      
            gene.frac.cutoff = 0.001
                #original was 0.01
            
            Receptor_ligand_inference_list <- c()
          
            for (r in 1:nrow(comparisons)) {

                comp = comparisons$comparison[r]
                ref = comparisons$reference[r]

              print(paste("Analyzing", ref, "vs", comp, sep = " "))
                
              if (str_detect(comparisons$reference[r], pattern = "_._") == TRUE) {
                
                ref = unlist(str_split(comparisons$reference[r], pattern = "_._"))
                
              }
                
              if (str_detect(comparisons$comparison[r], pattern = "_._") == TRUE) {
                  
                  comp = unlist(str_split(comparisons$comparison[r], pattern = "_._"))
                
              }
              
            Idents(seurat_obj) <- "Lesion"
            
            seurat_subset <- subset(seurat_obj, idents = c(comp, ref))
            
            celltypes.to.analyze <- unique(seurat_subset@meta.data$Semi_Detailed_Celltype)
            
            Idents(seurat_subset) <- "Semi_Detailed_Celltype"
          
                      for (c in 1:length(celltypes.to.analyze)) {
                        
                                ## receiver = individual celltype:
                                receiver = celltypes.to.analyze[c]
                                
                                print(paste(". . . . . . . . Analyzing", receiver, sep = " "))
                                
                                expressed_genes_receiver = get_expressed_genes(receiver, seurat_subset, pct = gene.frac.cutoff, assay_oi = "RNA")
                                
                                background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
                                
                                ## senders = all celltypes in data:
                                sender_celltypes = as.character(unique(Idents(seurat_subset)))
                                list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_subset, gene.frac.cutoff, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
                                expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
                                
                                #Find DE genes within Receiver Celltype:
                                seurat_obj_receiver= subset(seurat_subset, idents = receiver)
                                seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Lesion"]])
                                
                                subset.condition_oi = comp
                                subset.condition_ref = ref
                              
                                
                                if (all(subset.condition_oi %in% unique(seurat_obj_receiver$Lesion)) == FALSE & any(subset.condition_oi %in% unique(seurat_obj_receiver$Lesion)) == TRUE) {
                                  subset.condition_oi = subset.condition_oi[which(subset.condition_oi %in% unique(seurat_obj_receiver$Lesion))]
                                }
                                if (all(subset.condition_ref %in% unique(seurat_obj_receiver$Lesion)) == FALSE & any(subset.condition_ref %in% unique(seurat_obj_receiver$Lesion)) == TRUE) {
                                  subset.condition_ref = subset.condition_ref[which(subset.condition_ref %in% unique(seurat_obj_receiver$Lesion))]
                                }
                                
                                cells.oi.num <- seurat_obj_receiver@meta.data %>% filter(seurat_obj_receiver$Lesion %in% subset.condition_oi) %>% nrow() %>% as.numeric()
                                cells.ref.num <- seurat_obj_receiver@meta.data %>% filter(seurat_obj_receiver$Lesion %in% subset.condition_ref) %>% nrow() %>% as.numeric()
                                
                                if (cells.oi.num < 3 | cells.ref.num < 3) {
                                  active_ligand_target_links_df <- data.frame(ligand = NA,
                                                                              target = NA, 
                                                                              weight = NA,
                                                                              Receiver = receiver,
                                                                              Condition = comp)
                                  
                                  active_ligand_target_links_df[,1] <- as.character(active_ligand_target_links_df[,1])
                                  active_ligand_target_links_df[,2] <- as.character(active_ligand_target_links_df[,2])
                                  active_ligand_target_links_df[,3] <- as.numeric(active_ligand_target_links_df[,3])
                                }
                                
                                if (cells.oi.num >= 3 & cells.ref.num >= 3) {
                                  #Should we just load DE results from EdgeR DE results that incorporate batch-effect correction instead of recalculating here?
                                  if (length(ref) == 1 & length(comp) ==1) {
                                    DE_table_receiver = DE.table.input %>% filter(Reference == ref & Comparison == comp & CellType == receiver)
                                  }
                                  
                                  if (length(ref) > 1 | length(comp) > 1) {
                                    DE_table_receiver = DE.table.input %>% filter(Reference ==  str_c(ref, collapse = "/") & Comparison == comp & CellType == receiver)
                                  }
                                  
                                  geneset_oi = DE_table_receiver %>% filter(logFC >= 0.25 & padj <= 0.001) %>% arrange(desc(logFC)) %>% pull(Gene)
                                  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
                                  
                                  if (length(geneset_oi) == 0) {
                                    active_ligand_target_links_df <- data.frame(ligand = NA,
                                                                                target = NA, 
                                                                                weight = NA,
                                                                                Receiver = receiver,
                                                                                Condition = comp)
                                    
                                    active_ligand_target_links_df[,1] <- as.character(active_ligand_target_links_df[,1])
                                    active_ligand_target_links_df[,2] <- as.character(active_ligand_target_links_df[,2])
                                    active_ligand_target_links_df[,3] <- as.numeric(active_ligand_target_links_df[,3])
                        }
                                  
                                  if (length(geneset_oi) > 0) {
                                    
                                    #Define Potential Ligands for DE genes using previous reference networks:
                                    ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
                                    receptors = lr_network %>% dplyr::pull(receptor) %>% unique()
                                    
                                    expressed_ligands = intersect(ligands,expressed_genes_sender)
                                    expressed_receptors = intersect(receptors,expressed_genes_receiver)
                                    
                                    potential_ligands = lr_network %>% filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>% dplyr::pull(ligand) %>% unique()
                                    potential_ligands
                                    
                                    #Use Nichenet to estimate Ligand activities based on cell transcription profiles:
                                    ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
                                    
                                    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
                                    
                                    #Visualize top ranked ligands:
                                    
                                    #Should this be all ligands instead of top 40?
                                    #best_upstream_ligands = ligand_activities %>% top_n(40, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
                                    
                                    #Active Target Gene inference using receiver cell data:
                                    active_ligand_target_links_df = unique(ligand_activities$test_ligand) %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
                                    
                                    if (nrow(active_ligand_target_links_df) >0) {
                                      active_ligand_target_links_df$weight = as.numeric(active_ligand_target_links_df$weight)
                                      active_ligand_target_links_df$Receiver = receiver
                                      active_ligand_target_links_df$Condition = subset.condition_oi
                                    }
                                    if(nrow(active_ligand_target_links_df) == 0) {
                                      active_ligand_target_links_df <- data.frame(ligand = NA,
                                                                                  target = NA, 
                                                                                  weight = NA,
                                                                                  Receiver = receiver,
                                                                                  Condition = comp)
                                      
                                      active_ligand_target_links_df[,1] <- as.character(active_ligand_target_links_df[,1])
                                      active_ligand_target_links_df[,2] <- as.character(active_ligand_target_links_df[,2])
                                      active_ligand_target_links_df[,3] <- as.numeric(active_ligand_target_links_df[,3])
                                    } 
                                    
                                  }
                                }
                                
                                assign(paste(receiver,"_",comp, "_receptor_ligand_targets.df", sep = ""), active_ligand_target_links_df)     
                                
                                Receptor_ligand_inference_list <- c(Receptor_ligand_inference_list, paste(receiver,"_",comp, "_receptor_ligand_targets.df", sep = ""))           
            }
               
            
            }

          Receptor_ligand_inference_df <- purrr::reduce(mget(Receptor_ligand_inference_list), bind_rows)
         
          #Target Gene Fix:
                      gene.fix.df <- as.data.frame(unique(Receptor_ligand_inference_df$target))
                      colnames(gene.fix.df) = "Original_symbol"
                      
                      file <- system.file("extdata", "HGNC.txt", package = "aliases2entrez")
                      HGNC <- update_symbols()
                      symbols <- gene.fix.df$Original_symbol
                      ids <- convert_symbols(symbols, HGNC, c = 1)
                      
                      gene.fix.df <- cbind(gene.fix.df, ids)
                      names(gene.fix.df)[names(gene.fix.df) == "Symbols"] <- "A2E_Symbols"
                      
                      positive_ids <- ids[which(!is.na(ids$entrezID)),]
                      id.alias <- bitr(positive_ids$entrezID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
                      id.alias$ENTREZID <- as.numeric(id.alias$ENTREZID)
                      
                      gene.fix.df <- left_join(gene.fix.df, id.alias, by = c("entrezID" = "ENTREZID"))
                      names(gene.fix.df)[names(gene.fix.df) == "SYMBOL"] <- "Reformatted_symbol"
                      names(gene.fix.df)[names(gene.fix.df) == "entrezID"] <- "ENTREZID"
                      
                      entrez.to.ensembl <- bitr(gene.fix.df$ENTREZID, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")
                      entrez.to.ensembl$ENTREZID <- as.numeric(entrez.to.ensembl$ENTREZID) 
                      
                      gene.fix.df <- as.data.frame(full_join(gene.fix.df, entrez.to.ensembl))
                      gene.fix.df <- gene.fix.df[-which(duplicated(gene.fix.df$Original_symbol)),] 
                      
                      #Fix genes that need to be fixed here:
                      genes.to.fix <- gene.fix.df %>% filter(Original_symbol != Reformatted_symbol) %>% dplyr::select(Original_symbol) %>% unlist() %>% as.character()
                      
                      if (length(genes.to.fix) >0) {
                        
                          gene.fix.string <- c()
                          for (element in Receptor_ligand_inference_df$target[which(Receptor_ligand_inference_df$target %in% genes.to.fix)]) {
                            
                            row <- gene.fix.df[which(gene.fix.df$Original_symbol == element),]
                            gene.fix.string <- c(gene.fix.string, row$Reformatted_symbol)
                          }
                          
                          Receptor_ligand_inference_df$target[which(Receptor_ligand_inference_df$target %in% genes.to.fix)] <- gene.fix.string
                      }
                      
          
          
          Nichenet_receptor_induced_ligands_network <- Receptor_ligand_inference_df %>% filter(target %in% unique(c(CCI_network$Ligand, lr_network$ligand)))

          #Save each network:
          # save(CCI_network, file = "Data/RL_network2/cellphoneDB_network.Rdata")
          # save(Receptor_ligand_inference_df, file = "Data/RL_network2/nichenet_full_network_frac0.001.Rdata")
          # save(Nichenet_receptor_induced_ligands_network, file = "Data/RL_network2/nichenet_network_frac0.001.Rdata")

#########################################################
# Merge Networks:
#########################################################
          load("Data/RL_network2/cellphoneDB_network.Rdata")
          load("Data/RL_network2/nichenet_network_frac0.001.Rdata")

          #Reload network databases:
                network.file.path = "Data/RL_network2/out/"
                deconvoluted <- fread(paste(network.file.path, "deconvoluted.txt", sep = ""))
                means <- fread(paste(network.file.path, "means.txt", sep = ""))
                
                lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
                
                # Mike Edit: Only KEGG connections seem legit:
                lr_network = lr_network %>% mutate(bonafide = database %in% c("kegg"))
                lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
                lr_network = lr_network %>% filter(bonafide == TRUE)
                
                
# Correct both networks before merging:
          
          CCI_network <- CCI_network %>% mutate(interaction_pair = paste(Ligand, "_._", Receptor, sep = "")) 
          
          #Some ligand-receptor interactions in CellphoneDB network look incorrect: 
            #PECAM1 seems like it should only really homodimerize with PECAM1, but apparently this is curated?
            CCI_network <- CCI_network[-which(CCI_network$interaction_pair == "CD177_._PECAM1"),]

            #Everything about F11R (aka CD321) looks inconsistent, so i'm going to remove this from both networks:
            CCI_network <- CCI_network[-which(CCI_network$interaction_pair %in% c("BDNF_._F11R", "integrin_aLb2_complex_._F11R")),]
            
            if (any(Nichenet_receptor_induced_ligands_network$ligand == "F11R") == TRUE) {
              Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[-which(Nichenet_receptor_induced_ligands_network$ligand == "F11R"),]
            }
          
        
          #Both CellphoneDB database and Nichenet have some backwards receptor: ligand interactions that need to be corrected:
          ligands.not.in.CellphoneDB <- unique(Nichenet_receptor_induced_ligands_network$ligand)[-which(unique(Nichenet_receptor_induced_ligands_network$ligand) %in% CCI_network$Ligand)]
          
          #Some interactions are just mismatched backwards:
          mismatching_rls <- ligands.not.in.CellphoneDB[which(ligands.not.in.CellphoneDB %in% CCI_network$Receptor)]
          
              #Manually check here to see which order seems correct:
              true_receptors = c("CD40", "NECTIN1")
              true_ligands = mismatching_rls[-which(mismatching_rls %in% true_receptors)]
              
          #Correct Nichenet mislabeled Ligands:
          for (r in 1:length(true_receptors)) {
            
            backwards.gene <- true_receptors[r]
            
            backwards.data <- Nichenet_receptor_induced_ligands_network[which(Nichenet_receptor_induced_ligands_network$ligand == backwards.gene),]
            
            correct_ligand <- CCI_network %>% filter(Receptor == backwards.gene) %>% dplyr::select(Ligand) %>% unlist() %>% unique()
            
                  for (c in 1:length(correct_ligand)) {
                    ligand.data <- backwards.data
                    ligand.data$ligand <- correct_ligand[c]
                    
                    Nichenet_receptor_induced_ligands_network <- bind_rows(Nichenet_receptor_induced_ligands_network, ligand.data)
                  }
            
            Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[-which(Nichenet_receptor_induced_ligands_network$ligand == backwards.gene),]
            
            
          }
            
          
          #Correct CellphoneDB mislabeled Receptors:
          for (p in 1:length(true_ligands)) {
            
            backwards.gene <- true_ligands[p]
            
            backwards.rows <- CCI_network %>% filter(Receptor == backwards.gene)
            
            corrected.rows <- backwards.rows
            
            corrected.rows$Ligand <- backwards.rows$Receptor
            corrected.rows$Ligand_Source <- backwards.rows$Receptor_source
            corrected.rows$Receptor <- backwards.rows$Ligand
            corrected.rows$Receptor_source <- backwards.rows$Ligand_Source
            corrected.rows$interaction_pair <- paste(corrected.rows$Ligand, "_", corrected.rows$Receptor, sep = "")
            
            CCI_network[which(CCI_network$Receptor == backwards.gene),] <- corrected.rows
          }

          
          
          #Some genes are part of multi-gene complexes:        
          ligands.not.in.CellphoneDB <- unique(Nichenet_receptor_induced_ligands_network$ligand)[-which(unique(Nichenet_receptor_induced_ligands_network$ligand) %in% CCI_network$Ligand)]
          
          protein.subunit.genes <- ligands.not.in.CellphoneDB[which(ligands.not.in.CellphoneDB %in% deconvoluted$gene_name)]
          
          deconvoluted.subunits <- deconvoluted %>% filter(gene_name %in% protein.subunit.genes)
          deconvoluted.subunits <- deconvoluted.subunits[,c(1:5)]
          deconvoluted.subunits <- deconvoluted.subunits %>% filter(is_complex == TRUE)
          
          
          deconvoluted.subunits <- deconvoluted.subunits[-which(duplicated(deconvoluted.subunits$complex_name)),]
          
          deconvoluted.subunits$gene_complex <- paste(deconvoluted.subunits$gene_name, "_", deconvoluted.subunits$complex_name, sep = "")
          if (any(duplicated(deconvoluted.subunits$gene_complex))) {
            deconvoluted.subunits <- deconvoluted.subunits[-which(duplicated(deconvoluted.subunits$gene_complex)),]
          }

          #Convert Nichenet genes to CellphoneDB network protein complexes:
          #Ligand subunits:
          ligands <- unique(CCI_network$Ligand[which(CCI_network$Ligand %in% deconvoluted.subunits$complex_name)])
          
          if (length(ligands) > 0) {
            
            for (l in 1:length(unique(deconvoluted.subunits$gene_name))) {
              
              Nichenet.gene <- unique(deconvoluted.subunits$gene_name)[l]
              
              gene.complex <- deconvoluted.subunits %>% filter(gene_name == Nichenet.gene) %>% dplyr::select(complex_name) %>% unlist() %>% as.character() %>% unique()
              
              if (any(gene.complex %in% ligands)) {
                
                original.gene.data <- Nichenet_receptor_induced_ligands_network[which(Nichenet_receptor_induced_ligands_network$ligand == Nichenet.gene),]
                
                for (c in 1:length(gene.complex)) {
                  complex.data <- original.gene.data
                  complex.data$ligand <- gene.complex[c]
                  
                  Nichenet_receptor_induced_ligands_network <- bind_rows(Nichenet_receptor_induced_ligands_network, complex.data)
                }
                
                Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[-which(Nichenet_receptor_induced_ligands_network$ligand == Nichenet.gene),]
                
              }
            }
          }
          
          #Genes that are actually receptor subunits:
          receptor_complexes <- unique(CCI_network$Receptor[which(CCI_network$Receptor %in% deconvoluted.subunits$complex_name)])
          
          if (length(receptor_complexes) > 0) {
            
            #Correct mislabeled nichenet genes that are really mislabeled receptor subunits:
            for (l in 1:length(unique(deconvoluted.subunits$gene_name))) {
          
              Nichenet.gene <- unique(deconvoluted.subunits$gene_name)[l]
              
              gene.complex <- deconvoluted.subunits %>% filter(gene_name == Nichenet.gene) %>% dplyr::select(complex_name) %>% unlist() %>% as.character() %>% unique()
              
              if (any(gene.complex %in% receptor_complexes)) {
                
                #Correct nichenet data with proper ligands:
                backwards.data <- Nichenet_receptor_induced_ligands_network[which(Nichenet_receptor_induced_ligands_network$ligand == Nichenet.gene),]
                
                if (nrow(backwards.data)>0) {
                  
                    correct_ligand <- CCI_network %>% filter(Receptor %in% gene.complex) %>% dplyr::select(Ligand) %>% unlist() %>% unique()
                    
                    for (c in 1:length(correct_ligand)) {
                      ligand.data <- backwards.data
                      ligand.data$ligand <- correct_ligand[c]
                      
                      Nichenet_receptor_induced_ligands_network <- bind_rows(Nichenet_receptor_induced_ligands_network, ligand.data)
                }
                    
                    Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[-which(Nichenet_receptor_induced_ligands_network$ligand == Nichenet.gene),]
                }
              }
            }
          }
          
          #Some remaining Nichenet ligands look excluded from CellphoneDB database:
            ligands.not.in.CellphoneDB <- unique(Nichenet_receptor_induced_ligands_network$ligand)[-which(unique(Nichenet_receptor_induced_ligands_network$ligand) %in% c(CCI_network$Ligand, deconvoluted.subunits$complex_name))]
            
            ligs.not.in.cellphonedb.database <- unique(lr_network$ligand)[-which(unique(lr_network$ligand) %in% deconvoluted$gene_name)]
            ligs.not.in.cellphonedb.database <- intersect(ligands.not.in.CellphoneDB, ligs.not.in.cellphonedb.database)
            
            lr_interactions.not.in.CellphoneDB <- lr_network[which(lr_network$ligand %in% ligs.not.in.cellphonedb.database),]
              
            aggregate_bulk_data <- DE.input@assays$RNA@meta.features[,which(str_detect(colnames(DE.input@assays$RNA@meta.features), pattern = "bulk"))]
              
            lr_interactions.not.in.CellphoneDB <- lr_interactions.not.in.CellphoneDB %>% filter(ligand %in% rownames(aggregate_bulk_data) & receptor %in% rownames(aggregate_bulk_data))
            
            aggregate_bulk_data <- aggregate_bulk_data[which(rownames(aggregate_bulk_data) %in% c(lr_interactions.not.in.CellphoneDB$ligand, lr_interactions.not.in.CellphoneDB$receptor)),]
              
          #Manually create ligand-receptor network for these connections:
                aggregate_bulk_data <- t(aggregate_bulk_data) %>% as.data.frame() %>% rownames_to_column(var = "group") %>% pivot_longer(names_to = "Gene", values_to = "Count", cols = 2:ncol(.))
              
                aggregate_bulk_data$Lesion <- sapply(aggregate_bulk_data$group, FUN = function(x) unlist(str_split(x, pattern = "__"))[1])
                aggregate_bulk_data$Celltype <- sapply(aggregate_bulk_data$group, FUN = function(x) unlist(str_split(x, pattern = "__"))[2])
                
                aggregate_bulk_data <- aggregate_bulk_data[,-1]
                
                manual.network <- CCI_network[0,]
                
                for (r in 1:nrow(lr_interactions.not.in.CellphoneDB)) {
                  
                  lig <- lr_interactions.not.in.CellphoneDB$ligand[r]
                  rec <- lr_interactions.not.in.CellphoneDB$receptor[r]
                  
                  lr_data <- filter(aggregate_bulk_data, Gene %in% c(lig, rec))
                  
                    for (l in 1:length(unique(lr_data$Lesion))) {
                      
                      les <- unique(lr_data$Lesion)[l]
                    
                      les_network <- expand_grid(unique(lr_data$Celltype), unique(lr_data$Celltype)) %>% as.data.frame()
                      colnames(les_network) <- c("Ligand_Source", "Receptor_source")
                      
                      les_network$Ligand <- lig
                      les_network$Receptor <- rec
                      les_network$Lesion <- les
                      les_network$interaction_pair <- paste(les_network$Ligand, "_", les_network$Receptor, sep = "")
      
                      lr_les_data <- filter(lr_data, Lesion == les)
                    
                      lig_data <- lr_les_data %>% filter(Gene == lig) %>% rename("lig_expression" = Count, "Ligand" = Gene, "Ligand_Source" = Celltype)
                        
                      rec_data <- lr_les_data %>% filter(Gene == rec) %>% rename("rec_expression" = Count, "Receptor" = Gene, "Receptor_source" = Celltype)
                      
                        
                      les_network <- left_join(les_network, lig_data)
                      les_network <- left_join(les_network, rec_data)
                      
                      les_network <- les_network %>% mutate(Connection_product = round(lig_expression*rec_expression, digits = 3))
                    
                      les_network <- les_network[-which(les_network$Connection_product == 0),]
                      
                      les_network <- les_network %>% mutate(log10_Connection_product = round(log10(Connection_product), digits = 3))
                      
                      manual.network <- bind_rows(manual.network, les_network[,-which(colnames(les_network) %in% c("lig_expression", "rec_expression"))])
                      
                      }
                        
                }
                  
                
                #Scale manual network connection product to be similar to cellphoneDB connection product:
                manual.network$Connection_product[which(is.na(manual.network$Connection_product))] <- 0
                manual.network$log10_Connection_product[which(is.na(manual.network$log10_Connection_product))] <- min(manual.network$log10_Connection_product[-which(is.na(manual.network$log10_Connection_product))])
                
                ggplot(manual.network, aes(x = log10_Connection_product)) + geom_density() + theme_classic()
                ggplot(CCI_network, aes(x = log10_Connection_product)) + geom_density() + theme_classic()
                
                CellphoneDB.connection.peak <- density(CCI_network$log10_Connection_product)$x[which.max(density(CCI_network$log10_Connection_product)$y)]
                
                MaxY <- max(density(manual.network$log10_Connection_product)$y[density(manual.network$log10_Connection_product)$x > 2])
                manual.network.peak <- density(manual.network$log10_Connection_product)$x[which(density(manual.network$log10_Connection_product)$y == MaxY)]
                
                peak.difference <- CellphoneDB.connection.peak-manual.network.peak
                manual.network$log10_Connection_product <- manual.network$log10_Connection_product+peak.difference
                
                #Rescale to same standard deviation as CellphoneDB data:
                manual.network$log10_Connection_product <- ((manual.network$log10_Connection_product-CellphoneDB.connection.peak)*(sd(CCI_network$log10_Connection_product)/sd(manual.network$log10_Connection_product))+CellphoneDB.connection.peak)
              
                manual.network <- manual.network %>% mutate(Connection_product = 10^log10_Connection_product)
                  
                #Reformat CCI_network:
                CCI_network$Lesion <- sapply(CCI_network$Ligand_Source, FUN = function(x) unlist(str_split(x, pattern = "__"))[1])
                CCI_network$Ligand_Source <- sapply(CCI_network$Ligand_Source, FUN = function(x) unlist(str_split(x, pattern = "__"))[2])
                CCI_network$Receptor_source <- sapply(CCI_network$Receptor_source, FUN = function(x) unlist(str_split(x, pattern = "__"))[2])
                
                #Bind manual network to cellphoneDB network:
                CCI_network <- bind_rows(CCI_network, manual.network)
                
                
          # Last look at nichenet genes not in cellphonedb:
          ligands.not.in.CellphoneDB <- unique(Nichenet_receptor_induced_ligands_network$ligand)[-which(unique(Nichenet_receptor_induced_ligands_network$ligand) %in% c(CCI_network$Ligand, deconvoluted.subunits$complex_name))]
          
              #remove CDH2, CDH3,and FN1 -- all 3 are ligand-ligand adhesion interactions that were filtered from cellphoneDB network.
              Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[-which(Nichenet_receptor_induced_ligands_network$ligand %in% c("FN1", "CDH2", "CDH3")),]

              #Remaining genes look excluded because of lack of receptor expression for both receptor subunits:
              Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[-which(Nichenet_receptor_induced_ligands_network$ligand %in% ligands.not.in.CellphoneDB),]
              

  # Merge CellphoneDB Network with Nichenet second layer network:
                CCI_network$database_source <- "CellphoneDB"
                Nichenet_receptor_induced_ligands_network$database_source <- "Nichenet"
                
                # Reformat Nichenet dataframe:
                
                      # Receptors column
                      Nichenet_receptor_induced_ligands_network$Receptor <- NA
                      
                      ligands <- unique(Nichenet_receptor_induced_ligands_network$ligand)
                      
                      for (l in 1:length(ligands)) {
                        
                        nichenet.ligand <- ligands[l]
                        
                        nichenet.ligand.rows <- Nichenet_receptor_induced_ligands_network[which(Nichenet_receptor_induced_ligands_network$ligand == nichenet.ligand),]
                        
                        corrected.rows <- nichenet.ligand.rows
                        
                        receptors <- CCI_network %>% filter(Ligand == nichenet.ligand) %>% dplyr::select(Receptor) %>% unlist() %>% unique()
                        
                        for (r in 1:length(receptors)) {
                          
                          corrected.receptor.rows <- corrected.rows
                          corrected.receptor.rows$Receptor <- receptors[r]
                          
                          Nichenet_receptor_induced_ligands_network <- bind_rows(Nichenet_receptor_induced_ligands_network, corrected.receptor.rows)
                        }
                        
                      }
                      
                      Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[-which(is.na(Nichenet_receptor_induced_ligands_network$Receptor)),]
                      
                      # Remove Sender cell Ligand column:
                      Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[,-which(colnames(Nichenet_receptor_induced_ligands_network) == "ligand")]
                      
                      # Rename columns to Ligand:
                      colnames(Nichenet_receptor_induced_ligands_network)[which(colnames(Nichenet_receptor_induced_ligands_network) == "target")] <- "Ligand"
                      colnames(Nichenet_receptor_induced_ligands_network)[which(colnames(Nichenet_receptor_induced_ligands_network) == "Condition")] <- "Lesion"
                      
                      #Reformat Cell Sources:
                      Nichenet_receptor_induced_ligands_network$Ligand_Source <- Nichenet_receptor_induced_ligands_network$Receiver
                      Nichenet_receptor_induced_ligands_network$Receptor_source <- Nichenet_receptor_induced_ligands_network$Receiver
                      
                      Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network[,-which(colnames(Nichenet_receptor_induced_ligands_network) == "Receiver")]
                      
                #Scale Nichenet_weight to be similar to cellphoneDB connection product:
                colnames(Nichenet_receptor_induced_ligands_network)[which(colnames(Nichenet_receptor_induced_ligands_network) == "weight")] <- "Connection_product"
                Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network %>% mutate(log10_Connection_product = log10(Connection_product))
                
                
                
                ggplot(Nichenet_receptor_induced_ligands_network, aes(x = log10_Connection_product)) + geom_density() + theme_classic()
                ggplot(CCI_network, aes(x = log10_Connection_product)) + geom_density() + theme_classic()
                
                
                ###### Filter out weak Nichenet connections: (threshold chosen by looking at IL4/Il13 induction of CXCL9 which seems very weak)
                #Nichenet_receptor_induced_ligands_network <- Nichenet_receptor_induced_ligands_network %>% filter(log10_Connection_product >= -2.4)
                
                
                CellphoneDB.connection.peak <- density(CCI_network$log10_Connection_product)$x[which.max(density(CCI_network$log10_Connection_product)$y)]
                
                MaxY <- max(density(Nichenet_receptor_induced_ligands_network$log10_Connection_product)$y[density(Nichenet_receptor_induced_ligands_network$log10_Connection_product)$x > -2.4])
                #MaxY <- max(density(Nichenet_receptor_induced_ligands_network$log10_Connection_product)$y[density(Nichenet_receptor_induced_ligands_network$log10_Connection_product)$x > -2.1])
                Nichenet.strongest.connection.peak <- density(Nichenet_receptor_induced_ligands_network$log10_Connection_product)$x[which(density(Nichenet_receptor_induced_ligands_network$log10_Connection_product)$y == MaxY)]
                
                peak.difference <- CellphoneDB.connection.peak-Nichenet.strongest.connection.peak
                Nichenet_receptor_induced_ligands_network$log10_Connection_product <- Nichenet_receptor_induced_ligands_network$log10_Connection_product+peak.difference
                
                #Rescale to same standard deviation as CellphoneDB data:
                Nichenet_receptor_induced_ligands_network$log10_Connection_product <- ((Nichenet_receptor_induced_ligands_network$log10_Connection_product-CellphoneDB.connection.peak)*(sd(CCI_network$log10_Connection_product)/sd(Nichenet_receptor_induced_ligands_network$log10_Connection_product))+CellphoneDB.connection.peak)
                
                # Then bind_rows
                CCI_network <- CCI_network[,-which(colnames(CCI_network) == "interaction_pair")]
                
                merged.network <- bind_rows(CCI_network, Nichenet_receptor_induced_ligands_network)
                
                ggplot(merged.network, aes(x = log10_Connection_product, col = database_source)) + geom_density(size = 1) + theme_classic()
                
                #save(merged.network, file = "Data/RL_network2/merged.network_frac0.001.Rdata")
            
                
                
#######################################################
# Network Visualization:
#######################################################                        
library(Seurat)
library(tibble)
library(dplyr)
library(umap)
library(FNN)
library(igraph)
library(stringr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library (purrr)
library(data.table)
library(readxl)
library(tidygraph)
library(ggraph)
library(visNetwork)
library(graphlayouts)
library(ggforce)
library(nichenetr)     
library(ggnewscale)
                
    
      setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")
                
      #Load Single Cell Data:
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
                    
                    # Aggregate count data:
                    DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Lesion", "Semi_Detailed_Celltype"))
                    
      #Load CellphoneDB and Nichenet reference tables:
              
          #CellphoneDB Network:
          network.file.path = "Data/RL_network2/out/"
          deconvoluted <- fread(paste(network.file.path, "deconvoluted.txt", sep = ""))
            
          #Nichenet network:
          ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
          lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
          
              # Mike Edit: Only KEGG connections seem legit:
              lr_network = lr_network %>% mutate(bonafide = database %in% c("kegg"))
              lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
              lr_network = lr_network %>% filter(bonafide == TRUE)
          
          
      #Load merged network:          
      load("Data/RL_network2/merged.network_frac0.001_filtered.Rdata")
      
      #Remove homodimers because it's weird for a directed network:
      if (any(merged.network$Ligand == merged.network$Receptor)) {
        merged.network <- merged.network %>% filter(Ligand != Receptor)
      }

      network <- merged.network
      
      #Organize to nodes and edges:
      network <- network %>% mutate(from = case_when(database_source == "CellphoneDB" ~ paste(Ligand_Source, Ligand, sep = ":"),
                                                     database_source == "Nichenet" ~ paste(Receptor_source, Receptor, sep = ":")),
                                    to = case_when(database_source == "CellphoneDB" ~ paste(Receptor_source, Receptor, sep = ":"),
                                                   database_source == "Nichenet" ~ paste(Ligand_Source, Ligand, sep = ":")))
      
      network <- as.data.frame(network)
      
      network$edge_ID <- seq(1:nrow(network))
      rownames(network) <- seq(1:nrow(network))          
                
      

      #Create simple node df by celltype:
                        lesion.list <- unique(network$Lesion)                                                 #Edit for sample:
                        
                        node_list <- c()
                        
                        for (i in 1:length(lesion.list)) {
                          
                          #Filter Edges:
                          lesion.edges <- filter(network, Lesion == lesion.list[i])                     #Edit for sample:
                          
                          #Create Nodes DFs:
                          
                          nodes <- unique(c(lesion.edges$Ligand_Source, lesion.edges$Receptor_source))
                          
                          nodes.df <- data.frame(node_name = nodes,
                                                 CellType = nodes,
                                                 Lesion = lesion.list[i])
                          
                          #Assign Outputs:
                          assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
                          assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
                          
                          node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
                          
                        }                
                        
                        #Merge Nodes:
                        nodes.df <- purrr::reduce(mget(node_list), full_join)
                        
                        nodes.df$Lesion <- factor(nodes.df$Lesion, 
                                                  levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
                        
                        #Create Edges.df:
                        edges.df <- network                                                                  #Edit for sample:   
                        edges.df$from <- str_replace_all(edges.df$from, pattern = ":.*$", replacement = "")
                        edges.df$to <- str_replace_all(edges.df$to, pattern = ":.*$", replacement = "")
                        
                        edges.df$Lesion <- factor(edges.df$Lesion, 
                                                  levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
                        
                        edges.df <- edges.df[,which(colnames(edges.df) %in% c("from", "to", "Lesion", "Ligand_Source", "Ligand", "Receptor_source", "Receptor","log10_Connection_product", "database_source"))]
                        
                        
                        #Remove duplicate edges: (not sure why this happened. . . need to fix this eventually)
                        edges.df <- edges.df %>% mutate(edgeID = paste(from, "__", to, "__", Lesion, sep = ""))
                        
                        duplicates <-  edges.df$edgeID[which(duplicated(edges.df$edgeID))]         
                        duplicates <- edges.df[which(edges.df$edgeID %in% duplicates),]
                        duplicate.max.connections <- duplicates %>% group_by(edgeID) %>% summarize(max.connection = max(log10_Connection_product))
                        
                        duplicates <- duplicates[-which(duplicated(duplicates$edgeID)),]
                        duplicates <- full_join(duplicates, duplicate.max.connections)
                        duplicates <- duplicates[,-which(colnames(duplicates) == "log10_Connection_product")]
                        
                        duplicates <- duplicates %>% rename(log10_Connection_product = "max.connection")
                        
                        edges.df <- edges.df[-which(edges.df$edgeID %in% duplicates$edgeID),]
                        edges.df <- bind_rows(edges.df, duplicates)
                        
                        #RL_graph <- tbl_graph(nodes = nodes.df, edges = edges.df, directed = TRUE, node_key = "node_name")
                        
                        #Simple Graph:
                        simple.edges.df <- edges.df %>% filter(database_source == "CellphoneDB") %>% group_by(Ligand_Source, Receptor_source, from, to, Lesion) %>% summarize(count = n())
                        simple.edges.df$database_source <- "CellphoneDB"
                        Simple_RL_graph <- tbl_graph(nodes = nodes.df, edges = simple.edges.df, directed = TRUE, node_key = "node_name")
                        
                        
                              #Correct misconnected edges: 
                              graph.nodes <- Simple_RL_graph %>% activate(nodes) %>% as.data.frame()
                              graph.edges <- Simple_RL_graph %>% activate(edges) %>% as.data.frame()
                              
                              duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
                              
                              nodes.df <- rownames_to_column(nodes.df)
                              
                              
                              for (r in 1:nrow(duplicates)) {
                                
                                les <- as.character(duplicates[r,"Lesion"])
                                
                                if (duplicates[r,"database_source"] == "CellphoneDB") {
                                  from <- duplicates[r,"Ligand_Source"]
                                  
                                  to <- duplicates[r,"Receptor_source"]
                                  
                                  duplicates[r,"from"] <- nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.numeric()
                                  duplicates[r,"to"] <- nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.numeric()
                                  
                                }
                                
                                if (duplicates[r,"database_source"] == "Nichenet") {
                                  from <- duplicates[r,"Receptor_source"]
                                  
                                  to <- duplicates[r,"Ligand_Source"]
                                  
                                  duplicates[r,"from"] <- nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.numeric()
                                  duplicates[r,"to"] <- nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.numeric()
                                  
                                }
                              }
                              
                              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
                              
                              Simple_RL_graph <- Simple_RL_graph %>% activate(edges) %>% right_join(graph.edges)
                              
                              # Define Layout using merged data from all lesions:
                              nodes <- nodes.df[,which(colnames(nodes.df) %in% c("node_name", "CellType"))]
                              nodes <- nodes[-which(duplicated(nodes$node_name)),]
                              
                              merged.graph <- tbl_graph(nodes = nodes, edges = simple.edges.df, directed = TRUE, node_key = "node_name")
                              
                              layout <- layout_with_stress(merged.graph)
                              colnames(layout) <- c("x", "y")
                              
                              layout <- cbind(nodes, layout)
                              
                              #Copy coordinates to node df containing all faceted nodes:
                              layout_full <- full_join(nodes.df, layout)
                              
                        
                        #Faceted Graph:
                        ggraph(Simple_RL_graph, layout = "manual",
                               x = layout_full[, "x"],
                               y = layout_full[, "y"]) + 
                          geom_edge_fan(aes(width = count), alpha = 0.85, arrow = arrow(angle =20, length = unit(3, "mm"), type = "closed"), end_cap = circle(5,'mm')) + 
                          geom_edge_loop(aes(width = count, strength = 0.25), alpha = 0.85, arrow = arrow(angle =20, length = unit(3, "mm"), type = "closed"), end_cap = circle(5,'mm')) +
                          geom_node_point(aes(color = CellType), size = 7) + 
                          scale_edge_width(range = c(0.01, 3)) +
                          theme_graph() + 
                          theme(legend.position = "none") +
                          facet_nodes(~Lesion) +
                          geom_node_label(aes(label = as.character(node_name)), size = 5, repel = TRUE, point.padding = unit(0.2, "lines"), fill = "white")
                        
                           
                                          
      # Ligand-Receptor Celltype-Specific Node Plot:        

        contact_derm_DE_filtered <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE_filtered.txt")

        DE.table <- contact_derm_DE_filtered
              
        
              #Add Basic_Celltype to DE.table for network shapes
              DE.table <- DE.table %>% mutate(Basic_Celltype = case_when(CellType == "KRT" ~ "KRT",
                                                                         CellType == "MEL" ~ "MEL",
                                                                         CellType %in% c("Myeloid", "pDC", "LC") ~ "APC",
                                                                         CellType %in% c("CD4", "CD8", "NK", "Treg") ~ "LYMPH"))
              DE.table$Basic_Celltype <- factor(DE.table$Basic_Celltype, levels = c("KRT", "MEL", "APC", "LYMPH"))
              DE.table$CellType <- factor(DE.table$CellType, levels = c("KRT", "MEL", "LC", "Myeloid", "pDC", "CD4", "Treg", "CD8", "NK"))
              
        
        plot_lig_rec_network <- function(network, DE.table, seurat.obj, facet_by = "Lesion", node_shape = "Basic_Celltype", node_shape_vals = c(15, 8, 13, 16), 
                                         create_DAG = TRUE, plot_layout = "sugiyama") {
          
          # Ligand-Receptor-Celltype Specific Node Plot:        
          message("Making initial ligand-receptor celltype-specific network. . . ")
          #Create node df:
          lesion.list <- unique(network[, facet_by])                                                 
          
          node_list <- c()
          
          for (i in 1:length(lesion.list)) {
            
            #Filter Edges:
            lesion.edges <- filter(network, eval(as.symbol(facet_by)) == lesion.list[i])                     
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet_by = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          #Merge Nodes:
          nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            nodes.df[,facet_by] <- factor(nodes.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
          }
          
          #Create Edges.df:
          edges.df <- network                                                                          
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            edges.df[,facet_by] <- factor(edges.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          edges.df <- edges.df[,which(colnames(edges.df) %in% c("from", "to", facet_by, "Ligand_Source", "Ligand", "Receptor_source", "Receptor","log10_Connection_product", "database_source"))]
          
          #Remove duplicate nichenet edges: (not sure why this happened. . . need to fix this eventually)
          edges.df <- edges.df %>% mutate(edgeID = paste(from, "__", to, "__", eval(as.symbol(facet_by)), sep = ""))
          
          duplicates <-  edges.df$edgeID[which(duplicated(edges.df$edgeID))]         
          duplicates <- edges.df[which(edges.df$edgeID %in% duplicates),]
          duplicate.max.connections <- duplicates %>% group_by(edgeID) %>% summarize(max.connection = max(log10_Connection_product))
          
          duplicates <- duplicates[-which(duplicated(duplicates$edgeID)),]
          duplicates <- full_join(duplicates, duplicate.max.connections)
          duplicates <- duplicates[,-which(colnames(duplicates) == "log10_Connection_product")]
          
          duplicates <- duplicates %>% rename(log10_Connection_product = "max.connection")
          
          edges.df <- edges.df[-which(edges.df$edgeID %in% duplicates$edgeID),]
          edges.df <- bind_rows(edges.df, duplicates)
          
          # Scale edge width by percent of max connection strength:
          message("Scaling edges. . . ")
          
          ### Why is Day4 allergy Myeloid:CXCL9-KRT:CXCR3 almost as strong as Day4 Allergy lymphocyte CXCR3 connections?
          
          max.expression.vals <- edges.df %>% group_by(Ligand, Receptor) %>% summarize(max.expression = max(log10_Connection_product))
          
          scaled_edges <- edges.df
          
          for (r in 1:nrow(max.expression.vals)) {
            
            connection.row <- max.expression.vals[r,]
            connection.max <- max.expression.vals[r,"max.expression"]
            
            connection.edges <- edges.df %>% filter(Ligand == connection.row$Ligand & Receptor == connection.row$Receptor)
            connection.edges$log10_Connection_product <- sapply(connection.edges$log10_Connection_product, FUN = function(x) as.numeric(round(10^(x)/(10^(connection.max)), digits = 3)))
            
            scaled_edges[which(scaled_edges$edgeID %in% connection.edges$edgeID),] <- connection.edges
          }
          
          edges.df <- scaled_edges
          colnames(edges.df)[which(colnames(edges.df) == "log10_Connection_product")] <- "Scaled_Connection_Product"
          
          
          #Add edge DE status:
          edges.df$DE <- FALSE
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            les.celltypes <- unique(lesion.edges$Ligand_Source)
            
            for (lc in 1:length(les.celltypes)) {
              
              les.cell <- les.celltypes[lc]
              
              les.cell.edges <- lesion.edges %>% filter(Ligand_Source == les.cell)
              
              les.cell.DE <- DE.table %>% filter(Comparison == les & CellType == les.cell & Color_factor == "Significantly Up") %>% dplyr::pull(Gene)
              
              DE.complex.genes <- deconvoluted %>% filter(gene_name %in% les.cell.DE) %>% dplyr::pull(complex_name) %>% unique()
              
              #Only label DE CellphoneDB
              les.cell.edges.DE <- les.cell.edges %>% filter(Ligand %in% unique(les.cell.DE, DE.complex.genes) & database_source == "CellphoneDB") %>% dplyr::pull(edgeID)
              
              if (length(les.cell.edges.DE) >0) {
                edges.df[which(edges.df$edgeID %in% les.cell.edges.DE), "DE"] <- TRUE
              }
              
            }
          }
          
          
          #Prune edges to only interconnected CellphoneDB & Nichenet edges, plus all DE CellphoneDB edges:
          message("Pruning edges to interconnected cellphonedb and nichenet edges. . . ")
          
          filtered.edges.df <- as.data.frame(edges.df)[0,]
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            if (nrow(lesion.edges)>0) {
              
              #Prune nichenet edges:
              lesion.nichenet.edges <- lesion.edges %>% filter(database_source == "Nichenet" & Scaled_Connection_Product > 0) 
              
              if (nrow(lesion.nichenet.edges)>0){
                
                nichenet.edges.to.prune <- c()
                
                for (ne in 1:length(unique(lesion.nichenet.edges$from))) {
                  
                  nichenet.edge.from <- unique(lesion.nichenet.edges$from)[ne]
                  
                  upstream.edges <- lesion.edges %>% filter(to == nichenet.edge.from & Scaled_Connection_Product > 0 & database_source == "CellphoneDB")
                  
                  if (nrow(upstream.edges) == 0) {
                    nichenet.edges.to.prune <- c(nichenet.edges.to.prune, nichenet.edge.from)
                  }
                  
                }    
                
                lesion.nichenet.edges <- lesion.nichenet.edges %>% filter(!from %in% nichenet.edges.to.prune)
              }
              
              #Prune CellphoneDB edges:
              lesion.cellphonedb.edges <- lesion.edges %>% filter(database_source == "CellphoneDB" & Scaled_Connection_Product > 0) 
              
              if(nrow(lesion.cellphonedb.edges) >0) {
                
                cellphonedb.edges.to.prune <- c()
                
                for (ce in 1:length(unique(lesion.cellphonedb.edges$to))) {
                  
                  cellphonedb.edge.to <- unique(lesion.cellphonedb.edges$to)[ce]
                  
                  downstream.edges <- lesion.edges %>% filter(from == cellphonedb.edge.to & Scaled_Connection_Product > 0 & database_source == "Nichenet")
                  
                  if (nrow(downstream.edges) == 0) {
                    cellphonedb.edges.to.prune <- c(cellphonedb.edges.to.prune, cellphonedb.edge.to)
                  }
                  
                }
                
                lesion.cellphonedb.edges <- lesion.cellphonedb.edges %>% filter(!to %in% cellphonedb.edges.to.prune)
              }
              
              lesion.edges <- bind_rows(lesion.cellphonedb.edges, lesion.nichenet.edges)         
              
              filtered.edges.df <- bind_rows(filtered.edges.df, lesion.edges)
            }
          }
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.edges.df[,facet_by] <- factor(filtered.edges.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          #Duplicate node names columns so that each is preserved after graph formation:
          filtered.edges.df$from.named <- filtered.edges.df$from
          filtered.edges.df$to.named <- filtered.edges.df$to
          
          
          #Create filtered.nodes.df:
          node_list <- c()
          
          for (i in 1:length(unique(filtered.edges.df[,facet_by]))) {
            
            les <- as.character(unique(filtered.edges.df[,facet_by]))[i]
            #Filter Edges:
            lesion.edges <- filter(filtered.edges.df, eval(as.symbol(facet_by)) == les)                   
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          
          #Merge Nodes:
          filtered.nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.nodes.df[,facet_by] <- factor(filtered.nodes.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          
          #Make Initial Graph:    
          RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          #Correct misconnected edges:
          message("Correcting tbl_graph misconnected edges. . . ")
          
          graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
          graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
          
          duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
          
          filtered.nodes.df <- rownames_to_column(filtered.nodes.df)
          
          for (r in 1:nrow(duplicates)) {
            
            les <- as.character(duplicates[r,facet_by])
            
            if (duplicates[r,"database_source"] == "CellphoneDB") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
            
            if (duplicates[r,"database_source"] == "Nichenet") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
          }
          
          graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
          
          RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
          
          
          # Define Layout using merged data from all lesions:
          nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
          nodes <- nodes[-which(duplicated(nodes$node_name)),]  
          
          merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          layout <- layout_with_stress(merged.graph)
          colnames(layout) <- c("x", "y")
          
          layout <- cbind(nodes, layout)
          
          ## Ligand Edge Label mapping:
          label.mapping.df <- filtered.edges.df
          label.mapping.df$label_x <- NA
          label.mapping.df$label_y <- NA
          label.mapping.df$Edge_label <- ""
          
          for (l in 1:nrow(label.mapping.df)) {
            
            lig.from <- label.mapping.df$from.named[l]
            lig.to <- label.mapping.df$to.named[l]
            
            node1_info <- layout %>% filter(node_name == lig.from)
            node2_info <- layout %>% filter(node_name == lig.to)
            
            label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
            label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
          }
          
          #Only label ligands that are reasonably far away from each other, with DE prioritization:
          max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
          
          for (l in 1:length(unique(label.mapping.df$Ligand))) {
            lig <- unique(label.mapping.df$Ligand)[l]
            
            lig.edges <- label.mapping.df %>% filter(Ligand == lig)
            
            lig.edges$Edge_label <- "temp"
            
            while(any(lig.edges$Edge_label =="temp")){
              
              #temporarily make a column for fraction of max distance from first ungrouped edge:
              grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
              
              grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.035) 
              
              #Label closest to centroid edges grouped by 5% of max layout distance:
              centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
              centroid_x <- sum(centroid_x)/length(centroid_x)
              centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
              centroid_y <- sum(centroid_y)/length(centroid_y)
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
              
              group.edges <- grouped.lig.edges %>% arrange(desc(DE), centroid_dist) %>% dplyr::pull(edgeID)
              Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
              
              lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
              
              lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
              
            }
            
            label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
            
          }
          
          #Copy coordinates to node df containing all faceted nodes:
          layout_full <- full_join(filtered.nodes.df, layout)
          
          if (create_DAG == TRUE)  {
            message("Pruning to directed acyclic graph. . . ")
            
            #Prune cycles and feedback loops to create a directed acyclic graph:
            
            ## Find Cycles in network:
            FindCycles = function(g) {
              Cycles = NULL
              for(v1 in V(g)) {
                if(degree(g, v1, mode="in") == 0) { next }
                GoodNeighbors = neighbors(g, v1, mode="out")
                GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
                for(v2 in GoodNeighbors) {
                  TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
                  TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
                  TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
                  Cycles  = c(Cycles, TempCyc)
                }
              }
              Cycles
            }
            merged.igraph <- as.igraph(merged.graph)
            network.cycles <- FindCycles(merged.igraph)
            
            nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
            
            for (c in 1:length(network.cycles)) {
              cycle <- network.cycles[[c]]
              
              named.cycle <- c()  
              
              for (n in 1:length(cycle)) {
                node.num <- cycle[n]
                node.info <- nodes[as.numeric(node.num),]
                named.node <- node.info$node_name
                
                named.cycle <- c(named.cycle, named.node)
              }
              
              if (c == 1) {
                named.network.cycles <-list(named.cycle)
              }
              if (c > 1) {
                named.network.cycles <- append(named.network.cycles, list(named.cycle))
              }
            }
            
            #Edges causing cycles:
            
            for (c in 1:length(named.network.cycles)) {
              named.cycle <- named.network.cycles[[c]]
              
              first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
              last.edges <- filtered.edges.df %>% filter(from == named.cycle[length(named.cycle)-1] & to == named.cycle[length(named.cycle)])
              
              if (c ==1) {
                cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                  first.edge.from = unique(first.edges$from.named),
                                                  first.edge.to = unique(first.edges$to.named),
                                                  first.edges.num= nrow(first.edges),
                                                  first.edges.database = unique(first.edges$database_source),
                                                  last.edge.from = unique(last.edges$from.named),
                                                  last.edge.to = unique(last.edges$to.named),
                                                  last.edges.num = nrow(last.edges),
                                                  last.edges.database = unique(last.edges$database_source))
              }
              
              if (c > 1) {
                cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                               first.edge.from = unique(first.edges$from.named),
                                               first.edge.to = unique(first.edges$to.named),
                                               first.edges.num= nrow(first.edges),
                                               first.edges.database = unique(first.edges$database_source),
                                               last.edge.from = unique(last.edges$from.named),
                                               last.edge.to = unique(last.edges$to.named),
                                               last.edges.num = nrow(last.edges),
                                               last.edges.database = unique(last.edges$database_source))
                cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
              }
              
            }
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(first.edge.from, "__", first.edge.to, sep = ""),
                                          last.edge = paste(last.edge.from, "__", last.edge.to, sep = ""))
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, cycle_id = paste(first.edge, ".....",last.edge, sep = ""))
            
            #Cycles seem to be due to weak nichenet connections at last edge:
            cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$last.edge)),]
            cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("last.edge.from", "last.edge.to"))]
            
            if (nrow(cyclic_edges_to_prune) >0 ) {
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$last.edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$last.edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            
            #Re-make Graph without above pruned cycle edges:    
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, layout)
            
            ## Above function doesn't find all cycles though! Use igraph to find feedback edges requiring pruning for DAG creation:
            merged.igraph <- as.igraph(merged.graph)
            cycles <- feedback_arc_set(merged.igraph)
            cycle.edge.ids <- as_ids(cycles)
            
            if (length(cycle.edge.ids)>0){
              
              #convert to list:
              for (e in 1:length(cycle.edge.ids)) {
                if (e ==1){
                  cyclic.edges <- list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]])))
                  
                }
                
                if (e>1){
                  cyclic.edges <-  append(cyclic.edges, list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]]))))
                }
              }   
              
              #Add names:
              for (c in 1:length(cyclic.edges)) {
                cycle <- cyclic.edges[[c]]
                
                named.cycle <- c()  
                
                for (n in 1:length(cycle)) {
                  node.num <- cycle[n]
                  node.info <- nodes[as.numeric(node.num),]
                  named.node <- node.info$node_name
                  
                  named.cycle <- c(named.cycle, named.node)
                }
                
                if (c == 1) {
                  named.network.cycles <-list(named.cycle)
                }
                if (c > 1) {
                  named.network.cycles <- append(named.network.cycles, list(named.cycle))
                }
              }
              
              for (c in 1:length(named.network.cycles)) {
                named.cycle <- named.network.cycles[[c]]
                
                first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
                
                if (c ==1) {
                  cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                    edge.from = unique(first.edges$from.named),
                                                    edge.to = unique(first.edges$to.named),
                                                    edges.num= nrow(first.edges),
                                                    edges.database = unique(first.edges$database_source))
                  
                }
                
                if (c > 1) {
                  cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                 edge.from = unique(first.edges$from.named),
                                                 edge.to = unique(first.edges$to.named),
                                                 edges.num= nrow(first.edges),
                                                 edges.database = unique(first.edges$database_source))
                  
                  cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                }
                
              }
              
              cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(edge.from, "__", edge.to, sep = ""))
              
              #All cycles look to be due to weak nichenet connections at last edge:
              cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$first.edge)),]
              cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("edge.from", "edge.to"))]
              
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            # Make Final Graph:    
            message("Making final graph. . . ")
            
            #add shapes to nodes if you want:
            if (!is.null(node_shape)) {
              
              node.shape.df <- DE.table[,which(colnames(DE.table) %in% c("CellType", node_shape))]  
              node.shape.df <- node.shape.df[-which(duplicated(node.shape.df$CellType)),]
              
              filtered.nodes.df <- left_join(filtered.nodes.df, node.shape.df)
              
              if (class(seurat.obj@meta.data[,node_shape]) == "factor") {
                filtered.nodes.df[node_shape] <- factor(filtered.nodes.df[,node_shape], 
                                                        levels = levels(seurat.obj@meta.data[,node_shape]))
              }
            }
            
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType", node_shape, "node_label"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            ## Ligand Edge Label mapping:
            label.mapping.df <- filtered.edges.df
            label.mapping.df$label_x <- NA
            label.mapping.df$label_y <- NA
            label.mapping.df$Edge_label <- ""
            
            for (l in 1:nrow(label.mapping.df)) {
              
              lig.from <- label.mapping.df$from.named[l]
              lig.to <- label.mapping.df$to.named[l]
              
              node1_info <- layout %>% filter(node_name == lig.from)
              node2_info <- layout %>% filter(node_name == lig.to)
              
              label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
              label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
            }
            
            #Only label ligands that are reasonably far away from each other, with DE prioritization:
            max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
            
            for (l in 1:length(unique(label.mapping.df$Ligand))) {
              lig <- unique(label.mapping.df$Ligand)[l]
              
              lig.edges <- label.mapping.df %>% filter(Ligand == lig)
              
              lig.edges$Edge_label <- "temp"
              
              while(any(lig.edges$Edge_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.edges <- grouped.lig.edges %>% arrange(desc(DE),centroid_dist) %>% dplyr::pull(edgeID)
                Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
                
                lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                
                lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                
              }
              
              label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
              
            }
            
            ## Node Labels:
            node.labels.df <- layout
            node.labels.df$node_label <- ""
            
            node.labels.df$node_gene <- sapply(node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            #Only label nodes that are reasonably far away from each other:
            
            for (g in 1:length(unique(node.labels.df$node_gene))) {
              gene <- unique(node.labels.df$node_gene)[g]
              
              gene.nodes <- node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.035) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                
              }
              
              node.labels.df[which(node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, node.labels.df)
            
          }
          
          #Set point colors for plotting:
          GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
          
          ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
          node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
          
          if (plot_layout == "stress") {
            
            ##Plot With Facets:
            g <- ggraph(RL_graph, layout = "manual",
                        x = layout_full[, "x"],
                        y = layout_full[, "y"]) + 
              geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
              geom_node_point(aes(color = CellType), size = 0.5) + 
              geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
              scale_edge_width(range = c(0.1, 0.5)) +
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              facet_nodes(~Lesion) +
              scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
              new_scale_color()+
              scale_color_manual(values = c("#000000", "#750000")) +
              geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0)
            
            network.output.result <- list("Graph_object" = list(RL_graph),
                                          "Layout" = layout_full,
                                          "Plot" = g)
            
          }
          
          if (plot_layout == "sugiyama") {
            #Sugiyama Layout:
            sugi_layout <- layout_with_sugiyama(merged.graph, hgap = 500)
            sugi_layout <- as.data.frame(sugi_layout$layout)
            colnames(sugi_layout) <- c("x", "y")
            
            sugi_layout <- cbind(nodes, sugi_layout)
            
            ## Sugi Node Labels:
            sugi.node.labels.df <- sugi_layout
            sugi.node.labels.df$node_label <- ""
            
            sugi.node.labels.df$node_gene <- sapply(sugi.node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            max.sugi.layout.dist <- ((max(sugi_layout$x)-min(sugi_layout$x))^2+(max(sugi_layout$y)-min(sugi_layout$y))^2)^(0.5)
            
            #Only label nodes that are reasonably far away from each other:
            for (g in 1:length(unique(sugi.node.labels.df$node_gene))) {
              gene <- unique(sugi.node.labels.df$node_gene)[g]
              
              gene.nodes <- sugi.node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.sugi.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                
              }
              
              sugi.node.labels.df[which(sugi.node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            sugi_layout_full <- full_join(filtered.nodes.df, sugi.node.labels.df)
            
            g <- ggraph(RL_graph, layout = "manual",
                        x = sugi_layout_full[, "x"],
                        y = sugi_layout_full[, "y"]) + 
              geom_node_point(data = sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
              geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
              scale_shape_manual(values = node_shape_vals)+
              scale_color_manual(values = node.colors1) +
              geom_edge_diagonal(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
              scale_edge_width(range = c(0.1, 0.5)) +
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              facet_nodes(~Lesion) +
              scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
              expand_limits(y = c(min(sugi.node.labels.df$y), max(sugi.node.labels.df$y)+0.05*(range(sugi.node.labels.df$y)[2]-range(sugi.node.labels.df$y)[1]))) +
              geom_text_repel(data = sugi_layout_full, aes(x = x, y = y, label = node_label), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5)
            
            network.output.result <- list("Graph_object" = list(RL_graph),
                                          "Layout" = sugi_layout_full,
                                          "Plot" = g)
          }
          
          plot(g)
          return(network.output.result)
        }
        
        
        plot_lig_rec_merged_network <- function(network, DE.table, seurat.obj, condition_of_interest,
                                                facet_by = "Lesion", node_shape = "Basic_Celltype", node_shape_vals = c(15, 8, 13, 16), 
                                                create_DAG = TRUE, plot_layout = "sugiyama") {
          
          # Ligand-Receptor-Celltype Specific Node Plot:        
          message("Making initial ligand-receptor celltype-specific network. . . ")
          #Create node df:
          lesion.list <- unique(network[, facet_by])                                                 
          
          node_list <- c()
          
          for (i in 1:length(lesion.list)) {
            
            #Filter Edges:
            lesion.edges <- filter(network, eval(as.symbol(facet_by)) == lesion.list[i])                     
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet_by = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          #Merge Nodes:
          nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            nodes.df[,facet_by] <- factor(nodes.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
          }
          
          #Create Edges.df:
          edges.df <- network                                                                          
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            edges.df[,facet_by] <- factor(edges.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          edges.df <- edges.df[,which(colnames(edges.df) %in% c("from", "to", facet_by, "Ligand_Source", "Ligand", "Receptor_source", "Receptor","log10_Connection_product", "database_source"))]
          
          #Remove duplicate nichenet edges: (not sure why this happened. . . need to fix this eventually)
          edges.df <- edges.df %>% mutate(edgeID = paste(from, "__", to, "__", eval(as.symbol(facet_by)), sep = ""))
          
          duplicates <-  edges.df$edgeID[which(duplicated(edges.df$edgeID))]         
          duplicates <- edges.df[which(edges.df$edgeID %in% duplicates),]
          duplicate.max.connections <- duplicates %>% group_by(edgeID) %>% summarize(max.connection = max(log10_Connection_product))
          
          duplicates <- duplicates[-which(duplicated(duplicates$edgeID)),]
          duplicates <- full_join(duplicates, duplicate.max.connections)
          duplicates <- duplicates[,-which(colnames(duplicates) == "log10_Connection_product")]
          
          duplicates <- duplicates %>% rename(log10_Connection_product = "max.connection")
          
          edges.df <- edges.df[-which(edges.df$edgeID %in% duplicates$edgeID),]
          edges.df <- bind_rows(edges.df, duplicates)
          
          # Scale edge width by percent of max connection strength:
          message("Scaling edges. . . ")
          
          ### Why is Day4 allergy Myeloid:CXCL9-KRT:CXCR3 almost as strong as Day4 Allergy lymphocyte CXCR3 connections?
          
          max.expression.vals <- edges.df %>% group_by(Ligand, Receptor) %>% summarize(max.expression = max(log10_Connection_product))
          
          scaled_edges <- edges.df
          
          for (r in 1:nrow(max.expression.vals)) {
            
            connection.row <- max.expression.vals[r,]
            connection.max <- max.expression.vals[r,"max.expression"]
            
            connection.edges <- edges.df %>% filter(Ligand == connection.row$Ligand & Receptor == connection.row$Receptor)
            connection.edges$log10_Connection_product <- sapply(connection.edges$log10_Connection_product, FUN = function(x) as.numeric(round(10^(x)/(10^(connection.max)), digits = 3)))
            
            scaled_edges[which(scaled_edges$edgeID %in% connection.edges$edgeID),] <- connection.edges
          }
          
          edges.df <- scaled_edges
          colnames(edges.df)[which(colnames(edges.df) == "log10_Connection_product")] <- "Scaled_Connection_Product"
          
          
          #Add edge DE status:
          edges.df$DE <- FALSE
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            les.celltypes <- unique(lesion.edges$Ligand_Source)
            
            for (lc in 1:length(les.celltypes)) {
              
              les.cell <- les.celltypes[lc]
              
              les.cell.edges <- lesion.edges %>% filter(Ligand_Source == les.cell)
              
              les.cell.DE <- DE.table %>% filter(Comparison == les & CellType == les.cell & Color_factor == "Significantly Up") %>% dplyr::pull(Gene)
              
              DE.complex.genes <- deconvoluted %>% filter(gene_name %in% les.cell.DE) %>% dplyr::pull(complex_name) %>% unique()
              
              #Only label DE CellphoneDB
              les.cell.edges.DE <- les.cell.edges %>% filter(Ligand %in% unique(les.cell.DE, DE.complex.genes) & database_source == "CellphoneDB") %>% dplyr::pull(edgeID)
              
              if (length(les.cell.edges.DE) >0) {
                edges.df[which(edges.df$edgeID %in% les.cell.edges.DE), "DE"] <- TRUE
              }
              
            }
          }
          
          
          #Prune edges to only interconnected CellphoneDB & Nichenet edges:
          message("Pruning edges to interconnected cellphonedb and nichenet edges. . . ")
          
          filtered.edges.df <- as.data.frame(edges.df)[0,]
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            if (nrow(lesion.edges)>0) {
              
              #Prune nichenet edges:
              lesion.nichenet.edges <- lesion.edges %>% filter(database_source == "Nichenet" & Scaled_Connection_Product > 0) 
              
              if (nrow(lesion.nichenet.edges)>0){
                
                nichenet.edges.to.prune <- c()
                
                for (ne in 1:length(unique(lesion.nichenet.edges$from))) {
                  
                  nichenet.edge.from <- unique(lesion.nichenet.edges$from)[ne]
                  
                  upstream.edges <- lesion.edges %>% filter(to == nichenet.edge.from & Scaled_Connection_Product > 0 & database_source == "CellphoneDB")
                  
                  if (nrow(upstream.edges) == 0) {
                    nichenet.edges.to.prune <- c(nichenet.edges.to.prune, nichenet.edge.from)
                  }
                  
                }    
                
                lesion.nichenet.edges <- lesion.nichenet.edges %>% filter(!from %in% nichenet.edges.to.prune)
              }
              
              #Prune CellphoneDB edges:
              lesion.cellphonedb.edges <- lesion.edges %>% filter(database_source == "CellphoneDB" & Scaled_Connection_Product > 0) 
              
              if(nrow(lesion.cellphonedb.edges) >0) {
                
                cellphonedb.edges.to.prune <- c()
                
                for (ce in 1:length(unique(lesion.cellphonedb.edges$to))) {
                  
                  cellphonedb.edge.to <- unique(lesion.cellphonedb.edges$to)[ce]
                  
                  downstream.edges <- lesion.edges %>% filter(from == cellphonedb.edge.to & Scaled_Connection_Product > 0 & database_source == "Nichenet")
                  
                  if (nrow(downstream.edges) == 0) {
                    cellphonedb.edges.to.prune <- c(cellphonedb.edges.to.prune, cellphonedb.edge.to)
                  }
                  
                }
                
                lesion.cellphonedb.edges <- lesion.cellphonedb.edges %>% filter(!to %in% cellphonedb.edges.to.prune)
              }
              
              lesion.edges <- bind_rows(lesion.cellphonedb.edges, lesion.nichenet.edges)         
              
              filtered.edges.df <- bind_rows(filtered.edges.df, lesion.edges)
            }
          }
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.edges.df[,facet_by] <- factor(filtered.edges.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          #Duplicate node names columns so that each is preserved after graph formation:
          filtered.edges.df$from.named <- filtered.edges.df$from
          filtered.edges.df$to.named <- filtered.edges.df$to
          
          
          #Create filtered.nodes.df:
          node_list <- c()
          
          for (i in 1:length(unique(filtered.edges.df[,facet_by]))) {
            
            les <- as.character(unique(filtered.edges.df[,facet_by]))[i]
            #Filter Edges:
            lesion.edges <- filter(filtered.edges.df, eval(as.symbol(facet_by)) == les)                   
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          
          #Merge Nodes:
          filtered.nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.nodes.df[,facet_by] <- factor(filtered.nodes.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          
          #Make Initial Graph:    
          RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          #Correct misconnected edges:
          message("Correcting tbl_graph misconnected edges. . . ")
          
          graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
          graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
          
          duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
          
          filtered.nodes.df <- rownames_to_column(filtered.nodes.df)
          
          for (r in 1:nrow(duplicates)) {
            
            les <- as.character(duplicates[r,facet_by])
            
            if (duplicates[r,"database_source"] == "CellphoneDB") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
            
            if (duplicates[r,"database_source"] == "Nichenet") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
          }
          
          graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
          
          RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
          
          
          # Define Layout using merged data from all lesions:
          nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
          nodes <- nodes[-which(duplicated(nodes$node_name)),]  
          
          merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          layout <- layout_with_stress(merged.graph)
          colnames(layout) <- c("x", "y")
          
          layout <- cbind(nodes, layout)
          
          ## Ligand Edge Label mapping:
          label.mapping.df <- filtered.edges.df
          label.mapping.df$label_x <- NA
          label.mapping.df$label_y <- NA
          label.mapping.df$Edge_label <- ""
          
          for (l in 1:nrow(label.mapping.df)) {
            
            lig.from <- label.mapping.df$from.named[l]
            lig.to <- label.mapping.df$to.named[l]
            
            node1_info <- layout %>% filter(node_name == lig.from)
            node2_info <- layout %>% filter(node_name == lig.to)
            
            label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
            label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
          }
          
          #Only label ligands that are reasonably far away from each other, with DE prioritization:
          max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
          
          for (l in 1:length(unique(label.mapping.df$Ligand))) {
            lig <- unique(label.mapping.df$Ligand)[l]
            
            lig.edges <- label.mapping.df %>% filter(Ligand == lig)
            
            lig.edges$Edge_label <- "temp"
            
            while(any(lig.edges$Edge_label =="temp")){
              
              #temporarily make a column for fraction of max distance from first ungrouped edge:
              grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
              
              grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.035) 
              
              #Label closest to centroid edges grouped by 5% of max layout distance:
              centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
              centroid_x <- sum(centroid_x)/length(centroid_x)
              centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
              centroid_y <- sum(centroid_y)/length(centroid_y)
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
              
              group.edges <- grouped.lig.edges %>% arrange(desc(DE), centroid_dist) %>% dplyr::pull(edgeID)
              Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
              
              lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
              
              lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
              
            }
            
            label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
            
          }
          
          #Copy coordinates to node df containing all faceted nodes:
          layout_full <- full_join(filtered.nodes.df, layout)
          
          if (create_DAG == TRUE)  {
            message("Pruning to directed acyclic graph. . . ")
            
            #Prune cycles and feedback loops to create a directed acyclic graph:
            
            ## Find Cycles in network:
            FindCycles = function(g) {
              Cycles = NULL
              for(v1 in V(g)) {
                if(degree(g, v1, mode="in") == 0) { next }
                GoodNeighbors = neighbors(g, v1, mode="out")
                GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
                for(v2 in GoodNeighbors) {
                  TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
                  TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
                  TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
                  Cycles  = c(Cycles, TempCyc)
                }
              }
              Cycles
            }
            merged.igraph <- as.igraph(merged.graph)
            network.cycles <- FindCycles(merged.igraph)
            
            nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
            
            for (c in 1:length(network.cycles)) {
              cycle <- network.cycles[[c]]
              
              named.cycle <- c()  
              
              for (n in 1:length(cycle)) {
                node.num <- cycle[n]
                node.info <- nodes[as.numeric(node.num),]
                named.node <- node.info$node_name
                
                named.cycle <- c(named.cycle, named.node)
              }
              
              if (c == 1) {
                named.network.cycles <-list(named.cycle)
              }
              if (c > 1) {
                named.network.cycles <- append(named.network.cycles, list(named.cycle))
              }
            }
            
            #Edges causing cycles:
            
            for (c in 1:length(named.network.cycles)) {
              named.cycle <- named.network.cycles[[c]]
              
              first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
              last.edges <- filtered.edges.df %>% filter(from == named.cycle[length(named.cycle)-1] & to == named.cycle[length(named.cycle)])
              
              if (c ==1) {
                cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                  first.edge.from = unique(first.edges$from.named),
                                                  first.edge.to = unique(first.edges$to.named),
                                                  first.edges.num= nrow(first.edges),
                                                  first.edges.database = unique(first.edges$database_source),
                                                  last.edge.from = unique(last.edges$from.named),
                                                  last.edge.to = unique(last.edges$to.named),
                                                  last.edges.num = nrow(last.edges),
                                                  last.edges.database = unique(last.edges$database_source))
              }
              
              if (c > 1) {
                cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                               first.edge.from = unique(first.edges$from.named),
                                               first.edge.to = unique(first.edges$to.named),
                                               first.edges.num= nrow(first.edges),
                                               first.edges.database = unique(first.edges$database_source),
                                               last.edge.from = unique(last.edges$from.named),
                                               last.edge.to = unique(last.edges$to.named),
                                               last.edges.num = nrow(last.edges),
                                               last.edges.database = unique(last.edges$database_source))
                cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
              }
              
            }
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(first.edge.from, "__", first.edge.to, sep = ""),
                                          last.edge = paste(last.edge.from, "__", last.edge.to, sep = ""))
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, cycle_id = paste(first.edge, ".....",last.edge, sep = ""))
            
            #Cycles seem to be due to weak nichenet connections at last edge:
            cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$last.edge)),]
            cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("last.edge.from", "last.edge.to"))]
            
            if (nrow(cyclic_edges_to_prune) >0 ) {
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$last.edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$last.edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            
            #Re-make Graph without above pruned cycle edges:    
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, layout)
            
            ## Above function doesn't find all cycles though! Use igraph to find feedback edges requiring pruning for DAG creation:
            merged.igraph <- as.igraph(merged.graph)
            cycles <- feedback_arc_set(merged.igraph)
            cycle.edge.ids <- as_ids(cycles)
            
            if (length(cycle.edge.ids)>0){
              
              #convert to list:
              for (e in 1:length(cycle.edge.ids)) {
                if (e ==1){
                  cyclic.edges <- list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]])))
                  
                }
                
                if (e>1){
                  cyclic.edges <-  append(cyclic.edges, list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]]))))
                }
              }   
              
              #Add names:
              for (c in 1:length(cyclic.edges)) {
                cycle <- cyclic.edges[[c]]
                
                named.cycle <- c()  
                
                for (n in 1:length(cycle)) {
                  node.num <- cycle[n]
                  node.info <- nodes[as.numeric(node.num),]
                  named.node <- node.info$node_name
                  
                  named.cycle <- c(named.cycle, named.node)
                }
                
                if (c == 1) {
                  named.network.cycles <-list(named.cycle)
                }
                if (c > 1) {
                  named.network.cycles <- append(named.network.cycles, list(named.cycle))
                }
              }
              
              for (c in 1:length(named.network.cycles)) {
                named.cycle <- named.network.cycles[[c]]
                
                first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
                
                if (c ==1) {
                  cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                    edge.from = unique(first.edges$from.named),
                                                    edge.to = unique(first.edges$to.named),
                                                    edges.num= nrow(first.edges),
                                                    edges.database = unique(first.edges$database_source))
                  
                }
                
                if (c > 1) {
                  cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                 edge.from = unique(first.edges$from.named),
                                                 edge.to = unique(first.edges$to.named),
                                                 edges.num= nrow(first.edges),
                                                 edges.database = unique(first.edges$database_source))
                  
                  cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                }
                
              }
              
              cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(edge.from, "__", edge.to, sep = ""))
              
              #All cycles look to be due to weak nichenet connections at last edge:
              cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$first.edge)),]
              cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("edge.from", "edge.to"))]
              
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            # Make Final Graph:    
            message("Making final graph. . . ")
            
            #add shapes to nodes if you want:
            if (!is.null(node_shape)) {
              
              node.shape.df <- DE.table[,which(colnames(DE.table) %in% c("CellType", node_shape))]  
              node.shape.df <- node.shape.df[-which(duplicated(node.shape.df$CellType)),]
              
              filtered.nodes.df <- left_join(filtered.nodes.df, node.shape.df)
              
              if (class(seurat.obj@meta.data[,node_shape]) == "factor") {
                filtered.nodes.df[node_shape] <- factor(filtered.nodes.df[,node_shape], 
                                                        levels = levels(seurat.obj@meta.data[,node_shape]))
              }
            }
            
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType", node_shape, "node_label"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            ## Ligand Edge Label mapping:
            label.mapping.df <- filtered.edges.df
            label.mapping.df$label_x <- NA
            label.mapping.df$label_y <- NA
            label.mapping.df$Edge_label <- ""
            
            for (l in 1:nrow(label.mapping.df)) {
              
              lig.from <- label.mapping.df$from.named[l]
              lig.to <- label.mapping.df$to.named[l]
              
              node1_info <- layout %>% filter(node_name == lig.from)
              node2_info <- layout %>% filter(node_name == lig.to)
              
              label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
              label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
            }
            
            #Only label ligands that are reasonably far away from each other, with DE prioritization:
            max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
            
            for (l in 1:length(unique(label.mapping.df$Ligand))) {
              lig <- unique(label.mapping.df$Ligand)[l]
              
              lig.edges <- label.mapping.df %>% filter(Ligand == lig)
              
              lig.edges$Edge_label <- "temp"
              
              while(any(lig.edges$Edge_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.edges <- grouped.lig.edges %>% arrange(desc(DE),centroid_dist) %>% dplyr::pull(edgeID)
                Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
                
                lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                
                lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                
              }
              
              label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
              
            }
            
            ## Node Labels:
            node.labels.df <- layout
            node.labels.df$node_label <- ""
            
            node.labels.df$node_gene <- sapply(node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            #Only label nodes that are reasonably far away from each other:
            
            for (g in 1:length(unique(node.labels.df$node_gene))) {
              gene <- unique(node.labels.df$node_gene)[g]
              
              gene.nodes <- node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.035) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                
              }
              
              node.labels.df[which(node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, node.labels.df)
            
            
            ### Fold change comparisons for merged graph comparisons:
            filtered.edges.df$from_to_joint <- paste(filtered.edges.df$from.named, filtered.edges.df$to.named, sep = "__")
            
            filtered.edges.df <- filtered.edges.df[,-which(colnames(filtered.edges.df) %in% c("edgeID", "DE"))]
            filtered.edges.df <- filtered.edges.df %>% pivot_wider(values_from = "Scaled_Connection_Product", names_from = "Lesion")
            filtered.edges.df[is.na(filtered.edges.df)] <- 0
            
            lesions <- as.character(unique(filtered.nodes.df$Lesion))
            
            filtered.edges.df <- filtered.edges.df %>% mutate(fold_change = (eval(as.symbol(lesions[which(lesions == condition_of_interest)]))+0.1)/(eval(as.symbol(lesions[-which(lesions == condition_of_interest)]))+0.1))

            filtered.edges.df <- filtered.edges.df %>% mutate(foldchange_category = case_when(log2(fold_change) < -2 ~ "Down",
                                                                                              log2(fold_change) >= -2 & log2(fold_change) <= 2  ~ "Similar",
                                                                                              log2(fold_change) > 2 ~ "Up"))
            
            
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
          }
          
          
          #Set point colors for plotting:
          GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
          
          ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
          node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
          
          if (plot_layout == "stress") {
            
            ##Plot With Facets:
            g <- ggraph(RL_graph, layout = "manual",
                        x = layout_full[, "x"],
                        y = layout_full[, "y"]) + 
              geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
              geom_node_point(aes(color = CellType), size = 0.5) + 
              geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
              scale_edge_width(range = c(0.1, 0.5)) +
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              facet_nodes(~Lesion) +
              scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
              new_scale_color()+
              scale_color_manual(values = c("#000000", "#750000")) +
              geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0)
            
            network.output.result <- list("Graph_object" = list(RL_graph),
                                          "Layout" = layout_full,
                                          "Plot" = g)
            
          }
          
          if (plot_layout == "sugiyama") {
            #Sugiyama Layout:
            sugi_layout <- layout_with_sugiyama(merged.graph, hgap = 500)
            sugi_layout <- as.data.frame(sugi_layout$layout)
            colnames(sugi_layout) <- c("x", "y")
            
            sugi_layout <- cbind(nodes, sugi_layout)
            
            ## Sugi Node Labels:
            sugi.node.labels.df <- sugi_layout
            sugi.node.labels.df$node_label <- ""
            
            sugi.node.labels.df$node_gene <- sapply(sugi.node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            max.sugi.layout.dist <- ((max(sugi_layout$x)-min(sugi_layout$x))^2+(max(sugi_layout$y)-min(sugi_layout$y))^2)^(0.5)
            
            #Only label nodes that are reasonably far away from each other:
            for (g in 1:length(unique(sugi.node.labels.df$node_gene))) {
              gene <- unique(sugi.node.labels.df$node_gene)[g]
              
              gene.nodes <- sugi.node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.sugi.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                
              }
              
              sugi.node.labels.df[which(sugi.node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            g <- ggraph(merged.graph, layout = "manual",
                        x = sugi_layout[, "x"],
                        y = sugi_layout[, "y"]) + 
              geom_node_point(data = sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
              geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
              scale_shape_manual(values = node_shape_vals)+
              scale_color_manual(values = node.colors1) +
              geom_edge_diagonal(aes(color = foldchange_category), width = 0.5, alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              scale_edge_colour_manual(values = c("#0F1BDF", "#787878","#DA1B1B")) +
              expand_limits(y = c(min(sugi.node.labels.df$y), max(sugi.node.labels.df$y)+0.05*(range(sugi.node.labels.df$y)[2]-range(sugi.node.labels.df$y)[1]))) +
              geom_text_repel(data = sugi.node.labels.df, aes(x = x, y = y, label = node_label), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5)+
              labs(title = paste(condition_of_interest, " vs ", str_c(as.character(unique(filtered.nodes.df$Lesion))[-which(as.character(unique(filtered.nodes.df$Lesion)) == condition_of_interest)],collapse = ", "), sep = ""))+
              theme(plot.title = element_text(hjust = 0.5))
            
            network.output.result <- list("Graph_object" = list(merged.graph),
                                          "Layout" = sugi_layout,
                                          "Plot" = g)
          }
          
          plot(g)
          return(network.output.result)
        }
        
        ##Edited to highlight DE genes in red:
        plot_lig_rec_merged_network <- function(network, DE.table, seurat.obj, condition_of_interest,
                                                facet_by = "Lesion", node_shape = "Basic_Celltype", node_shape_vals = c(15, 8, 13, 16), 
                                                create_DAG = TRUE, plot_layout = "sugiyama") {
          
          # Ligand-Receptor-Celltype Specific Node Plot:        
          message("Making initial ligand-receptor celltype-specific network. . . ")
          #Create node df:
          lesion.list <- unique(network[, facet_by])                                                 
          
          node_list <- c()
          
          for (i in 1:length(lesion.list)) {
            
            #Filter Edges:
            lesion.edges <- filter(network, eval(as.symbol(facet_by)) == lesion.list[i])                     
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet_by = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          #Merge Nodes:
          nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            nodes.df[,facet_by] <- factor(nodes.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
          }
          
          #Create Edges.df:
          edges.df <- network                                                                          
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            edges.df[,facet_by] <- factor(edges.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          edges.df <- edges.df[,which(colnames(edges.df) %in% c("from", "to", facet_by, "Ligand_Source", "Ligand", "Receptor_source", "Receptor","log10_Connection_product", "database_source"))]
          
          #Remove duplicate nichenet edges: (not sure why this happened. . . need to fix this eventually)
          edges.df <- edges.df %>% mutate(edgeID = paste(from, "__", to, "__", eval(as.symbol(facet_by)), sep = ""))
          
          duplicates <-  edges.df$edgeID[which(duplicated(edges.df$edgeID))]         
          duplicates <- edges.df[which(edges.df$edgeID %in% duplicates),]
          duplicate.max.connections <- duplicates %>% group_by(edgeID) %>% summarize(max.connection = max(log10_Connection_product))
          
          duplicates <- duplicates[-which(duplicated(duplicates$edgeID)),]
          duplicates <- full_join(duplicates, duplicate.max.connections)
          duplicates <- duplicates[,-which(colnames(duplicates) == "log10_Connection_product")]
          
          duplicates <- duplicates %>% rename(log10_Connection_product = "max.connection")
          
          edges.df <- edges.df[-which(edges.df$edgeID %in% duplicates$edgeID),]
          edges.df <- bind_rows(edges.df, duplicates)
          
          # Scale edge width by percent of max connection strength:
          message("Scaling edges. . . ")
          
          ### Why is Day4 allergy Myeloid:CXCL9-KRT:CXCR3 almost as strong as Day4 Allergy lymphocyte CXCR3 connections?
          
          max.expression.vals <- edges.df %>% group_by(Ligand, Receptor) %>% summarize(max.expression = max(log10_Connection_product))
          
          scaled_edges <- edges.df
          
          for (r in 1:nrow(max.expression.vals)) {
            
            connection.row <- max.expression.vals[r,]
            connection.max <- max.expression.vals[r,"max.expression"]
            
            connection.edges <- edges.df %>% filter(Ligand == connection.row$Ligand & Receptor == connection.row$Receptor)
            connection.edges$log10_Connection_product <- sapply(connection.edges$log10_Connection_product, FUN = function(x) as.numeric(round(10^(x)/(10^(connection.max)), digits = 3)))
            
            scaled_edges[which(scaled_edges$edgeID %in% connection.edges$edgeID),] <- connection.edges
          }
          
          edges.df <- scaled_edges
          colnames(edges.df)[which(colnames(edges.df) == "log10_Connection_product")] <- "Scaled_Connection_Product"
          
          
          #Add edge DE status:
          edges.df$DE <- FALSE
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            les.celltypes <- unique(lesion.edges$Ligand_Source)
            
            for (lc in 1:length(les.celltypes)) {
              
              les.cell <- les.celltypes[lc]
              
              les.cell.edges <- lesion.edges %>% filter(Ligand_Source == les.cell)
              
              les.cell.DE <- DE.table %>% filter(Comparison == les & CellType == les.cell & Color_factor == "Significantly Up") %>% dplyr::pull(Gene)
              
              DE.complex.genes <- deconvoluted %>% filter(gene_name %in% les.cell.DE) %>% dplyr::pull(complex_name) %>% unique()
              
              #Only label DE CellphoneDB
              les.cell.edges.DE <- les.cell.edges %>% filter(Ligand %in% unique(les.cell.DE, DE.complex.genes) & database_source == "CellphoneDB") %>% dplyr::pull(edgeID)
              
              if (length(les.cell.edges.DE) >0) {
                edges.df[which(edges.df$edgeID %in% les.cell.edges.DE), "DE"] <- TRUE
              }
              
            }
          }
          
          
          #Prune edges to only interconnected CellphoneDB & Nichenet edges:
          message("Pruning edges to interconnected cellphonedb and nichenet edges. . . ")
          
          filtered.edges.df <- as.data.frame(edges.df)[0,]
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            if (nrow(lesion.edges)>0) {
              
              #Prune nichenet edges:
              lesion.nichenet.edges <- lesion.edges %>% filter(database_source == "Nichenet" & Scaled_Connection_Product > 0) 
              
              if (nrow(lesion.nichenet.edges)>0){
                
                nichenet.edges.to.prune <- c()
                
                for (ne in 1:length(unique(lesion.nichenet.edges$from))) {
                  
                  nichenet.edge.from <- unique(lesion.nichenet.edges$from)[ne]
                  
                  upstream.edges <- lesion.edges %>% filter(to == nichenet.edge.from & Scaled_Connection_Product > 0 & database_source == "CellphoneDB")
                  
                  if (nrow(upstream.edges) == 0) {
                    nichenet.edges.to.prune <- c(nichenet.edges.to.prune, nichenet.edge.from)
                  }
                  
                }    
                
                lesion.nichenet.edges <- lesion.nichenet.edges %>% filter(!from %in% nichenet.edges.to.prune)
              }
              
              #Prune CellphoneDB edges:
              lesion.cellphonedb.edges <- lesion.edges %>% filter(database_source == "CellphoneDB" & Scaled_Connection_Product > 0) 
              
              if(nrow(lesion.cellphonedb.edges) >0) {
                
                cellphonedb.edges.to.prune <- c()
                
                for (ce in 1:length(unique(lesion.cellphonedb.edges$to))) {
                  
                  cellphonedb.edge.to <- unique(lesion.cellphonedb.edges$to)[ce]
                  
                  downstream.edges <- lesion.edges %>% filter(from == cellphonedb.edge.to & Scaled_Connection_Product > 0 & database_source == "Nichenet")
                  
                  if (nrow(downstream.edges) == 0) {
                    cellphonedb.edges.to.prune <- c(cellphonedb.edges.to.prune, cellphonedb.edge.to)
                  }
                  
                }
                
                lesion.cellphonedb.edges <- lesion.cellphonedb.edges %>% filter(!to %in% cellphonedb.edges.to.prune)
              }
              
              lesion.edges <- bind_rows(lesion.cellphonedb.edges, lesion.nichenet.edges)         
              
              filtered.edges.df <- bind_rows(filtered.edges.df, lesion.edges)
            }
          }
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.edges.df[,facet_by] <- factor(filtered.edges.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          #Duplicate node names columns so that each is preserved after graph formation:
          filtered.edges.df$from.named <- filtered.edges.df$from
          filtered.edges.df$to.named <- filtered.edges.df$to
          
          
          #Create filtered.nodes.df:
          node_list <- c()
          
          for (i in 1:length(unique(filtered.edges.df[,facet_by]))) {
            
            les <- as.character(unique(filtered.edges.df[,facet_by]))[i]
            #Filter Edges:
            lesion.edges <- filter(filtered.edges.df, eval(as.symbol(facet_by)) == les)                   
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          
          #Merge Nodes:
          filtered.nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.nodes.df[,facet_by] <- factor(filtered.nodes.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          
          #Make Initial Graph:    
          RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          #Correct misconnected edges:
          message("Correcting tbl_graph misconnected edges. . . ")
          
          graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
          graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
          
          duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
          
          filtered.nodes.df <- rownames_to_column(filtered.nodes.df)
          
          for (r in 1:nrow(duplicates)) {
            
            les <- as.character(duplicates[r,facet_by])
            
            if (duplicates[r,"database_source"] == "CellphoneDB") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
            
            if (duplicates[r,"database_source"] == "Nichenet") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
          }
          
          graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
          
          RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
          
          
          # Define Layout using merged data from all lesions:
          nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
          nodes <- nodes[-which(duplicated(nodes$node_name)),]  
          
          merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          layout <- layout_with_stress(merged.graph)
          colnames(layout) <- c("x", "y")
          
          layout <- cbind(nodes, layout)
          
          ## Ligand Edge Label mapping:
          label.mapping.df <- filtered.edges.df
          label.mapping.df$label_x <- NA
          label.mapping.df$label_y <- NA
          label.mapping.df$Edge_label <- ""
          
          for (l in 1:nrow(label.mapping.df)) {
            
            lig.from <- label.mapping.df$from.named[l]
            lig.to <- label.mapping.df$to.named[l]
            
            node1_info <- layout %>% filter(node_name == lig.from)
            node2_info <- layout %>% filter(node_name == lig.to)
            
            label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
            label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
          }
          
          #Only label ligands that are reasonably far away from each other, with DE prioritization:
          max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
          
          for (l in 1:length(unique(label.mapping.df$Ligand))) {
            lig <- unique(label.mapping.df$Ligand)[l]
            
            lig.edges <- label.mapping.df %>% filter(Ligand == lig)
            
            lig.edges$Edge_label <- "temp"
            
            while(any(lig.edges$Edge_label =="temp")){
              
              #temporarily make a column for fraction of max distance from first ungrouped edge:
              grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
              
              grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.035) 
              
              #Label closest to centroid edges grouped by 5% of max layout distance:
              centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
              centroid_x <- sum(centroid_x)/length(centroid_x)
              centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
              centroid_y <- sum(centroid_y)/length(centroid_y)
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
              
              group.edges <- grouped.lig.edges %>% arrange(desc(DE), centroid_dist) %>% dplyr::pull(edgeID)
              Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
              
              lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
              
              lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
              
            }
            
            label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
            
          }
          
          #Copy coordinates to node df containing all faceted nodes:
          layout_full <- full_join(filtered.nodes.df, layout)
          
          if (create_DAG == TRUE)  {
            message("Pruning to directed acyclic graph. . . ")
            
            #Prune cycles and feedback loops to create a directed acyclic graph:
            
            ## Find Cycles in network:
            FindCycles = function(g) {
              Cycles = NULL
              for(v1 in V(g)) {
                if(degree(g, v1, mode="in") == 0) { next }
                GoodNeighbors = neighbors(g, v1, mode="out")
                GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
                for(v2 in GoodNeighbors) {
                  TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
                  TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
                  TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
                  Cycles  = c(Cycles, TempCyc)
                }
              }
              Cycles
            }
            merged.igraph <- as.igraph(merged.graph)
            network.cycles <- FindCycles(merged.igraph)
            
            nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
            
            for (c in 1:length(network.cycles)) {
              cycle <- network.cycles[[c]]
              
              named.cycle <- c()  
              
              for (n in 1:length(cycle)) {
                node.num <- cycle[n]
                node.info <- nodes[as.numeric(node.num),]
                named.node <- node.info$node_name
                
                named.cycle <- c(named.cycle, named.node)
              }
              
              if (c == 1) {
                named.network.cycles <-list(named.cycle)
              }
              if (c > 1) {
                named.network.cycles <- append(named.network.cycles, list(named.cycle))
              }
            }
            
            #Edges causing cycles:
            
            for (c in 1:length(named.network.cycles)) {
              named.cycle <- named.network.cycles[[c]]
              
              first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
              last.edges <- filtered.edges.df %>% filter(from == named.cycle[length(named.cycle)-1] & to == named.cycle[length(named.cycle)])
              
              if (c ==1) {
                cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                  first.edge.from = unique(first.edges$from.named),
                                                  first.edge.to = unique(first.edges$to.named),
                                                  first.edges.num= nrow(first.edges),
                                                  first.edges.database = unique(first.edges$database_source),
                                                  last.edge.from = unique(last.edges$from.named),
                                                  last.edge.to = unique(last.edges$to.named),
                                                  last.edges.num = nrow(last.edges),
                                                  last.edges.database = unique(last.edges$database_source))
              }
              
              if (c > 1) {
                cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                               first.edge.from = unique(first.edges$from.named),
                                               first.edge.to = unique(first.edges$to.named),
                                               first.edges.num= nrow(first.edges),
                                               first.edges.database = unique(first.edges$database_source),
                                               last.edge.from = unique(last.edges$from.named),
                                               last.edge.to = unique(last.edges$to.named),
                                               last.edges.num = nrow(last.edges),
                                               last.edges.database = unique(last.edges$database_source))
                cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
              }
              
            }
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(first.edge.from, "__", first.edge.to, sep = ""),
                                          last.edge = paste(last.edge.from, "__", last.edge.to, sep = ""))
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, cycle_id = paste(first.edge, ".....",last.edge, sep = ""))
            
            #Cycles seem to be due to weak nichenet connections at last edge:
            cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$last.edge)),]
            cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("last.edge.from", "last.edge.to"))]
            
            if (nrow(cyclic_edges_to_prune) >0 ) {
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$last.edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$last.edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            
            #Re-make Graph without above pruned cycle edges:    
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, layout)
            
            ## Above function doesn't find all cycles though! Use igraph to find feedback edges requiring pruning for DAG creation:
            merged.igraph <- as.igraph(merged.graph)
            cycles <- feedback_arc_set(merged.igraph)
            cycle.edge.ids <- as_ids(cycles)
            
            if (length(cycle.edge.ids)>0){
              
              #convert to list:
              for (e in 1:length(cycle.edge.ids)) {
                if (e ==1){
                  cyclic.edges <- list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]])))
                  
                }
                
                if (e>1){
                  cyclic.edges <-  append(cyclic.edges, list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]]))))
                }
              }   
              
              #Add names:
              for (c in 1:length(cyclic.edges)) {
                cycle <- cyclic.edges[[c]]
                
                named.cycle <- c()  
                
                for (n in 1:length(cycle)) {
                  node.num <- cycle[n]
                  node.info <- nodes[as.numeric(node.num),]
                  named.node <- node.info$node_name
                  
                  named.cycle <- c(named.cycle, named.node)
                }
                
                if (c == 1) {
                  named.network.cycles <-list(named.cycle)
                }
                if (c > 1) {
                  named.network.cycles <- append(named.network.cycles, list(named.cycle))
                }
              }
              
              for (c in 1:length(named.network.cycles)) {
                named.cycle <- named.network.cycles[[c]]
                
                first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
                
                if (c ==1) {
                  cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                    edge.from = unique(first.edges$from.named),
                                                    edge.to = unique(first.edges$to.named),
                                                    edges.num= nrow(first.edges),
                                                    edges.database = unique(first.edges$database_source))
                  
                }
                
                if (c > 1) {
                  cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                 edge.from = unique(first.edges$from.named),
                                                 edge.to = unique(first.edges$to.named),
                                                 edges.num= nrow(first.edges),
                                                 edges.database = unique(first.edges$database_source))
                  
                  cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                }
                
              }
              
              cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(edge.from, "__", edge.to, sep = ""))
              
              #All cycles look to be due to weak nichenet connections at last edge:
              cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$first.edge)),]
              cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("edge.from", "edge.to"))]
              
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            # Make Final Graph:    
            message("Making final graph. . . ")
            
            #add shapes to nodes if you want:
            if (!is.null(node_shape)) {
              
              node.shape.df <- DE.table[,which(colnames(DE.table) %in% c("CellType", node_shape))]  
              node.shape.df <- node.shape.df[-which(duplicated(node.shape.df$CellType)),]
              
              filtered.nodes.df <- left_join(filtered.nodes.df, node.shape.df)
              
              if (class(seurat.obj@meta.data[,node_shape]) == "factor") {
                filtered.nodes.df[node_shape] <- factor(filtered.nodes.df[,node_shape], 
                                                        levels = levels(seurat.obj@meta.data[,node_shape]))
              }
            }
            
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType", node_shape, "node_label"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            ## Ligand Edge Label mapping:
            label.mapping.df <- filtered.edges.df
            label.mapping.df$label_x <- NA
            label.mapping.df$label_y <- NA
            label.mapping.df$Edge_label <- ""
            
            for (l in 1:nrow(label.mapping.df)) {
              
              lig.from <- label.mapping.df$from.named[l]
              lig.to <- label.mapping.df$to.named[l]
              
              node1_info <- layout %>% filter(node_name == lig.from)
              node2_info <- layout %>% filter(node_name == lig.to)
              
              label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
              label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
            }
            
            #Only label ligands that are reasonably far away from each other, with DE prioritization:
            max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
            
            for (l in 1:length(unique(label.mapping.df$Ligand))) {
              lig <- unique(label.mapping.df$Ligand)[l]
              
              lig.edges <- label.mapping.df %>% filter(Ligand == lig)
              
              lig.edges$Edge_label <- "temp"
              
              while(any(lig.edges$Edge_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.edges <- grouped.lig.edges %>% arrange(desc(DE),centroid_dist) %>% dplyr::pull(edgeID)
                Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
                
                lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                
                lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                
              }
              
              label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
              
            }
            
            ## Node Labels:
            node.labels.df <- layout
            node.labels.df$node_label <- ""
            
            node.labels.df$node_gene <- sapply(node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            #Only label nodes that are reasonably far away from each other:
            
            for (g in 1:length(unique(node.labels.df$node_gene))) {
              gene <- unique(node.labels.df$node_gene)[g]
              
              gene.nodes <- node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.035) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                
              }
              
              node.labels.df[which(node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, node.labels.df)
            
            
            ### Fold change comparisons for merged graph comparisons:
            filtered.edges.df$from_to_joint <- paste(filtered.edges.df$from.named, filtered.edges.df$to.named, sep = "__")
            
            filtered.edges.df <- filtered.edges.df[,-which(colnames(filtered.edges.df) %in% c("edgeID", "DE"))]
            filtered.edges.df <- filtered.edges.df %>% pivot_wider(values_from = "Scaled_Connection_Product", names_from = "Lesion")
            filtered.edges.df[is.na(filtered.edges.df)] <- 0
            
            lesions <- as.character(unique(filtered.nodes.df$Lesion))
            
            filtered.edges.df <- filtered.edges.df %>% mutate(fold_change = (eval(as.symbol(lesions[which(lesions == condition_of_interest)]))+0.1)/(eval(as.symbol(lesions[-which(lesions == condition_of_interest)]))+0.1))
            
            filtered.edges.df <- filtered.edges.df %>% mutate(foldchange_category = case_when(log2(fold_change) < -2 ~ "Down",
                                                                                              log2(fold_change) >= -2 & log2(fold_change) <= 2  ~ "Similar",
                                                                                              log2(fold_change) > 2 ~ "Up"))
            
            
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
          }
          
          
          #Set point colors for plotting:
          GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
          
          ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
          node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
          
          if (plot_layout == "stress") {
            
            ##Plot With Facets:
            g <- ggraph(RL_graph, layout = "manual",
                        x = layout_full[, "x"],
                        y = layout_full[, "y"]) + 
              geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
              geom_node_point(aes(color = CellType), size = 0.5) + 
              geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
              scale_edge_width(range = c(0.1, 0.5)) +
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              facet_nodes(~Lesion) +
              scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
              new_scale_color()+
              scale_color_manual(values = c("#000000", "#750000")) +
              geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0)
            
            network.output.result <- list("Graph_object" = list(RL_graph),
                                          "Layout" = layout_full,
                                          "Plot" = g)
            
          }
          
          if (plot_layout == "sugiyama") {
            #Sugiyama Layout:
            sugi_layout <- layout_with_sugiyama(merged.graph, hgap = 500)
            sugi_layout <- as.data.frame(sugi_layout$layout)
            colnames(sugi_layout) <- c("x", "y")
            
            sugi_layout <- cbind(nodes, sugi_layout)
            
            ## Sugi Node Labels:
            sugi.node.labels.df <- sugi_layout
            sugi.node.labels.df$node_label <- ""
            
            sugi.node.labels.df$node_gene <- sapply(sugi.node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            max.sugi.layout.dist <- ((max(sugi_layout$x)-min(sugi_layout$x))^2+(max(sugi_layout$y)-min(sugi_layout$y))^2)^(0.5)
            
            
            #Color Significantly Up DE Genes:
            DE.genes.to.color <- day2_DE.table %>% filter(Comparison == condition_of_interest & Color_factor == "Significantly Up" & Gene %in% unique(sugi.node.labels.df$node_gene))
            
            DE.genes.to.color %>% filter(Gene == "IFNG")
            
            sugi.node.labels.df$node_label_color <- "black"
            
            for (r in 1:nrow(sugi.node.labels.df)) {
              
              cell <- sugi.node.labels.df$CellType[r]
              gene <- sugi.node.labels.df$node_gene[r]
              
              filtered.gene.result <- DE.genes.to.color %>% filter(Gene == gene & CellType == cell)
              
              if (nrow(filtered.gene.result) >0) {
                sugi.node.labels.df$node_label_color[r]<- "red"
              }
              
            }
            
            #Label nodes that are reasonably far away from each other, with preference for DE genes that will be labeled in red:
            
            for (g in 1:length(unique(sugi.node.labels.df$node_gene))) {
              gene <- unique(sugi.node.labels.df$node_gene)[g]
              
              gene.nodes <- sugi.node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.sugi.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(desc(node_label_color),centroid_dist) %>% dplyr::pull(node_name)
                node.to.label <- group.nodes[1]
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% node.to.label),"node_label"] <- gene
                
              }
              
              sugi.node.labels.df[which(sugi.node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            g <- ggraph(merged.graph, layout = "manual",
                        x = sugi_layout[, "x"],
                        y = sugi_layout[, "y"]) + 
              geom_node_point(data = sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
              geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
              scale_shape_manual(values = node_shape_vals)+
              scale_color_manual(values = node.colors1) +
              geom_edge_diagonal(aes(color = foldchange_category), width = 0.5, alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              scale_edge_colour_manual(values = c("#0F1BDF", "#787878","#DA1B1B")) +
              expand_limits(y = c(min(sugi.node.labels.df$y), max(sugi.node.labels.df$y)+0.05*(range(sugi.node.labels.df$y)[2]-range(sugi.node.labels.df$y)[1]))) +
              ggnewscale::new_scale_color() +
              geom_text_repel(data = sugi.node.labels.df, aes(x = x, y = y, label = node_label, color = node_label_color), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5) +
              scale_color_manual(values = c("black", "red")) +
              labs(title = paste(condition_of_interest, " vs ", str_c(as.character(unique(filtered.nodes.df$Lesion))[-which(as.character(unique(filtered.nodes.df$Lesion)) == condition_of_interest)],collapse = ", "), sep = ""))+
              theme(plot.title = element_text(hjust = 0.5))
            
            network.output.result <- list("Graph_object" = list(merged.graph),
                                          "Layout" = sugi_layout,
                                          "Plot" = g)
          }
          
          plot(g)
          return(network.output.result)
        }
        
        ##Edited to highlight both DE and Trending Ligands:
        plot_lig_rec_merged_network <- function(network, DE.table, seurat.obj, condition_of_interest, trend_ligands_to_label = NULL,
                                                facet_by = "Lesion", node_shape = "Basic_Celltype", node_shape_vals = c(15, 8, 13, 16), 
                                                create_DAG = TRUE, plot_layout = "sugiyama") {
          
          # Ligand-Receptor-Celltype Specific Node Plot:        
          message("Making initial ligand-receptor celltype-specific network. . . ")
          #Create node df:
          lesion.list <- unique(network[, facet_by])                                                 
          
          node_list <- c()
          
          for (i in 1:length(lesion.list)) {
            
            #Filter Edges:
            lesion.edges <- filter(network, eval(as.symbol(facet_by)) == lesion.list[i])                     
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet_by = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          #Merge Nodes:
          nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            nodes.df[,facet_by] <- factor(nodes.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
          }
          
          #Create Edges.df:
          edges.df <- network                                                                          
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            edges.df[,facet_by] <- factor(edges.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          edges.df <- edges.df[,which(colnames(edges.df) %in% c("from", "to", facet_by, "Ligand_Source", "Ligand", "Receptor_source", "Receptor","log10_Connection_product", "database_source"))]
          
          #Remove duplicate nichenet edges: (not sure why this happened. . . need to fix this eventually)
          edges.df <- edges.df %>% mutate(edgeID = paste(from, "__", to, "__", eval(as.symbol(facet_by)), sep = ""))
          
          duplicates <-  edges.df$edgeID[which(duplicated(edges.df$edgeID))]         
          duplicates <- edges.df[which(edges.df$edgeID %in% duplicates),]
          duplicate.max.connections <- duplicates %>% group_by(edgeID) %>% summarize(max.connection = max(log10_Connection_product))
          
          duplicates <- duplicates[-which(duplicated(duplicates$edgeID)),]
          duplicates <- full_join(duplicates, duplicate.max.connections)
          duplicates <- duplicates[,-which(colnames(duplicates) == "log10_Connection_product")]
          
          duplicates <- duplicates %>% rename(log10_Connection_product = "max.connection")
          
          edges.df <- edges.df[-which(edges.df$edgeID %in% duplicates$edgeID),]
          edges.df <- bind_rows(edges.df, duplicates)
          
          # Scale edge width by percent of max connection strength:
          message("Scaling edges. . . ")
          
          max.expression.vals <- edges.df %>% group_by(Ligand, Receptor) %>% summarize(max.expression = max(log10_Connection_product))
          
          scaled_edges <- edges.df
          
          for (r in 1:nrow(max.expression.vals)) {
            
            connection.row <- max.expression.vals[r,]
            connection.max <- max.expression.vals[r,"max.expression"]
            
            connection.edges <- edges.df %>% filter(Ligand == connection.row$Ligand & Receptor == connection.row$Receptor)
            connection.edges$log10_Connection_product <- sapply(connection.edges$log10_Connection_product, FUN = function(x) as.numeric(round(10^(x)/(10^(connection.max)), digits = 3)))
            
            scaled_edges[which(scaled_edges$edgeID %in% connection.edges$edgeID),] <- connection.edges
          }
          
          edges.df <- scaled_edges
          colnames(edges.df)[which(colnames(edges.df) == "log10_Connection_product")] <- "Scaled_Connection_Product"
          
          
          #Add edge DE status:
          edges.df$DE <- FALSE
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            les.celltypes <- unique(lesion.edges$Ligand_Source)
            
            for (lc in 1:length(les.celltypes)) {
              
              les.cell <- les.celltypes[lc]
              
              les.cell.edges <- lesion.edges %>% filter(Ligand_Source == les.cell)
              
              les.cell.DE <- DE.table %>% filter(Comparison == les & CellType == les.cell & Color_factor == "Significantly Up" & logFC >= 1 & Reference != "Acetone") %>% dplyr::pull(Gene)
              
              DE.complex.genes <- deconvoluted %>% filter(gene_name %in% les.cell.DE) %>% dplyr::pull(complex_name) %>% unique()
              
              #Only label DE CellphoneDB
              les.cell.edges.DE <- les.cell.edges %>% filter(Ligand %in% unique(les.cell.DE, DE.complex.genes) & database_source == "CellphoneDB") %>% dplyr::pull(edgeID)
              
              if (length(les.cell.edges.DE) >0) {
                edges.df[which(edges.df$edgeID %in% les.cell.edges.DE), "DE"] <- TRUE
              }
              
            }
          }
          
          
          #Prune edges to only interconnected CellphoneDB & Nichenet edges:
          message("Pruning edges to interconnected cellphonedb and nichenet edges. . . ")
          
          filtered.edges.df <- as.data.frame(edges.df)[0,]
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            if (nrow(lesion.edges)>0) {
              
              #Prune nichenet edges:
              lesion.nichenet.edges <- lesion.edges %>% filter(database_source == "Nichenet" & Scaled_Connection_Product > 0) 
              
              if (nrow(lesion.nichenet.edges)>0){
                
                nichenet.edges.to.prune <- c()
                
                for (ne in 1:length(unique(lesion.nichenet.edges$from))) {
                  
                  nichenet.edge.from <- unique(lesion.nichenet.edges$from)[ne]
                  
                  upstream.edges <- lesion.edges %>% filter(to == nichenet.edge.from & Scaled_Connection_Product > 0 & database_source == "CellphoneDB")
                  
                  if (nrow(upstream.edges) == 0) {
                    nichenet.edges.to.prune <- c(nichenet.edges.to.prune, nichenet.edge.from)
                  }
                  
                }    
                
                lesion.nichenet.edges <- lesion.nichenet.edges %>% filter(!from %in% nichenet.edges.to.prune)
              }
              
              #Prune CellphoneDB edges:
              lesion.cellphonedb.edges <- lesion.edges %>% filter(database_source == "CellphoneDB" & Scaled_Connection_Product > 0) 
              
              if(nrow(lesion.cellphonedb.edges) >0) {
                
                cellphonedb.edges.to.prune <- c()
                
                for (ce in 1:length(unique(lesion.cellphonedb.edges$to))) {
                  
                  cellphonedb.edge.to <- unique(lesion.cellphonedb.edges$to)[ce]
                  
                  downstream.edges <- lesion.edges %>% filter(from == cellphonedb.edge.to & Scaled_Connection_Product > 0 & database_source == "Nichenet")
                  
                  if (nrow(downstream.edges) == 0) {
                    cellphonedb.edges.to.prune <- c(cellphonedb.edges.to.prune, cellphonedb.edge.to)
                  }
                  
                }
                
                lesion.cellphonedb.edges <- lesion.cellphonedb.edges %>% filter(!to %in% cellphonedb.edges.to.prune)
              }
              
              lesion.edges <- bind_rows(lesion.cellphonedb.edges, lesion.nichenet.edges)         
              
              filtered.edges.df <- bind_rows(filtered.edges.df, lesion.edges)
            }
          }
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.edges.df[,facet_by] <- factor(filtered.edges.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          #Duplicate node names columns so that each is preserved after graph formation:
          filtered.edges.df$from.named <- filtered.edges.df$from
          filtered.edges.df$to.named <- filtered.edges.df$to
          
          #Create filtered.nodes.df:
          node_list <- c()
          
          for (i in 1:length(unique(filtered.edges.df[,facet_by]))) {
            
            les <- as.character(unique(filtered.edges.df[,facet_by]))[i]
            #Filter Edges:
            lesion.edges <- filter(filtered.edges.df, eval(as.symbol(facet_by)) == les)                   
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          
          #Merge Nodes:
          filtered.nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.nodes.df[,facet_by] <- factor(filtered.nodes.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          
          #Make Initial Graph:    
          RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          #Correct misconnected edges:
          message("Correcting tbl_graph misconnected edges. . . ")
          
          graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
          graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
          
          duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
          
          filtered.nodes.df <- rownames_to_column(filtered.nodes.df)
          
          for (r in 1:nrow(duplicates)) {
            
            les <- as.character(duplicates[r,facet_by])
            
            if (duplicates[r,"database_source"] == "CellphoneDB") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
            
            if (duplicates[r,"database_source"] == "Nichenet") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
          }
          
          graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
          
          RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
          
          
          # Define Layout using merged data from all lesions:
          nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
          nodes <- nodes[-which(duplicated(nodes$node_name)),]  
          
          merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          layout <- layout_with_stress(merged.graph)
          colnames(layout) <- c("x", "y")
          
          layout <- cbind(nodes, layout)
          
          ## Ligand Edge Label mapping:
          label.mapping.df <- filtered.edges.df
          label.mapping.df$label_x <- NA
          label.mapping.df$label_y <- NA
          label.mapping.df$Edge_label <- ""
          
          for (l in 1:nrow(label.mapping.df)) {
            
            lig.from <- label.mapping.df$from.named[l]
            lig.to <- label.mapping.df$to.named[l]
            
            node1_info <- layout %>% filter(node_name == lig.from)
            node2_info <- layout %>% filter(node_name == lig.to)
            
            label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
            label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
          }
          
          #Only label ligands that are reasonably far away from each other, with DE prioritization:
          max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
          
          for (l in 1:length(unique(label.mapping.df$Ligand))) {
            lig <- unique(label.mapping.df$Ligand)[l]
            
            lig.edges <- label.mapping.df %>% filter(Ligand == lig)
            
            lig.edges$Edge_label <- "temp"
            
            while(any(lig.edges$Edge_label =="temp")){
              
              #temporarily make a column for fraction of max distance from first ungrouped edge:
              grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
              
              grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.035) 
              
              #Label closest to centroid edges grouped by 5% of max layout distance:
              centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
              centroid_x <- sum(centroid_x)/length(centroid_x)
              centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
              centroid_y <- sum(centroid_y)/length(centroid_y)
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
              
              group.edges <- grouped.lig.edges %>% arrange(desc(DE), centroid_dist) %>% dplyr::pull(edgeID)
              Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
              
              lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
              
              lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
              
            }
            
            label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
            
          }
          
          #Copy coordinates to node df containing all faceted nodes:
          layout_full <- full_join(filtered.nodes.df, layout)
          
          if (create_DAG == TRUE)  {
            message("Pruning to directed acyclic graph. . . ")
            
            #Prune cycles and feedback loops to create a directed acyclic graph:
            
            ## Find Cycles in network:
            FindCycles = function(g) {
              Cycles = NULL
              for(v1 in V(g)) {
                if(degree(g, v1, mode="in") == 0) { next }
                GoodNeighbors = neighbors(g, v1, mode="out")
                GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
                for(v2 in GoodNeighbors) {
                  TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
                  TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
                  TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
                  Cycles  = c(Cycles, TempCyc)
                }
              }
              Cycles
            }
            merged.igraph <- as.igraph(merged.graph)
            network.cycles <- FindCycles(merged.igraph)
            
            nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
            
            for (c in 1:length(network.cycles)) {
              cycle <- network.cycles[[c]]
              
              named.cycle <- c()  
              
              for (n in 1:length(cycle)) {
                node.num <- cycle[n]
                node.info <- nodes[as.numeric(node.num),]
                named.node <- node.info$node_name
                
                named.cycle <- c(named.cycle, named.node)
              }
              
              if (c == 1) {
                named.network.cycles <-list(named.cycle)
              }
              if (c > 1) {
                named.network.cycles <- append(named.network.cycles, list(named.cycle))
              }
            }
            
            #Edges causing cycles:
            
            for (c in 1:length(named.network.cycles)) {
              named.cycle <- named.network.cycles[[c]]
              
              first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
              last.edges <- filtered.edges.df %>% filter(from == named.cycle[length(named.cycle)-1] & to == named.cycle[length(named.cycle)])
              
              if (c ==1) {
                cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                  first.edge.from = unique(first.edges$from.named),
                                                  first.edge.to = unique(first.edges$to.named),
                                                  first.edges.num= nrow(first.edges),
                                                  first.edges.database = unique(first.edges$database_source),
                                                  last.edge.from = unique(last.edges$from.named),
                                                  last.edge.to = unique(last.edges$to.named),
                                                  last.edges.num = nrow(last.edges),
                                                  last.edges.database = unique(last.edges$database_source))
              }
              
              if (c > 1) {
                cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                               first.edge.from = unique(first.edges$from.named),
                                               first.edge.to = unique(first.edges$to.named),
                                               first.edges.num= nrow(first.edges),
                                               first.edges.database = unique(first.edges$database_source),
                                               last.edge.from = unique(last.edges$from.named),
                                               last.edge.to = unique(last.edges$to.named),
                                               last.edges.num = nrow(last.edges),
                                               last.edges.database = unique(last.edges$database_source))
                cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
              }
              
            }
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(first.edge.from, "__", first.edge.to, sep = ""),
                                          last.edge = paste(last.edge.from, "__", last.edge.to, sep = ""))
            cyclic.edge.info.df <- mutate(cyclic.edge.info.df, cycle_id = paste(first.edge, ".....",last.edge, sep = ""))
            
            #Cycles seem to be due to weak nichenet connections at last edge:
            cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$last.edge)),]
            cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("last.edge.from", "last.edge.to"))]
            
            if (nrow(cyclic_edges_to_prune) >0 ) {
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$last.edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$last.edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            
            #Re-make Graph without above pruned cycle edges:    
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, layout)
            
            ## Above function doesn't find all cycles though! Use igraph to find feedback edges requiring pruning for DAG creation:
            merged.igraph <- as.igraph(merged.graph)
            cycles <- feedback_arc_set(merged.igraph)
            cycle.edge.ids <- as_ids(cycles)
            
            if (length(cycle.edge.ids)>0){
              
              #convert to list:
              for (e in 1:length(cycle.edge.ids)) {
                if (e ==1){
                  cyclic.edges <- list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]])))
                  
                }
                
                if (e>1){
                  cyclic.edges <-  append(cyclic.edges, list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]]))))
                }
              }   
              
              #Add names:
              for (c in 1:length(cyclic.edges)) {
                cycle <- cyclic.edges[[c]]
                
                named.cycle <- c()  
                
                for (n in 1:length(cycle)) {
                  node.num <- cycle[n]
                  node.info <- nodes[as.numeric(node.num),]
                  named.node <- node.info$node_name
                  
                  named.cycle <- c(named.cycle, named.node)
                }
                
                if (c == 1) {
                  named.network.cycles <-list(named.cycle)
                }
                if (c > 1) {
                  named.network.cycles <- append(named.network.cycles, list(named.cycle))
                }
              }
              
              for (c in 1:length(named.network.cycles)) {
                named.cycle <- named.network.cycles[[c]]
                
                first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
                
                if (c ==1) {
                  cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                    edge.from = unique(first.edges$from.named),
                                                    edge.to = unique(first.edges$to.named),
                                                    edges.num= nrow(first.edges),
                                                    edges.database = unique(first.edges$database_source))
                  
                }
                
                if (c > 1) {
                  cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                 edge.from = unique(first.edges$from.named),
                                                 edge.to = unique(first.edges$to.named),
                                                 edges.num= nrow(first.edges),
                                                 edges.database = unique(first.edges$database_source))
                  
                  cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                }
                
              }
              
              cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(edge.from, "__", edge.to, sep = ""))
              
              #All cycles look to be due to weak nichenet connections at last edge:
              cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$first.edge)),]
              cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("edge.from", "edge.to"))]
              
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            # Make Final Graph:    
            message("Making final graph. . . ")
            
            #add shapes to nodes if you want:
            if (!is.null(node_shape)) {
              
              node.shape.df <- DE.table[,which(colnames(DE.table) %in% c("CellType", node_shape))]  
              node.shape.df <- node.shape.df[-which(duplicated(node.shape.df$CellType)),]
              
              filtered.nodes.df <- left_join(filtered.nodes.df, node.shape.df)
              
              if (class(seurat.obj@meta.data[,node_shape]) == "factor") {
                filtered.nodes.df[node_shape] <- factor(filtered.nodes.df[,node_shape], 
                                                        levels = levels(seurat.obj@meta.data[,node_shape]))
              }
            }
            
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType", node_shape, "node_label"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            ## Ligand Edge Label mapping:
            label.mapping.df <- filtered.edges.df
            label.mapping.df$label_x <- NA
            label.mapping.df$label_y <- NA
            label.mapping.df$Edge_label <- ""
            
            for (l in 1:nrow(label.mapping.df)) {
              
              lig.from <- label.mapping.df$from.named[l]
              lig.to <- label.mapping.df$to.named[l]
              
              node1_info <- layout %>% filter(node_name == lig.from)
              node2_info <- layout %>% filter(node_name == lig.to)
              
              label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
              label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
            }
            
            #Only label ligands that are reasonably far away from each other, with DE prioritization:
            max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
            
            for (l in 1:length(unique(label.mapping.df$Ligand))) {
              lig <- unique(label.mapping.df$Ligand)[l]
              
              lig.edges <- label.mapping.df %>% filter(Ligand == lig)
              
              lig.edges$Edge_label <- "temp"
              
              while(any(lig.edges$Edge_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.edges <- grouped.lig.edges %>% arrange(desc(DE),centroid_dist) %>% dplyr::pull(edgeID)
                Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
                
                lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                
                lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                
              }
              
              label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
              
            }
            
            ## Node Labels:
            node.labels.df <- layout
            node.labels.df$node_label <- ""
            
            node.labels.df$node_gene <- sapply(node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            #Only label nodes that are reasonably far away from each other:
            
            for (g in 1:length(unique(node.labels.df$node_gene))) {
              gene <- unique(node.labels.df$node_gene)[g]
              
              gene.nodes <- node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.035) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                
              }
              
              node.labels.df[which(node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, node.labels.df)
            
            
            ### Fold change comparisons for merged graph comparisons:
            filtered.edges.df$from_to_joint <- paste(filtered.edges.df$from.named, filtered.edges.df$to.named, sep = "__")
            
            filtered.edges.df <- filtered.edges.df[,-which(colnames(filtered.edges.df) %in% c("edgeID", "DE"))]
            filtered.edges.df <- filtered.edges.df %>% pivot_wider(values_from = "Scaled_Connection_Product", names_from = "Lesion")
            filtered.edges.df[is.na(filtered.edges.df)] <- 0
            
            lesions <- as.character(unique(filtered.nodes.df$Lesion))
            
            filtered.edges.df <- filtered.edges.df %>% mutate(fold_change = (eval(as.symbol(lesions[which(lesions == condition_of_interest)]))+0.1)/(eval(as.symbol(lesions[-which(lesions == condition_of_interest)]))+0.1))
            
            filtered.edges.df <- filtered.edges.df %>% mutate(foldchange_category = case_when(log2(fold_change) < -2 ~ "Down",
                                                                                              log2(fold_change) >= -2 & log2(fold_change) <= 2  ~ "Similar",
                                                                                              log2(fold_change) > 2 ~ "Up"))
            
            
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
          }
          
          
          #Set point colors for plotting:
          GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
          
          ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
          node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
          
          if (plot_layout == "stress") {
            
            ##Plot With Facets:
            g <- ggraph(RL_graph, layout = "manual",
                        x = layout_full[, "x"],
                        y = layout_full[, "y"]) + 
              geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
              geom_node_point(aes(color = CellType), size = 0.5) + 
              geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
              scale_edge_width(range = c(0.1, 0.5)) +
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              facet_nodes(~Lesion) +
              scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
              new_scale_color()+
              scale_color_manual(values = c("#000000", "#750000")) +
              geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0)
            
            network.output.result <- list("Graph_object" = list(RL_graph),
                                          "Layout" = layout_full,
                                          "Plot" = g)
            
          }
          
          if (plot_layout == "sugiyama") {
            #Sugiyama Layout:
            sugi_layout <- layout_with_sugiyama(merged.graph, hgap = 500)
            sugi_layout <- as.data.frame(sugi_layout$layout)
            colnames(sugi_layout) <- c("x", "y")
            
            sugi_layout <- cbind(nodes, sugi_layout)
            
            ## Sugi Node Labels:
            sugi.node.labels.df <- sugi_layout
            sugi.node.labels.df$node_label <- ""
            
            sugi.node.labels.df$node_gene <- sapply(sugi.node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            max.sugi.layout.dist <- ((max(sugi_layout$x)-min(sugi_layout$x))^2+(max(sugi_layout$y)-min(sugi_layout$y))^2)^(0.5)
            
            
            #Color DE and Trending Genes in each condition:
            
            
            #Compare conditions to label DE and trending Ligands:
            #DE:
            DE.labeling.df <- edges.df %>% filter(DE == TRUE & database_source == "CellphoneDB")
            DE.labeling.df <- DE.labeling.df[,which(colnames(DE.labeling.df) %in% c("Ligand", "from", "Lesion", "edgeID"))]
            DE.labeling.df <- DE.labeling.df %>% pivot_wider(names_from = "Lesion", values_from = "edgeID") %>% as.data.frame() %>% suppressWarnings()
            
            labeling.data <- DE.labeling.df[,-which(colnames(DE.labeling.df) %in% c("Ligand", "from"))] 
            
            for (col in 1:ncol(labeling.data)) {
              
              labeling.data[,col] <- sapply(labeling.data[,col], FUN = function (x) {
                
                x <- str_c(unlist(x), collapse = "___")
                
                
                return(x)
              })
            }
            for (col in 1:ncol(labeling.data)) {
              
              labeling.data[,col] <- sapply(labeling.data[,col], FUN = function (x) {
                if (x != "") {
                  x <- "DE"
                }
                
                if (x == "") {
                  x <- "Non-DE"
                }
                
                
                return(x)
              })
            }
            DE.labeling.df[,-which(colnames(DE.labeling.df) %in% c("Ligand", "from"))]  <- labeling.data
            
            
            DE.labeling.df <- DE.labeling.df %>% mutate(DE_color_factor = case_when(eval(as.symbol(condition_of_interest)) == "DE" & eval(as.symbol(colnames(DE.labeling.df)[-which(colnames(DE.labeling.df) %in% c(condition_of_interest, "Ligand", "from"))])) == "Non-DE" ~ "Unique DE",
                                                                                    eval(as.symbol(condition_of_interest)) == "DE" & any(eval(as.symbol(colnames(DE.labeling.df)[-which(colnames(DE.labeling.df) %in% c(condition_of_interest, "Ligand", "from"))])) == "DE") ~ "Shared DE",
                                                                                    eval(as.symbol(condition_of_interest)) == "Non-DE" & any(eval(as.symbol(colnames(DE.labeling.df)[-which(colnames(DE.labeling.df) %in% c(condition_of_interest, "Ligand", "from"))])) == "DE") ~ "Comparison Unique DE",
                                                                                    TRUE ~ "non DE"))
            #Trending:
            if (is.null(trend_ligands_to_label) == FALSE) {
              
              bulk.data <- as.data.frame(t(seurat.obj@assays$RNA@meta.features[, grep("bulk", colnames(seurat.obj@assays$RNA@meta.features))]))
              
              if (ncol(bulk.data) == 0 | nrow(bulk.data) ==0) {
                stop( "Error: calc aggregate bulk data first before plotting with this function")
              }
              
              bulk.data$Celltype <- sapply(rownames(bulk.data), FUN = function(x) unlist(str_split(x, pattern = "__"))[2])
              bulk.data$Lesion <- sapply(rownames(bulk.data), FUN = function(x) unlist(str_split(x, pattern = "__"))[1])
              
              bulk.data <- pivot_longer(bulk.data, cols = colnames(bulk.data)[-which(colnames(bulk.data) %in% c("Celltype", "Lesion"))], names_to = "Gene", values_to = "CPM")
              
              trend.data <- bulk.data %>% filter(Gene %in% trend.genes)
              
              baseline.lesion <- as.character(levels(seurat.obj$Lesion)[1])       
              comparison.lesions <- as.character(unique(edges.df$Lesion))[-which(as.character(unique(edges.df$Lesion)) == condition_of_interest)]
              condition_of_interest = condition_of_interest
              
              #Only define trends as relevant to baseline and conditions in network:
              trend.data <- trend.data %>% filter(Lesion %in% c(baseline.lesion, condition_of_interest, comparison.lesions))
              
              max.gene.data <- trend.data %>% group_by(Gene, Lesion, Celltype) %>% summarize(mean(CPM)) %>% group_by(Gene) %>% summarize(max.lesion.average = max(`mean(CPM)`))
              mean.gene.data <- trend.data %>% group_by(Gene, Lesion, Celltype) %>% summarize(mean.lesion = mean(CPM))
              
              Trend.labeling.df <- as.data.frame(max.gene.data)
              Trend.labeling.df$Max.cell <- rep("", times = nrow(Trend.labeling.df))
              Trend.labeling.df$Trend <- rep("", times = nrow(Trend.labeling.df))
              for (g in 1:nrow(max.gene.data)) {
                
                gene <- as.character(max.gene.data[g,"Gene"])
                max.gene.average = as.numeric(max.gene.data[g,"max.lesion.average"])
                
                if (max.gene.average == 0) {
                  
                  Trend.labeling.df$Trend[g] <- "Not in data subset"
                  
                  next()
                }
                
                max.celltype <- mean.gene.data %>% filter(Gene == gene & mean.lesion == max.gene.average) %>% dplyr::pull(Celltype)
                max.lesion <- mean.gene.data %>% filter(Gene == gene & mean.lesion == max.gene.average) %>% dplyr::pull(Lesion)
                
                Trend.labeling.df$Max.cell[g] <- max.celltype
                
                #Define a few simple trends using cell type of max average expression:
                gene.trend.data <- mean.gene.data %>% filter(Gene == gene & Celltype == max.celltype)                  
                gene.trend.data <- gene.trend.data %>% pivot_wider(names_from = "Lesion", values_from = "mean.lesion")
                gene.trend.data <- gene.trend.data[,-which(colnames(gene.trend.data) %in% c("Gene", "Celltype"))]
                
                #Convert to percent of max:
                gene.trend.data <- as.data.frame(t(apply(gene.trend.data, MARGIN = 1, FUN = function(x) (x/max.gene.average)*100)))
                
                #Genes that go down relative to baseline:
                if (max.lesion == baseline.lesion) {
                  
                  if (gene.trend.data[,eval(condition_of_interest)] < 50 & any(gene.trend.data[,eval(comparison.lesions)] < 50)) {
                    Trend.labeling.df$Trend[g] <- "Common Downregulated"
                  }
                  
                  if (gene.trend.data[,eval(condition_of_interest)] >= 50 & any(gene.trend.data[,eval(comparison.lesions)] < 50)) {
                    Trend.labeling.df$Trend[g] <- "Comparison Downregulated"
                  }
                  
                  if (gene.trend.data[,eval(condition_of_interest)] < 50 & all(gene.trend.data[,eval(comparison.lesions)] >= 50)) {
                    Trend.labeling.df$Trend[g] <- "Unique Downregulated"
                  }
                  
                }
                
                #Genes that go up relative to baseline:
                if (max.lesion != baseline.lesion) {
                  
                  if (gene.trend.data[,eval(condition_of_interest)] >= 50 & all(gene.trend.data[,eval(comparison.lesions)] < 50)) {
                    Trend.labeling.df$Trend[g] <- "Unique Upregulated"
                  }
                  
                  if (gene.trend.data[,eval(condition_of_interest)] >= 50 & any(gene.trend.data[,eval(comparison.lesions)] >= 50)) {
                    Trend.labeling.df$Trend[g] <- "Common Upregulated"
                  }
                  
                  if (gene.trend.data[,eval(condition_of_interest)] < 50 & any(gene.trend.data[,eval(comparison.lesions)] >= 50)) {
                    Trend.labeling.df$Trend[g] <- "Comparison Upregulated"
                  }
                  
                  
                  
                }
                
              }
              
              
              Trend.labeling.df <- inner_join(Trend.labeling.df[,which(colnames(Trend.labeling.df) %in% c("Gene", "Max.cell", "Trend"))], edges.df, by = c("Gene" = "Ligand"))
              Trend.labeling.df <- Trend.labeling.df[,which(colnames(Trend.labeling.df) %in% c("Gene", "Trend", "from"))]
              Trend.labeling.df <- Trend.labeling.df[-which(duplicated(Trend.labeling.df$from)),]
              colnames(Trend.labeling.df) <- c("Ligand", "Label", "from")
              
              
              
              #Ligand Labeling:
              DE.labeling.df <- DE.labeling.df[,which(colnames(DE.labeling.df) %in% c("Ligand", "from", "DE_color_factor"))]
              colnames(DE.labeling.df) <- c("Ligand", "from", "Label")
              
              Ligand.labeling.df <- bind_rows(DE.labeling.df, Trend.labeling.df)
              Ligand.labeling.df$CellType <- sapply(Ligand.labeling.df$from, FUN = function(x) unlist(str_split(x, pattern = ":"))[1])
            }
            
            if (is.null(trend_ligands_to_label) == TRUE) {
              Ligand.labeling.df <- DE.labeling.df
              
            }
            
            #Ligand Labeling:
            sugi.node.labels.df$node_category <- ""
            for (r in 1:nrow(sugi.node.labels.df)) {
              
              cell <- sugi.node.labels.df$CellType[r]
              gene <- sugi.node.labels.df$node_gene[r]
              
              filtered.gene.result <- Ligand.labeling.df %>% filter(Ligand == gene & CellType == cell)
              
              if (nrow(filtered.gene.result) >0) {
                sugi.node.labels.df$node_category[r]<- filtered.gene.result$Label[1]
              }
            }
            
            #Remove Common Downregulated trend from labels to simplify interpretation:
            sugi.node.labels.df$node_category[which(sugi.node.labels.df$node_category %in% c("Common Downregulated"))] <- ""
            
            #Label nodes that are reasonably far away from each other, with preference for special categorized ligands:
            
            for (g in 1:length(unique(sugi.node.labels.df$node_gene))) {
              gene <- unique(sugi.node.labels.df$node_gene)[g]
              
              gene.nodes <- sugi.node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.sugi.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(desc(node_category),centroid_dist) %>% dplyr::pull(node_name)
                node.to.label <- group.nodes[1]
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% node.to.label),"node_label"] <- gene
                
              }
              
              sugi.node.labels.df[which(sugi.node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            unique(sugi.node.labels.df$node_category)
            sugi.node.labels.df$node_category <- factor(sugi.node.labels.df$node_category,
                                                        levels = c("",  "Unique DE", "Unique Upregulated", "Shared DE", "Common Upregulated", "Comparison Unique DE", "Comparison Upregulated"))
            
            g <- ggraph(merged.graph, layout = "manual",
                        x = sugi_layout[, "x"],
                        y = sugi_layout[, "y"]) + 
              geom_node_point(data = sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
              geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
              scale_shape_manual(values = node_shape_vals)+
              scale_color_manual(values = node.colors1) +
              geom_edge_diagonal(aes(color = foldchange_category), width = 0.5, alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              scale_edge_colour_manual(values = c("#0F1BDF", "#787878","#DA1B1B")) +
              expand_limits(y = c(min(sugi.node.labels.df$y), max(sugi.node.labels.df$y)+0.05*(range(sugi.node.labels.df$y)[2]-range(sugi.node.labels.df$y)[1]))) +
              ggnewscale::new_scale_color() +
              geom_text_repel(data = sugi.node.labels.df, aes(x = x, y = y, label = node_label, color = node_category), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5) +
              scale_color_manual(values = c("black", "red", "#FF7D7D",           
                                            "#781774", "#E378DF",
                                            "#2200FA", "#00D4FA")) +
              labs(title = paste(condition_of_interest, " vs ", str_c(as.character(unique(filtered.nodes.df$Lesion))[-which(as.character(unique(filtered.nodes.df$Lesion)) == condition_of_interest)],collapse = ", "), sep = ""))+
              theme(plot.title = element_text(hjust = 0.5))
            
            network.output.result <- list("Graph_object" = list(merged.graph),
                                          "Layout" = sugi_layout,
                                          "Node_label_df" = sugi.node.labels.df,
                                          "Plot" = g)
          }
          
          plot(g)
          return(network.output.result)
        }
        
        
        ##Edited to highlight both DE and Trending Ligands:
        plot_lig_rec_merged_network <- function(network, DE.table, seurat.obj, condition_of_interest,
                                                facet_by = "Lesion", node_shape = "Basic_Celltype", node_shape_vals = c(15, 8, 13, 16), node_size = 2, label_size = 2, node_colors = NA,
                                                create_DAG = TRUE, plot_layout = "sugiyama") {
          
          # Ligand-Receptor-Celltype Specific Node Plot:        
          message("Making initial ligand-receptor celltype-specific network. . . ")
          #Create node df:
          lesion.list <- unique(network[, facet_by])                                                 
          
          node_list <- c()
          
          for (i in 1:length(lesion.list)) {
            
            #Filter Edges:
            lesion.edges <- filter(network, eval(as.symbol(facet_by)) == lesion.list[i])                     
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet_by = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          #Merge Nodes:
          nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            nodes.df[,facet_by] <- factor(nodes.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
          }
          
          #Create Edges.df:
          edges.df <- network                                                                          
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            edges.df[,facet_by] <- factor(edges.df[,facet_by], 
                                          levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          edges.df <- edges.df[,which(colnames(edges.df) %in% c("from", "to", facet_by, "Ligand_Source", "Ligand", "Receptor_source", "Receptor","log10_Connection_product", "database_source"))]
          
          #Remove duplicate nichenet edges: (not sure why this happened. . . need to fix this eventually)
          edges.df <- edges.df %>% mutate(edgeID = paste(from, "__", to, "__", eval(as.symbol(facet_by)), sep = ""))
          
          duplicates <-  edges.df$edgeID[which(duplicated(edges.df$edgeID))]         
          if (length(duplicates) > 0 ) {
            duplicates <- edges.df[which(edges.df$edgeID %in% duplicates),]
            duplicate.max.connections <- duplicates %>% group_by(edgeID) %>% summarize(max.connection = max(log10_Connection_product))
            
            duplicates <- duplicates[-which(duplicated(duplicates$edgeID)),]
            duplicates <- full_join(duplicates, duplicate.max.connections)
            duplicates <- duplicates[,-which(colnames(duplicates) == "log10_Connection_product")]
            
            duplicates <- duplicates %>% rename(log10_Connection_product = "max.connection")
            
            edges.df <- edges.df[-which(edges.df$edgeID %in% duplicates$edgeID),]
            edges.df <- bind_rows(edges.df, duplicates)
          }
          
          
          # Scale edge width by percent of max connection strength:
          message("Scaling edges. . . ")
          
          max.expression.vals <- edges.df %>% group_by(Ligand, Receptor) %>% summarize(max.expression = max(log10_Connection_product))
          
          scaled_edges <- edges.df
          
          for (r in 1:nrow(max.expression.vals)) {
            
            connection.row <- max.expression.vals[r,]
            connection.max <- max.expression.vals[r,"max.expression"]
            
            connection.edges <- edges.df %>% filter(Ligand == connection.row$Ligand & Receptor == connection.row$Receptor)
            connection.edges$log10_Connection_product <- sapply(connection.edges$log10_Connection_product, FUN = function(x) as.numeric(round(10^(x)/(10^(connection.max)), digits = 3)))
            
            scaled_edges[which(scaled_edges$edgeID %in% connection.edges$edgeID),] <- connection.edges
          }
          
          edges.df <- scaled_edges
          colnames(edges.df)[which(colnames(edges.df) == "log10_Connection_product")] <- "Scaled_Connection_Product"
          
          
          #Add edge DE status:
          edges.df$DE <- FALSE
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            les.celltypes <- unique(lesion.edges$Ligand_Source)
            
            for (lc in 1:length(les.celltypes)) {
              
              les.cell <- les.celltypes[lc]
              
              les.cell.edges <- lesion.edges %>% filter(Ligand_Source == les.cell)
              
              les.cell.DE <- DE.table %>% filter(Comparison == les & CellType == les.cell & Color_factor == "Significantly Up" & logFC >= 1 & Reference != "Acetone") %>% dplyr::pull(Gene)
              
              DE.complex.genes <- deconvoluted %>% filter(gene_name %in% les.cell.DE) %>% dplyr::pull(complex_name) %>% unique()
              
              #Only label DE CellphoneDB
              les.cell.edges.DE <- les.cell.edges %>% filter(Ligand %in% unique(les.cell.DE, DE.complex.genes) & database_source == "CellphoneDB") %>% dplyr::pull(edgeID)
              
              if (length(les.cell.edges.DE) >0) {
                edges.df[which(edges.df$edgeID %in% les.cell.edges.DE), "DE"] <- TRUE
              }
              
            }
          }
          
          
          #Prune edges to only interconnected CellphoneDB & Nichenet edges:
          message("Pruning edges to interconnected cellphonedb and nichenet edges. . . ")
          
          filtered.edges.df <- as.data.frame(edges.df)[0,]
          
          for (l in 1:length(unique(edges.df[,facet_by]))) {
            
            les <- unique(edges.df[,facet_by])[l]
            
            lesion.edges <- edges.df %>% filter(eval(as.symbol(facet_by)) == les & Scaled_Connection_Product > 0)
            
            if (nrow(lesion.edges)>0) {
              
              #Prune nichenet edges:
              lesion.nichenet.edges <- lesion.edges %>% filter(database_source == "Nichenet" & Scaled_Connection_Product > 0) 
              
              if (nrow(lesion.nichenet.edges)>0){
                
                nichenet.edges.to.prune <- c()
                
                for (ne in 1:length(unique(lesion.nichenet.edges$from))) {
                  
                  nichenet.edge.from <- unique(lesion.nichenet.edges$from)[ne]
                  
                  upstream.edges <- lesion.edges %>% filter(to == nichenet.edge.from & Scaled_Connection_Product > 0 & database_source == "CellphoneDB")
                  
                  if (nrow(upstream.edges) == 0) {
                    nichenet.edges.to.prune <- c(nichenet.edges.to.prune, nichenet.edge.from)
                  }
                  
                }    
                
                lesion.nichenet.edges <- lesion.nichenet.edges %>% filter(!from %in% nichenet.edges.to.prune)
              }
              
              #Prune CellphoneDB edges:
              lesion.cellphonedb.edges <- lesion.edges %>% filter(database_source == "CellphoneDB" & Scaled_Connection_Product > 0) 
              
              if(nrow(lesion.cellphonedb.edges) >0) {
                
                cellphonedb.edges.to.prune <- c()
                
                for (ce in 1:length(unique(lesion.cellphonedb.edges$to))) {
                  
                  cellphonedb.edge.to <- unique(lesion.cellphonedb.edges$to)[ce]
                  
                  downstream.edges <- lesion.edges %>% filter(from == cellphonedb.edge.to & Scaled_Connection_Product > 0 & database_source == "Nichenet")
                  
                  if (nrow(downstream.edges) == 0) {
                    cellphonedb.edges.to.prune <- c(cellphonedb.edges.to.prune, cellphonedb.edge.to)
                  }
                  
                }
                
                lesion.cellphonedb.edges <- lesion.cellphonedb.edges %>% filter(!to %in% cellphonedb.edges.to.prune)
              }
              
              lesion.edges <- bind_rows(lesion.cellphonedb.edges, lesion.nichenet.edges)         
              
              filtered.edges.df <- bind_rows(filtered.edges.df, lesion.edges)
            }
          }
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.edges.df[,facet_by] <- factor(filtered.edges.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          #Duplicate node names columns so that each is preserved after graph formation:
          filtered.edges.df$from.named <- filtered.edges.df$from
          filtered.edges.df$to.named <- filtered.edges.df$to
          
          #Create filtered.nodes.df:
          node_list <- c()
          
          for (i in 1:length(unique(filtered.edges.df[,facet_by]))) {
            
            les <- as.character(unique(filtered.edges.df[,facet_by]))[i]
            #Filter Edges:
            lesion.edges <- filter(filtered.edges.df, eval(as.symbol(facet_by)) == les)                   
            
            #Create Nodes DFs:
            
            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
            
            nodes.df <- data.frame(node_name = nodes,
                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                   Facet = lesion.list[i])
            colnames(nodes.df)[3] <- facet_by
            
            #Assign Outputs:
            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
            
            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
            
          }
          
          
          #Merge Nodes:
          filtered.nodes.df <- purrr::reduce(mget(node_list), full_join)
          
          if (class(seurat.obj@meta.data[,facet_by]) == "factor") {
            filtered.nodes.df[,facet_by] <- factor(filtered.nodes.df[,facet_by], 
                                                   levels = levels(seurat.obj@meta.data[,facet_by]))
            
          }
          
          
          #Make Initial Graph:    
          RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          #Correct misconnected edges:
          message("Correcting tbl_graph misconnected edges. . . ")
          
          graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
          graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
          
          duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
          
          filtered.nodes.df <- rownames_to_column(filtered.nodes.df)
          
          for (r in 1:nrow(duplicates)) {
            
            les <- as.character(duplicates[r,facet_by])
            
            if (duplicates[r,"database_source"] == "CellphoneDB") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
            
            if (duplicates[r,"database_source"] == "Nichenet") {
              from <- duplicates[r,"from.named"]
              
              to <- duplicates[r,"to.named"]
              
              duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
              duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
              
            }
          }
          
          graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
          
          RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
          
          
          # Define Layout using merged data from all lesions:
          nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
          nodes <- nodes[-which(duplicated(nodes$node_name)),]  
          
          merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
          
          layout <- layout_with_stress(merged.graph)
          colnames(layout) <- c("x", "y")
          
          layout <- cbind(nodes, layout)
          
          ## Ligand Edge Label mapping:
          label.mapping.df <- filtered.edges.df
          label.mapping.df$label_x <- NA
          label.mapping.df$label_y <- NA
          label.mapping.df$Edge_label <- ""
          
          for (l in 1:nrow(label.mapping.df)) {
            
            lig.from <- label.mapping.df$from.named[l]
            lig.to <- label.mapping.df$to.named[l]
            
            node1_info <- layout %>% filter(node_name == lig.from)
            node2_info <- layout %>% filter(node_name == lig.to)
            
            label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
            label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
          }
          
          #Only label ligands that are reasonably far away from each other, with DE prioritization:
          max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
          
          for (l in 1:length(unique(label.mapping.df$Ligand))) {
            lig <- unique(label.mapping.df$Ligand)[l]
            
            lig.edges <- label.mapping.df %>% filter(Ligand == lig)
            
            lig.edges$Edge_label <- "temp"
            
            while(any(lig.edges$Edge_label =="temp")){
              
              #temporarily make a column for fraction of max distance from first ungrouped edge:
              grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
              
              grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.035) 
              
              #Label closest to centroid edges grouped by 5% of max layout distance:
              centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
              centroid_x <- sum(centroid_x)/length(centroid_x)
              centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
              centroid_y <- sum(centroid_y)/length(centroid_y)
              
              grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
              
              group.edges <- grouped.lig.edges %>% arrange(desc(DE), centroid_dist) %>% dplyr::pull(edgeID)
              Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
              
              lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
              
              lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
              
            }
            
            label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
            
          }
          
          #Copy coordinates to node df containing all faceted nodes:
          layout_full <- full_join(filtered.nodes.df, layout)
          
          if (create_DAG == TRUE)  {
            message("Pruning to directed acyclic graph. . . ")
            
            #Prune cycles and feedback loops to create a directed acyclic graph:
            
            ## Find Cycles in network:
            FindCycles = function(g) {
              Cycles = NULL
              for(v1 in V(g)) {
                if(degree(g, v1, mode="in") == 0) { next }
                GoodNeighbors = neighbors(g, v1, mode="out")
                GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
                for(v2 in GoodNeighbors) {
                  TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
                  TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
                  TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
                  Cycles  = c(Cycles, TempCyc)
                }
              }
              Cycles
            }
            merged.igraph <- as.igraph(merged.graph)
            network.cycles <- FindCycles(merged.igraph)
            
            if (length(network.cycles)>0) {
              
              
              nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
              
              for (c in 1:length(network.cycles)) {
                cycle <- network.cycles[[c]]
                
                named.cycle <- c()  
                
                for (n in 1:length(cycle)) {
                  node.num <- cycle[n]
                  node.info <- nodes[as.numeric(node.num),]
                  named.node <- node.info$node_name
                  
                  named.cycle <- c(named.cycle, named.node)
                }
                
                if (c == 1) {
                  named.network.cycles <-list(named.cycle)
                }
                if (c > 1) {
                  named.network.cycles <- append(named.network.cycles, list(named.cycle))
                }
              }
              
              #Edges causing cycles:
              
              for (c in 1:length(named.network.cycles)) {
                named.cycle <- named.network.cycles[[c]]
                
                first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
                last.edges <- filtered.edges.df %>% filter(from == named.cycle[length(named.cycle)-1] & to == named.cycle[length(named.cycle)])
                
                if (c ==1) {
                  cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                    first.edge.from = unique(first.edges$from.named),
                                                    first.edge.to = unique(first.edges$to.named),
                                                    first.edges.num= nrow(first.edges),
                                                    first.edges.database = unique(first.edges$database_source),
                                                    last.edge.from = unique(last.edges$from.named),
                                                    last.edge.to = unique(last.edges$to.named),
                                                    last.edges.num = nrow(last.edges),
                                                    last.edges.database = unique(last.edges$database_source))
                }
                
                if (c > 1) {
                  cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                 first.edge.from = unique(first.edges$from.named),
                                                 first.edge.to = unique(first.edges$to.named),
                                                 first.edges.num= nrow(first.edges),
                                                 first.edges.database = unique(first.edges$database_source),
                                                 last.edge.from = unique(last.edges$from.named),
                                                 last.edge.to = unique(last.edges$to.named),
                                                 last.edges.num = nrow(last.edges),
                                                 last.edges.database = unique(last.edges$database_source))
                  cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                }
                
              }
              cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(first.edge.from, "__", first.edge.to, sep = ""),
                                            last.edge = paste(last.edge.from, "__", last.edge.to, sep = ""))
              cyclic.edge.info.df <- mutate(cyclic.edge.info.df, cycle_id = paste(first.edge, ".....",last.edge, sep = ""))
              
              #Cycles seem to be due to weak nichenet connections at last edge:
              cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$last.edge)),]
              cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("last.edge.from", "last.edge.to"))]
              
              if (nrow(cyclic_edges_to_prune) >0 ) {
                for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                  
                  cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$last.edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$last.edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                  
                  if (length(cyclic.edgeID) >0) {
                    filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                  }
                  rm(cyclic.edgeID)
                }
              }
              
              #Re-make Graph without above pruned cycle edges:    
              RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
              
              #Correct misconnected edges:
              graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
              graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
              
              duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
              
              if (nrow(duplicates)>0) {
                for (r in 1:nrow(duplicates)) {
                  
                  les <- as.character(duplicates[r,facet_by])
                  
                  if (duplicates[r,"database_source"] == "CellphoneDB") {
                    from <- duplicates[r,"from.named"]
                    
                    to <- duplicates[r,"to.named"]
                    
                    duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                    duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                    
                  }
                  
                  if (duplicates[r,"database_source"] == "Nichenet") {
                    from <- duplicates[r,"from.named"]
                    
                    to <- duplicates[r,"to.named"]
                    
                    duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                    duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                    
                  }
                }
                
                graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
              }
              
              RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
              
              # Define Layout using merged data from all lesions:
              nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
              nodes <- nodes[-which(duplicated(nodes$node_name)),]  
              
              merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
              
              layout <- layout_with_stress(merged.graph)
              colnames(layout) <- c("x", "y")
              
              layout <- cbind(nodes, layout)
              
              #Copy coordinates to node df containing all faceted nodes:
              layout_full <- full_join(filtered.nodes.df, layout)
            } 
            
            ## Above function doesn't find all cycles though! Use igraph to find feedback edges requiring pruning for DAG creation:
            merged.igraph <- as.igraph(merged.graph)
            cycles <- feedback_arc_set(merged.igraph)
            cycle.edge.ids <- as_ids(cycles)
            
            if (length(cycle.edge.ids)>0){
              
              #convert to list:
              for (e in 1:length(cycle.edge.ids)) {
                if (e ==1){
                  cyclic.edges <- list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]])))
                  
                }
                
                if (e>1){
                  cyclic.edges <-  append(cyclic.edges, list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]]))))
                }
              }   
              
              #Add names:
              for (c in 1:length(cyclic.edges)) {
                cycle <- cyclic.edges[[c]]
                
                named.cycle <- c()  
                
                for (n in 1:length(cycle)) {
                  node.num <- cycle[n]
                  node.info <- nodes[as.numeric(node.num),]
                  named.node <- node.info$node_name
                  
                  named.cycle <- c(named.cycle, named.node)
                }
                
                if (c == 1) {
                  named.network.cycles <-list(named.cycle)
                }
                if (c > 1) {
                  named.network.cycles <- append(named.network.cycles, list(named.cycle))
                }
              }
              
              for (c in 1:length(named.network.cycles)) {
                named.cycle <- named.network.cycles[[c]]
                
                first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
                
                if (c ==1) {
                  cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                    edge.from = unique(first.edges$from.named),
                                                    edge.to = unique(first.edges$to.named),
                                                    edges.num= nrow(first.edges),
                                                    edges.database = unique(first.edges$database_source))
                  
                }
                
                if (c > 1) {
                  cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                 edge.from = unique(first.edges$from.named),
                                                 edge.to = unique(first.edges$to.named),
                                                 edges.num= nrow(first.edges),
                                                 edges.database = unique(first.edges$database_source))
                  
                  cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                }
                
              }
              
              cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(edge.from, "__", edge.to, sep = ""))
              
              #All cycles look to be due to weak nichenet connections at last edge:
              cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$first.edge)),]
              cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("edge.from", "edge.to"))]
              
              for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                
                cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                
                if (length(cyclic.edgeID) >0) {
                  filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                }
                rm(cyclic.edgeID)
              }
            }
            # Make Final Graph:    
            message("Making final graph. . . ")
            
            #add shapes to nodes if you want:
            if (!is.null(node_shape)) {
              
              node.shape.df <- DE.table[,which(colnames(DE.table) %in% c("CellType", node_shape))]  
              node.shape.df <- node.shape.df[-which(duplicated(node.shape.df$CellType)),]
              
              filtered.nodes.df <- left_join(filtered.nodes.df, node.shape.df)
              
              if (class(seurat.obj@meta.data[,node_shape]) == "factor") {
                filtered.nodes.df[node_shape] <- factor(filtered.nodes.df[,node_shape], 
                                                        levels = levels(seurat.obj@meta.data[,node_shape]))
              }
            }
            
            RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            #Correct misconnected edges:
            graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
            graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
            
            duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
            
            if (nrow(duplicates)>0) {
              
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,facet_by])
                
                if (duplicates[r,"database_source"] == "CellphoneDB") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
                
                if (duplicates[r,"database_source"] == "Nichenet") {
                  from <- duplicates[r,"from.named"]
                  
                  to <- duplicates[r,"to.named"]
                  
                  duplicates[r,"from"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                  duplicates[r,"to"] <- filtered.nodes.df %>% filter(eval(as.symbol(facet_by)) == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                  
                }
              }
              
              graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
            }
            
            RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
            
            
            # Define Layout using merged data from all lesions:
            nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType", node_shape, "node_label"))]
            nodes <- nodes[-which(duplicated(nodes$node_name)),]  
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
            layout <- layout_with_stress(merged.graph)
            colnames(layout) <- c("x", "y")
            
            layout <- cbind(nodes, layout)
            
            ## Ligand Edge Label mapping:
            label.mapping.df <- filtered.edges.df
            label.mapping.df$label_x <- NA
            label.mapping.df$label_y <- NA
            label.mapping.df$Edge_label <- ""
            
            for (l in 1:nrow(label.mapping.df)) {
              
              lig.from <- label.mapping.df$from.named[l]
              lig.to <- label.mapping.df$to.named[l]
              
              node1_info <- layout %>% filter(node_name == lig.from)
              node2_info <- layout %>% filter(node_name == lig.to)
              
              label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
              label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
            }
            
            #Only label ligands that are reasonably far away from each other, with DE prioritization:
            max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
            
            for (l in 1:length(unique(label.mapping.df$Ligand))) {
              lig <- unique(label.mapping.df$Ligand)[l]
              
              lig.edges <- label.mapping.df %>% filter(Ligand == lig)
              
              lig.edges$Edge_label <- "temp"
              
              while(any(lig.edges$Edge_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.edges <- grouped.lig.edges %>% arrange(desc(DE),centroid_dist) %>% dplyr::pull(edgeID)
                Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
                
                lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                
                lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                
              }
              
              label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
              
            }
            
            ## Node Labels:
            node.labels.df <- layout
            node.labels.df$node_label <- ""
            
            node.labels.df$node_gene <- sapply(node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            #Only label nodes that are reasonably far away from each other:
            
            for (g in 1:length(unique(node.labels.df$node_gene))) {
              gene <- unique(node.labels.df$node_gene)[g]
              
              gene.nodes <- node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.035) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                
              }
              
              node.labels.df[which(node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            #Copy coordinates to node df containing all faceted nodes:
            layout_full <- full_join(filtered.nodes.df, node.labels.df)
            
            
            ### Fold change comparisons for merged graph comparisons:
            filtered.edges.df$from_to_joint <- paste(filtered.edges.df$from.named, filtered.edges.df$to.named, sep = "__")
            
            filtered.edges.df <- filtered.edges.df[,-which(colnames(filtered.edges.df) %in% c("edgeID", "DE"))]
            filtered.edges.df <- filtered.edges.df %>% pivot_wider(values_from = "Scaled_Connection_Product", names_from = "Lesion")
            filtered.edges.df[is.na(filtered.edges.df)] <- 0
            
            lesions <- as.character(unique(filtered.nodes.df$Lesion))
            
            filtered.edges.df <- filtered.edges.df %>% mutate(fold_change = (eval(as.symbol(lesions[which(lesions == condition_of_interest)]))+0.05)/(eval(as.symbol(lesions[-which(lesions == condition_of_interest)]))+0.05))
            
            filtered.edges.df <- filtered.edges.df %>% mutate(foldchange_category = case_when(log2(fold_change) < -2 ~ "Down",
                                                                                              log2(fold_change) >= -2 & log2(fold_change) <= 2  ~ "Similar",
                                                                                              log2(fold_change) > 2 ~ "Up"))
            
            
            
            merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
            
          }
          
          
          #Set point colors for plotting:
          GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
          
          ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
          
          if (all(is.na(node_colors)) == TRUE) {
            node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
          }
          if (all(is.na(node_colors)) == FALSE) {
            node.colors1 <- node_colors
          }
          
          if (plot_layout == "stress") {
            
            ##Plot With Facets:
            g <- ggraph(RL_graph, layout = "manual",
                        x = layout_full[, "x"],
                        y = layout_full[, "y"]) + 
              geom_node_point(data = layout, color = "gray", alpha = 0.2, size = node_size) + 
              geom_node_point(aes(color = CellType), size = node_size) + 
              geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
              scale_edge_width(range = c(0.1, 0.5)) +
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              facet_nodes(~Lesion) +
              scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
              new_scale_color()+
              scale_color_manual(values = c("#000000", "#750000")) +
              geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = label_size, alpha =1, box.padding = 0)
            
            network.output.result <- list("Graph_object" = list(RL_graph),
                                          "Layout" = layout_full,
                                          "Plot" = g)
            
          }
          
          if (plot_layout == "sugiyama") {
            #Sugiyama Layout:
            sugi_layout <- layout_with_sugiyama(merged.graph, hgap = 500)
            sugi_layout <- as.data.frame(sugi_layout$layout)
            colnames(sugi_layout) <- c("x", "y")
            
            sugi_layout <- cbind(nodes, sugi_layout)
            
            ## Sugi Node Labels:
            sugi.node.labels.df <- sugi_layout
            sugi.node.labels.df$node_label <- ""
            
            sugi.node.labels.df$node_gene <- sapply(sugi.node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
            
            max.sugi.layout.dist <- ((max(sugi_layout$x)-min(sugi_layout$x))^2+(max(sugi_layout$y)-min(sugi_layout$y))^2)^(0.5)
            
            
            #Color DE and Trending Genes in each condition:
            
            
            #Compare conditions to label DE and trending Ligands:
            #DE:
            DE.labeling.df <- edges.df %>% filter(DE == TRUE & database_source == "CellphoneDB")
            DE.labeling.df <- DE.labeling.df[,which(colnames(DE.labeling.df) %in% c("Ligand", "from", "Lesion", "edgeID"))]
            
            if (nrow(DE.labeling.df)>0) {
              
              DE.labeling.df <- DE.labeling.df %>% pivot_wider(names_from = "Lesion", values_from = "edgeID") %>% as.data.frame() %>% suppressWarnings()
              
              labeling.data <- DE.labeling.df[,-which(colnames(DE.labeling.df) %in% c("Ligand", "from"))] 
              
              for (col in 1:ncol(labeling.data)) {
                
                labeling.data[,col] <- sapply(labeling.data[,col], FUN = function (x) {
                  
                  x <- str_c(unlist(x), collapse = "___")
                  
                  
                  return(x)
                })
              }
              for (col in 1:ncol(labeling.data)) {
                
                labeling.data[,col] <- sapply(labeling.data[,col], FUN = function (x) {
                  if (x != "") {
                    x <- "DE"
                  }
                  
                  if (x == "") {
                    x <- "Non-DE"
                  }
                  
                  
                  return(x)
                })
              }
              DE.labeling.df[,-which(colnames(DE.labeling.df) %in% c("Ligand", "from"))]  <- labeling.data
              
              
              DE.labeling.df <- DE.labeling.df %>% mutate(DE_color_factor = case_when(eval(as.symbol(condition_of_interest)) == "DE" & eval(as.symbol(colnames(DE.labeling.df)[-which(colnames(DE.labeling.df) %in% c(condition_of_interest, "Ligand", "from"))])) == "Non-DE" ~ "Unique DE",
                                                                                      eval(as.symbol(condition_of_interest)) == "DE" & any(eval(as.symbol(colnames(DE.labeling.df)[-which(colnames(DE.labeling.df) %in% c(condition_of_interest, "Ligand", "from"))])) == "DE") ~ "Shared DE",
                                                                                      eval(as.symbol(condition_of_interest)) == "Non-DE" & any(eval(as.symbol(colnames(DE.labeling.df)[-which(colnames(DE.labeling.df) %in% c(condition_of_interest, "Ligand", "from"))])) == "DE") ~ "Comparison Unique DE",
                                                                                      TRUE ~ "non DE"))
            }
            if (nrow(DE.labeling.df)==0) {
              colnames(DE.labeling.df)[4] <- "DE_color_factor"
            } 
            #Trending:
            bulk.data <- as.data.frame(t(seurat.obj@assays$RNA@meta.features[, grep("bulk", colnames(seurat.obj@assays$RNA@meta.features))]))
            
            if (ncol(bulk.data) == 0 | nrow(bulk.data) ==0) {
              stop( "Error: calc aggregate bulk data first before plotting with this function")
            }
            
            bulk.data$Celltype <- sapply(rownames(bulk.data), FUN = function(x) unlist(str_split(x, pattern = "__"))[2])
            bulk.data$Lesion <- sapply(rownames(bulk.data), FUN = function(x) unlist(str_split(x, pattern = "__"))[1])
            
            bulk.data <- pivot_longer(bulk.data, cols = colnames(bulk.data)[-which(colnames(bulk.data) %in% c("Celltype", "Lesion"))], names_to = "Gene", values_to = "CPM")
            
            trend.data <- bulk.data %>% filter(Gene %in% unique(sugi.node.labels.df$node_gene))
            
            
            baseline.lesion <- as.character(levels(seurat.obj$Lesion)[1])       
            comparison.lesions <- as.character(unique(edges.df$Lesion))[-which(as.character(unique(edges.df$Lesion)) == condition_of_interest)]
            condition_of_interest = condition_of_interest
            
            #Only define trends as relevant to baseline and conditions in network:
            trend.data <- trend.data %>% filter(Lesion %in% c(baseline.lesion, condition_of_interest, comparison.lesions))
            
            max.gene.data <- trend.data %>% group_by(Gene, Lesion, Celltype) %>% summarize(mean(CPM)) %>% group_by(Gene) %>% summarize(max.lesion.average = max(`mean(CPM)`))
            mean.gene.data <- trend.data %>% group_by(Gene, Lesion, Celltype) %>% summarize(mean.lesion = mean(CPM))
            
            Trend.labeling.df <- as.data.frame(max.gene.data)
            Trend.labeling.df$Max.cell <- rep("", times = nrow(Trend.labeling.df))
            Trend.labeling.df$Trend <- rep("", times = nrow(Trend.labeling.df))
            
            for (g in 1:nrow(max.gene.data)) {
              
              gene <- as.character(max.gene.data[g,"Gene"])
              max.gene.average = as.numeric(max.gene.data[g,"max.lesion.average"])
              
              if (max.gene.average == 0) {
                
                Trend.labeling.df$Trend[g] <- "Not in data subset"
                
                next()
              }
              
              max.celltype <- mean.gene.data %>% filter(Gene == gene & mean.lesion == max.gene.average) %>% dplyr::pull(Celltype)
              max.lesion <- mean.gene.data %>% filter(Gene == gene & mean.lesion == max.gene.average) %>% dplyr::pull(Lesion)
              
              Trend.labeling.df$Max.cell[g] <- max.celltype
              
              #Define a few simple trends using cell type of max average expression:
              gene.trend.data <- mean.gene.data %>% filter(Gene == gene & Celltype == max.celltype)                  
              gene.trend.data <- gene.trend.data %>% pivot_wider(names_from = "Lesion", values_from = "mean.lesion")
              gene.trend.data <- gene.trend.data[,-which(colnames(gene.trend.data) %in% c("Gene", "Celltype"))]
              
              #Convert to percent of max:
              gene.trend.data <- as.data.frame(t(apply(gene.trend.data, MARGIN = 1, FUN = function(x) (x/max.gene.average)*100)))
              
              #Genes that go down relative to baseline:
              if (max.lesion == baseline.lesion) {
                
                if (gene.trend.data[,eval(condition_of_interest)] < 25 & any(gene.trend.data[,eval(comparison.lesions)] < 25)) {
                  Trend.labeling.df$Trend[g] <- "Common Downregulated"
                }
                
                if (gene.trend.data[,eval(condition_of_interest)] >= 50 & any(gene.trend.data[,eval(comparison.lesions)] < 25)) {
                  Trend.labeling.df$Trend[g] <- "Comparison Downregulated"
                }
                
                if (gene.trend.data[,eval(condition_of_interest)] < 25 & all(gene.trend.data[,eval(comparison.lesions)] >= 50)) {
                  Trend.labeling.df$Trend[g] <- "Unique Downregulated"
                }
                if (gene.trend.data[,eval(condition_of_interest)] >= 25 & all(gene.trend.data[,eval(comparison.lesions)] >= 25)) {
                  Trend.labeling.df$Trend[g] <- "No Significant Trend"
                }
                
              }
              
              #Genes that go up relative to baseline:
              if (max.lesion != baseline.lesion) {
                
                if (gene.trend.data[,eval(condition_of_interest)] >= 75 & all(gene.trend.data[,eval(comparison.lesions)] < 50) & all(gene.trend.data[,eval(baseline.lesion)] < 50)) {
                  Trend.labeling.df$Trend[g] <- "Unique Upregulated"
                }
                
                if (gene.trend.data[,eval(condition_of_interest)] >= 75 & any(gene.trend.data[,eval(comparison.lesions)] >= 75) & all(gene.trend.data[,eval(baseline.lesion)] < 50)) {
                  Trend.labeling.df$Trend[g] <- "Common Upregulated"
                }
                
                if (gene.trend.data[,eval(condition_of_interest)] < 50 & any(gene.trend.data[,eval(comparison.lesions)] >= 75) & all(gene.trend.data[,eval(baseline.lesion)] < 50)) {
                  Trend.labeling.df$Trend[g] <- "Comparison Upregulated"
                }
                if (gene.trend.data[,eval(condition_of_interest)] <= 75 & all(gene.trend.data[,eval(comparison.lesions)] <= 75) | all(gene.trend.data[,eval(baseline.lesion)] > 50)) {
                  Trend.labeling.df$Trend[g] <- "No Significant Trend"
                }
                
                
              }
              
            }
            
            
            Trend.labeling.df <- inner_join(Trend.labeling.df[,which(colnames(Trend.labeling.df) %in% c("Gene", "Max.cell", "Trend"))], edges.df, by = c("Gene" = "Ligand"))
            Trend.labeling.df <- Trend.labeling.df %>% filter(Max.cell == Ligand_Source)
            
            Trend.labeling.df <- Trend.labeling.df[,which(colnames(Trend.labeling.df) %in% c("Gene", "Trend", "from"))]
            Trend.labeling.df <- Trend.labeling.df[-which(duplicated(Trend.labeling.df$from)),]
            colnames(Trend.labeling.df) <- c("Ligand", "Label", "from")
            
            
            #Ligand Labeling:
            DE.labeling.df <- DE.labeling.df[,which(colnames(DE.labeling.df) %in% c("Ligand", "from", "DE_color_factor"))]
            colnames(DE.labeling.df) <- c("Ligand", "from", "Label")
            
            Ligand.labeling.df <- bind_rows(DE.labeling.df, Trend.labeling.df)
            Ligand.labeling.df$CellType <- sapply(Ligand.labeling.df$from, FUN = function(x) unlist(str_split(x, pattern = ":"))[1])
            
            
            
            
            #Ligand Labeling:
            sugi.node.labels.df$node_category <- ""
            for (r in 1:nrow(sugi.node.labels.df)) {
              
              cell <- sugi.node.labels.df$CellType[r]
              gene <- sugi.node.labels.df$node_gene[r]
              
              filtered.gene.result <- Ligand.labeling.df %>% filter(Ligand == gene & CellType == cell)
              
              if (nrow(filtered.gene.result) >0) {
                sugi.node.labels.df$node_category[r]<- filtered.gene.result$Label[1]
              }
            }
            
            #Remove Common Downregulated trend from labels to simplify interpretation:
            sugi.node.labels.df$node_category[which(sugi.node.labels.df$node_category %in% c("Common Downregulated", "No Significant Trend"))] <- ""
            
            #Label nodes that are reasonably far away from each other, with preference for special categorized ligands:
            for (g in 1:length(unique(sugi.node.labels.df$node_gene))) {
              gene <- unique(sugi.node.labels.df$node_gene)[g]
              
              gene.nodes <- sugi.node.labels.df %>% filter(node_gene == gene)
              
              gene.nodes$node_label <- "temp"
              
              while(any(gene.nodes$node_label =="temp")){
                
                #temporarily make a column for fraction of max distance from first ungrouped edge:
                grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.sugi.layout.dist), digits = 3))
                
                grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.05) 
                
                #Label closest to centroid edges grouped by 5% of max layout distance:
                centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                centroid_x <- sum(centroid_x)/length(centroid_x)
                centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                centroid_y <- sum(centroid_y)/length(centroid_y)
                
                grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                
                group.nodes <- grouped.gene.nodes %>% arrange(desc(node_category),centroid_dist) %>% dplyr::pull(node_name)
                node.to.label <- group.nodes[1]
                
                gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                
                gene.nodes[which(gene.nodes$node_name %in% node.to.label),"node_label"] <- gene
                
              }
              
              sugi.node.labels.df[which(sugi.node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
              
            }
            
            #Make sure all nodes with special DE/Trend status are labeled:
            sugi.node.labels.df[which(sugi.node.labels.df$node_category != ""),] <- sugi.node.labels.df %>% filter(node_category != "") %>% mutate(node_label = node_gene)
            
            sugi.node.labels.df$node_category <- factor(sugi.node.labels.df$node_category,
                                                        levels = c("",  "Unique DE", "Unique Upregulated", "Shared DE", "Common Upregulated", "Comparison Unique DE", "Comparison Upregulated"))
            
            g <- ggraph(merged.graph, layout = "manual",
                        x = sugi_layout[, "x"],
                        y = sugi_layout[, "y"]) + 
              geom_node_point(data = sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = node_size) + 
              geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = node_size) + 
              scale_shape_manual(values = node_shape_vals)+
              scale_color_manual(values = node.colors1) +
              geom_edge_diagonal(aes(color = foldchange_category), width = 0.5, alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
              theme_graph(base_family = "AvantGarde") + 
              theme(legend.position = "bottom") +
              scale_edge_colour_manual(values = c("#0F1BDF", "#787878","#DA1B1B")) +
              expand_limits(y = c(min(sugi.node.labels.df$y), max(sugi.node.labels.df$y)+0.05*(range(sugi.node.labels.df$y)[2]-range(sugi.node.labels.df$y)[1]))) +
              ggnewscale::new_scale_color() +
              scale_color_manual(values = c("black", "red", "#D68800",           
                                            "#781774", "#E378DF",
                                            "#2200FA", "#00D4FA")) +
              geom_text_repel(data = sugi.node.labels.df, aes(x = x, y = y, label = node_label, color = node_category), nudge_y = 0.15, size = label_size, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5) +
              labs(title = paste(condition_of_interest, " vs ", str_c(as.character(unique(filtered.nodes.df$Lesion))[-which(as.character(unique(filtered.nodes.df$Lesion)) == condition_of_interest)],collapse = ", "), sep = ""))+
              theme(plot.title = element_text(hjust = 0.5))
            
            network.output.result <- list("Graph_object" = list(merged.graph),
                                          "Layout" = sugi_layout,
                                          "Node_label_df" = sugi.node.labels.df,
                                          "Plot" = g)
          }
          
          plot(g)
          return(network.output.result)
        }
        
        
        
        
        set.seed(100)
        network_4lesion.result <- plot_lig_rec_network(network, DE.table, DE.input)

        ##Repeat just for Irritant and Day 2 Allergy lesions:
            day2_network <- network %>% filter(Lesion %in% c("Irritant","Day2_Allergy"))
            day2_DE.table <- DE.table %>% filter(Comparison %in% c("Irritant", "Day2_Allergy"))
                  
            set.seed(100)
            network_2lesion.result <- plot_lig_rec_network(day2_network, day2_DE.table, DE.input)
            
            set.seed(100)
            network_2lesion.merged.result <- plot_lig_rec_merged_network(day2_network, day2_DE.table, DE.input, condition_of_interest = "Day2_Allergy", 
                                                                         node_colors = c("#FF820D", "#C80707", "#FF4646", "#267605", "#876241", "#436CEB", "#0CBB06", "#8AAF58", "#9D0BE1"))
            
            network_2lesion.merged.result$Plot

            #With every node labeled:
            node_shape_vals = c(15, 8, 13, 16)
            GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
            
            ColorCount = as.numeric(length(unique(network_2lesion.merged.result$Layout$CellType)))
            node.colors1 <- c("#FF820D", "#C80707", "#FF4646", "#267605", "#876241", "#436CEB", "#0CBB06", "#8AAF58", "#9D0BE1")
            
            set.seed(100)
            network_2lesion.merged.all.labeled <- ggraph(network_2lesion.merged.result$Graph_object[[1]], layout = "manual",
                                                         x = network_2lesion.merged.result$Layout[, "x"],
                                                         y = network_2lesion.merged.result$Layout[, "y"]) + 
                                                    geom_node_point(data = network_2lesion.merged.result$Layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
                                                    geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
                                                    scale_shape_manual(values = node_shape_vals)+
                                                    scale_color_manual(values = node.colors1) +
                                                    geom_edge_diagonal(aes(color = foldchange_category), width = 0.5, alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
                                                    theme_graph(base_family = "AvantGarde") + 
                                                    theme(legend.position = "bottom") +
                                                    scale_edge_colour_manual(values = c("#0F1BDF", "#787878","#DA1B1B")) +
                                                    expand_limits(y = c(min(network_2lesion.merged.result$Node_label_df$y), max(network_2lesion.merged.result$Node_label_df$y)+0.05*(range(network_2lesion.merged.result$Node_label_df$y)[2]-range(network_2lesion.merged.result$Node_label_df$y)[1]))) +
                                                    ggnewscale::new_scale_color() +
                                                    geom_text_repel(data = network_2lesion.merged.result$Node_label_df, aes(x = x, y = y, label = node_gene, color = node_category), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5) +
                                                    scale_color_manual(values =  c("black", "red", "#D68800",           
                                                                                   "#781774", "#E378DF",
                                                                                   "#2200FA", "#00D4FA")) +
                                                    labs(title = "Day 2 Allergy vs Irritant" ) +
                                                    theme(plot.title = element_text(hjust = 0.5))
            
            #With only DE Trend nodes labeled:
            node_shape_vals = c(15, 8, 13, 16)
            GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
            
            ColorCount = as.numeric(length(unique(network_2lesion.merged.result$Layout$CellType)))
            node.colors1 <- c("#FF820D", "#C80707", "#FF4646", "#267605", "#876241", "#436CEB", "#0CBB06", "#8AAF58", "#9D0BE1")
            
            DE.trending.nodes <- network_2lesion.merged.result$Node_label_df %>% filter(node_category != "")
            
            set.seed(100)
            network_2lesion.merged.DE.trend.labeled <- ggraph(network_2lesion.merged.result$Graph_object[[1]], layout = "manual",
                                                         x = network_2lesion.merged.result$Layout[, "x"],
                                                         y = network_2lesion.merged.result$Layout[, "y"]) + 
                                                    geom_node_point(data = network_2lesion.merged.result$Layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
                                                    geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
                                                    scale_shape_manual(values = node_shape_vals)+
                                                    scale_color_manual(values = node.colors1) +
                                                    geom_edge_diagonal(aes(color = foldchange_category), width = 0.5, alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
                                                    theme_graph(base_family = "AvantGarde") + 
                                                    theme(legend.position = "bottom") +
                                                    scale_edge_colour_manual(values = c("#0F1BDF", "#787878","#DA1B1B")) +
                                                    expand_limits(y = c(min(network_2lesion.merged.result$Node_label_df$y), max(network_2lesion.merged.result$Node_label_df$y)+0.05*(range(network_2lesion.merged.result$Node_label_df$y)[2]-range(network_2lesion.merged.result$Node_label_df$y)[1]))) +
                                                    ggnewscale::new_scale_color() +
                                                    geom_text_repel(data = DE.trending.nodes, aes(x = x, y = y, label = node_gene, color = node_category), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5) +
                                                    scale_color_manual(values =  c("red", "#D68800",           
                                                                                   "#781774", "#E378DF",
                                                                                   "#2200FA", "#00D4FA")) +
                                                    labs(title = "Day 2 Allergy vs Irritant" ) +
                                                    theme(plot.title = element_text(hjust = 0.5))
            
            network_2lesion.merged.DE.trend.labeled
            
            Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
            
            plot_violin_bar_srt(DE.input, gene = "CCL16", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype")
            plot_violin_srt(DE.input, gene = "ANXA1", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype")
            plot_barplot_srt(DE.input, gene = "CCL13", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
            

            
        setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/rl_network")
        #ggsave(network_4lesion.result$Plot, width = 80, height = 40,units = "cm",filename = "network_4lesion.pdf", limitsize = FALSE)  
        #ggsave(network_2lesion.result$Plot, width = 80, height = 40,units = "cm",filename = "network_2lesion.new.png", limitsize = FALSE)  
        #ggsave(network_2lesion.merged.result$Plot, width = 80, height = 40,units = "cm",filename = "network_2lesion_merged.new.png", limitsize = FALSE)  
        #ggsave(network_2lesion.merged.result$Plot, width = 45, height = 35,units = "cm",filename = "network_2lesion_merged.new.pdf", limitsize = FALSE)  
        #ggsave(network_2lesion.merged.DE.trend.labeled, width = 45, height = 35,units = "cm",filename = "network_2lesion_merged.de.trend.labled.new.pdf", limitsize = FALSE)  
        #ggsave(network_2lesion.merged.all.labeled, width = 200, height = 200,units = "cm",filename = "network_2lesion_merged.large.new.pdf", limitsize = FALSE)  
        
        Day2.edges <- network_2lesion.merged.result$Graph_object[[1]] %>% activate(edges) %>% as.data.frame()
        
        
        ####### Validate by showing IL4 network signals:
        contact_derm_DE <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE.txt")
        
        IL4.edges <- network_2lesion.merged.result$Graph_object[[1]] %>% activate(edges) %>% as.data.frame() %>% filter(Ligand %in% c("IL4", "IL13"))
        IL4.induced.edges <- network_2lesion.merged.result$Graph_object[[1]] %>% activate(edges) %>% as.data.frame() %>% filter(Receptor %in% c("IL13_receptor", "IL4_receptor") & database_source == "Nichenet")
        
        Day2_IL4.IL13_network <- network %>% filter(Lesion %in% c("Irritant","Day2_Allergy")) %>% filter(Receptor %in% c("IL13_receptor", "IL4_receptor") | Ligand %in% c("IL4", "IL13"))
        
        set.seed(100)
        network_IL4.IL13.merged.result <- plot_lig_rec_merged_network(Day2_IL4.IL13_network, day2_DE.table, DE.input, condition_of_interest = "Day2_Allergy", label_size = 5, node_size = 5, node_shape_vals = c(15, 13, 16))
        
        network_IL4.IL13.merged.result$Plot
        
        #ggsave(network_IL4.IL13.merged.result$Plot, width = 20, height = 15,units = "cm",filename = "network_IL4.IL13_merged.new.pdf", limitsize = FALSE)  
          
        load("Data/RL_network2/nichenet_full_network_frac0.001.Rdata.Rdata")
        
        #Examine Myeloid Cells by Regular nichenet:
              seurat_obj <- DE.input
              DE.table.input <- contact_derm_DE
              gene.frac.cutoff = 0.001
              comp = "Day2_Allergy"
              ref = "Nonlesional"
              
              Idents(seurat_obj) <- "Lesion"
              seurat_subset <- subset(seurat_obj, idents = c(comp, ref))
              celltypes.to.analyze <- unique(seurat_subset@meta.data$Semi_Detailed_Celltype)
              Idents(seurat_subset) <- "Semi_Detailed_Celltype"
                
              ## receiver = individual celltype:
                  receiver = "Myeloid"
                   
                  print(paste(". . . . . . . . Analyzing", receiver, sep = " "))
                  
                  expressed_genes_receiver = get_expressed_genes(receiver, seurat_subset, pct = gene.frac.cutoff, assay_oi = "RNA")
                  
                  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
                  
                  ## senders = all celltypes in data:
                  sender_celltypes = as.character(unique(Idents(seurat_subset)))
                  list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_subset, gene.frac.cutoff, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
                  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
                  
                  #Find DE genes within Receiver Celltype:
                  seurat_obj_receiver= subset(seurat_subset, idents = receiver)
                  seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Lesion"]])
                  
                  subset.condition_oi = comp
                  subset.condition_ref = ref
                  
                  
                  if (all(subset.condition_oi %in% unique(seurat_obj_receiver$Lesion)) == FALSE & any(subset.condition_oi %in% unique(seurat_obj_receiver$Lesion)) == TRUE) {
                    subset.condition_oi = subset.condition_oi[which(subset.condition_oi %in% unique(seurat_obj_receiver$Lesion))]
                  }
                  if (all(subset.condition_ref %in% unique(seurat_obj_receiver$Lesion)) == FALSE & any(subset.condition_ref %in% unique(seurat_obj_receiver$Lesion)) == TRUE) {
                    subset.condition_ref = subset.condition_ref[which(subset.condition_ref %in% unique(seurat_obj_receiver$Lesion))]
                  }
                  
                  cells.oi.num <- seurat_obj_receiver@meta.data %>% filter(seurat_obj_receiver$Lesion %in% subset.condition_oi) %>% nrow() %>% as.numeric()
                  cells.ref.num <- seurat_obj_receiver@meta.data %>% filter(seurat_obj_receiver$Lesion %in% subset.condition_ref) %>% nrow() %>% as.numeric()
                  
                  if (cells.oi.num < 3 | cells.ref.num < 3) {
                    active_ligand_target_links_df <- data.frame(ligand = NA,
                                                                target = NA, 
                                                                weight = NA,
                                                                Receiver = receiver,
                                                                Condition = comp)
                    
                    active_ligand_target_links_df[,1] <- as.character(active_ligand_target_links_df[,1])
                    active_ligand_target_links_df[,2] <- as.character(active_ligand_target_links_df[,2])
                    active_ligand_target_links_df[,3] <- as.numeric(active_ligand_target_links_df[,3])
                  }
                  
                    #Should we just load DE results from EdgeR DE results that incorporate batch-effect correction instead of recalculating here?
                    DE_table_receiver = DE.table.input %>% filter(Reference == ref & Comparison == comp & CellType == receiver)
                    
                    geneset_oi = DE_table_receiver %>% filter(logFC >= 0.25 & padj <= 0.001) %>% arrange(desc(logFC)) %>% pull(Gene)
                    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
                    
                      #Define Potential Ligands for DE genes using previous reference networks:
                      ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
                      receptors = lr_network %>% dplyr::pull(receptor) %>% unique()
                      
                      expressed_ligands = intersect(ligands,expressed_genes_sender)
                      expressed_receptors = intersect(receptors,expressed_genes_receiver)
                      
                      potential_ligands = lr_network %>% filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>% dplyr::pull(ligand) %>% unique()
                      potential_ligands
                      
                      #Use Nichenet to estimate Ligand activities based on cell transcription profiles:
                      ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
                      
                      ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
                      
                      #Visualize top ranked ligands:
                    
                          #Should this be all ligands instead of top 40?
                          best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% dplyr::pull(test_ligand) %>% unique()
                          
                          DotPlot(seurat_obj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
                          
                          #Active Target Gene inference using receiver cell data:
                          active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
                          
                          active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.90)
                          
                          order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
                          order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
                          rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
                          colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
                          
                          vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
                          p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
                          p_ligand_target_network
                          
                    
                     Myeloid.nichenet.plot <-p_ligand_target_network
                     

                     DE.table.all <- DE.table.input %>% filter(Reference == ref & Comparison == comp)
                     DE.table.all <- DE.table.all[,which(colnames(DE.table.all) %in% c("Gene", "logFC", "CellType"))]            
                     DE.table.all <- DE.table.all %>% pivot_wider(values_from = "logFC", names_from = "CellType")
                     DE.table.all[is.na(DE.table.all)] = 0
                     
                     if (any(rowSums(DE.table.all[,-1]) == 0)) {
                       DE.table.all <- DE.table.all[-which(rowSums(DE.table.all[,-1]) == 0),]
                      }
                     
                     colnames(DE.table.all)[1] <- "gene"
                     
                     # Combine ligand activities with DE information
                     ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE.table.all %>% rename(ligand = gene))
                     ligand_activities_de[is.na(ligand_activities_de)] = 0
                     
                     # make LFC heatmap
                     lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
                     rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
                     
                     order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
                     vis_ligand_lfc = lfc_matrix[order_ligands,] %>% as.matrix()
                     
                     colnames(vis_ligand_lfc) = DE_table_all[,-1] %>% colnames() %>% make.names()
                     
                     p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
                     p_ligand_lfc
                     
                     #Summary combined Plots:
                             # ligand activity heatmap
                             ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
                             
                             rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
                             colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
                             
                             vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
                             p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
                             # ligand expression Seurat dotplot
                             order_ligands_adapted = order_ligands
                             #order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
                             #order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
                             
                             rotated_dotplot = DotPlot(seurat_obj %>% subset(Semi_Detailed_Celltype %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
                             figures_without_legend = cowplot::plot_grid(
                               p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                               rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
                               p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                               p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                               align = "hv",
                               nrow = 1,
                               rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))
                             
                             legends = cowplot::plot_grid(
                               ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
                               ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
                               ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
                               ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
                               nrow = 1,
                               align = "h", rel_widths = c(1.5, 1, 1, 1))
                             
                             combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
                             combined_plot
                                  
                             myeloid.combined_plot <- combined_plot
                          
                             #ggsave(myeloid.combined_plot, width = 60, height = 35,units = "cm",filename = "myeloid_D2allergy_nichenetplot.pdf", limitsize = FALSE)  
                             
                          
                          
                          
        Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
                     
        plot_violin_bar_srt(DE.input, gene = "CD274", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
        

        plot_barplot_srt(DE.input, gene = "CCL16", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
        
        
        plot_barplot_srt(DE.input, gene = "IL4", color_by = "Lesion", number_labels = F, colors = Lesion.colors, cell_proportion_weighted = FALSE)
        
        LTB.barplot <- last_plot()
        
        plot_violin_bar_srt(DE.input, gene = "ALOX15", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors, plot_patient_pos_frac = FALSE)
        plot_violin_bar_srt(DE.input, gene = "TIMP1", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors, plot_patient_pos_frac = FALSE)

        IL13.violin.barplot <- last_plot()
        
        setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")
        
        itch.genes <- c("IL31", "OSM", "LIF", "IL4", "IL13", "IL6")
        
        set.seed(100)
        itch.heatmap <- plot_heatmap_srt(DE.input, 
                                        genes = itch.genes, 
                                        type = "bulk", 
                                        cluster_by = "row", 
                                        pdf_format = "tile", 
                                        scale_by = "row",
                                        text_angle = 90,
                                        color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                        log_scale_base = 10,
                                        ceiling = 2,
                                        floor = -2)
        
        Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
        
        
        #DE Cytokines with lots of connections at bottom of network:
              plot_violin_bar_srt(DE.input, gene = "CCL2", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              CCL2.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CCL2", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CCL2.barplot <- last_plot()
              
              plot_violin_bar_srt(DE.input, gene = "CXCL9", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              CXCL9.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CXCL9", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CXCL9.barplot <- last_plot()
              
              plot_violin_bar_srt(DE.input, gene = "CXCL10", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              CXCL10.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CXCL10", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CXCL10.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "SEMA7A", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              SEMA7A.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "EREG", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              EREG.barplot <- last_plot()
              
              
              plot_barplot_srt(DE.input, gene = "AREG", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              AREG.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "VEGFB", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              VEGFB.barplot <- last_plot()
              
              
              
              plot_barplot_srt(DE.input, gene = "CD274", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CD274.barplot <- last_plot()
              
              
              plot_barplot_srt(DE.input, gene = "CXCL11", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CXCL11.barplot <- last_plot()
              
              
              
              plot_barplot_srt(DE.input, gene = "CXCL1", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CXCL1.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CXCL2", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CXCL2.barplot <- last_plot()
              
              
              plot_barplot_srt(DE.input, gene = "CXCL3", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CXCL3.barplot <- last_plot()
              
              
              plot_barplot_srt(DE.input, gene = "IFNG", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              IFNG.barplot <- last_plot()
              
              
              plot_barplot_srt(DE.input, gene = "IL13", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              IL13.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "IL4", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              IL4.barplot <- last_plot()
              
        #Trending Irritant Ligands not found by DE:
              plot_violin_bar_srt(DE.input, gene = "CD70", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              CD70.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CD70", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CD70.barplot <- last_plot()
        
              plot_violin_bar_srt(DE.input, gene = "TNFSF4", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              TNFSF4.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "TNFSF4", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              TNFSF4.barplot <- last_plot()
        
              plot_violin_bar_srt(DE.input, gene = "VEGFA", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              VEGFA.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "VEGFA", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              VEGFA.barplot <- last_plot()
              
              
        #Trending Allergy Ligands not found by DE:
              plot_violin_bar_srt(DE.input, gene = "IL4", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              IL4.violin.barplot <- last_plot()
              
              
              plot_violin_bar_srt(DE.input, gene = "LIF", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              LIF.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "LIF", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              LIF.barplot <- last_plot()
              
              plot_violin_bar_srt(DE.input, gene = "OSM", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              OSM.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "OSM", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              OSM.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CNTF", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CNTF.barplot <- last_plot()
              
              plot_violin_bar_srt(DE.input, gene = "LTA", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              LTA.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "LTA", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              LTA.barplot <- last_plot()
              
              plot_violin_bar_srt(DE.input, gene = "CCL7", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              CCL7.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CCL7", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CCL7.barplot <- last_plot()
              
              plot_violin_bar_srt(DE.input, gene = "CCL8", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              CCL8.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "CCL8", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CCL8.barplot <- last_plot()
              
              plot_violin_bar_srt(DE.input, gene = "CCL24", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
              CCL24.violin.barplot <- last_plot()
              
              plot_barplot_srt(DE.input, gene = "IL4R", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
              CCL24.barplot <- last_plot()
              
              
              
              
        setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")
        
        
        #ggsave(IFNG.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/IFNG.violin.barplot.png", limitsize = FALSE)  
        #ggsave(IL4.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/IL4.violin.barplot.png", limitsize = FALSE)  
        #ggsave(IL13.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/IL13.violin.barplot.png", limitsize = FALSE)  
        
              
        #ggsave(CCL2.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL2.violin.barplot.png", limitsize = FALSE)  
        #ggsave(CCL2.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL2.barplot.pdf", limitsize = FALSE)  
              
        #ggsave(CXCL9.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL9.violin.barplot.png", limitsize = FALSE)  
        #ggsave(CXCL9.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL9.barplot.pdf", limitsize = FALSE)  
      
        #ggsave(CXCL10.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL10.violin.barplot.png", limitsize = FALSE)  
        #ggsave(CXCL10.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL10.barplot.pdf", limitsize = FALSE)  
        
              
              
        #ggsave(CD70.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CD70.violin.barplot.png", limitsize = FALSE)  
        #ggsave(CD70.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CD70.barplot.pdf", limitsize = FALSE)  
      
              
        #ggsave(TNFSF4.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/TNFSF4.violin.barplot.png", limitsize = FALSE)  
        #ggsave(TNFSF4.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/TNFSF4.barplot.pdf", limitsize = FALSE)  
        
        #ggsave(VEGFA.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/VEGFA.violin.barplot.png", limitsize = FALSE)  
        #ggsave(VEGFA.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/VEGFA.barplot.pdf", limitsize = FALSE)  
        
              
        
        #ggsave(LIF.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/LIF.violin.barplot.png", limitsize = FALSE)  
        #ggsave(LIF.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/LIF.barplot.pdf", limitsize = FALSE)  
              
        #ggsave(OSM.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/OSM.violin.barplot.png", limitsize = FALSE)  
        #ggsave(OSM.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/OSM.barplot.pdf", limitsize = FALSE)  
        
        #ggsave(CNTF.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CNTF.barplot.pdf", limitsize = FALSE)  
              
        #ggsave(LTA.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/LTA.violin.barplot.png", limitsize = FALSE)  
        #ggsave(LTA.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/LTA.barplot.pdf", limitsize = FALSE)  
              
        #ggsave(CCL7.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL7.violin.barplot.png", limitsize = FALSE)  
        #ggsave(CCL7.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL7.barplot.pdf", limitsize = FALSE)  

        #ggsave(CCL8.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL8.violin.barplot.png", limitsize = FALSE)  
        #ggsave(CCL8.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL8.barplot.pdf", limitsize = FALSE)  
        
        #ggsave(CCL24.violin.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL24.violin.barplot.png", limitsize = FALSE)  
        #ggsave(CCL24.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CCL24.barplot.pdf", limitsize = FALSE)  
        
        #ggsave(SEMA7A.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/SEMA7A.barplot.pdf", limitsize = FALSE)  
        #ggsave(AREG.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/AREG.barplot.pdf", limitsize = FALSE)  
        #ggsave(EREG.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/EREG.barplot.pdf", limitsize = FALSE)  
        #ggsave(CD274.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CD274.barplot.pdf", limitsize = FALSE)  
        #ggsave(CXCL11.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL11.barplot.pdf", limitsize = FALSE)  
        #ggsave(VEGFB.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/VEGFB.barplot.pdf", limitsize = FALSE)  
        #ggsave(CXCL1.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL1.barplot.pdf", limitsize = FALSE)  
        #ggsave(CXCL2.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL2.barplot.pdf", limitsize = FALSE)  
        #ggsave(CXCL3.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/CXCL3.barplot.pdf", limitsize = FALSE)  
        #ggsave(IFNG.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/IFNG.barplot.pdf", limitsize = FALSE)  
        #ggsave(IL4.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/IL4.barplot.pdf", limitsize = FALSE)  
        #ggsave(IL13.barplot, width = 60, height = 32,units = "cm",filename = "Plots/barplots/IL13.barplot.pdf", limitsize = FALSE)  
        
              
                
        #edges <- network_2lesion.result$Graph_object[[1]] %>% activate(edges) %>% as.data.frame()
        
        
        #Receptor ligand heatmaps:
                    ligs <- network %>% dplyr::pull(Ligand) %>% unique()
                    
                    complex.ligs <- deconvoluted %>% filter(complex_name %in% ligs) %>% dplyr::pull(complex_name) %>% unique()
                    complex.lig.genes <- deconvoluted %>% filter(complex_name %in% ligs) %>% dplyr::pull(gene_name) %>% unique()
                    
                    if (length(complex.ligs)>0) {
                      ligs <- ligs[-which(ligs %in% complex.ligs)]
                      ligs <- unique(c(ligs, complex.lig.genes))
                    }
                    
                    recs <- network %>% dplyr::pull(Receptor) %>% unique()
                    
                    complex.recs <- deconvoluted %>% filter(complex_name %in% recs) %>% dplyr::pull(complex_name) %>% unique()
                    complex.rec.genes <- deconvoluted %>% filter(complex_name %in% recs) %>% dplyr::pull(gene_name) %>% unique()
                    
                    if (length(complex.recs)>0) {
                      recs <- recs[-which(recs %in% complex.recs)]
                      recs <- unique(c(recs, complex.rec.genes))
                    }
                    
                    #By Lesion and Celltype:
                    DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Semi_Detailed_Celltype", "Lesion"))
                    
                          #Filter to DE ligs and receptors:
                            contact_derm_DE_filtered <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE_filtered.txt")
                            contact_derm_DE <- read.table("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_DE.txt")
                            
                          DE.ligs <- c()
                          for (l in 1:length(ligs)) {
                            lig <- ligs[l]

                            lig.data <- contact_derm_DE_filtered %>% filter(Gene == lig) %>% 
                                                          filter(Reference %in% c("Irritant/Day2_Allergy/Day4_Allergy", "Nonlesional") & Color_factor == "Significantly Up" & Comparison != "Acetone_Vehicle" & logFC > 1) 
                            
                            if (nrow(lig.data) > 0) {
                              DE.ligs <- c(DE.ligs, lig)
                            }
                            rm(lig.data)
                          }
                          
                          DE.recs <- c()
                          for (r in 1:length(recs)) {
                            rec <- recs[r]
                            
                            rec.data <- contact_derm_DE_filtered %>% filter(Gene == rec) %>% 
                              filter(Reference %in% c("Irritant/Day2_Allergy/Day4_Allergy", "Nonlesional") & Color_factor == "Significantly Up" & Comparison != "Acetone_Vehicle" & logFC > 1) 
                            
                            if (nrow(rec.data) > 0) {
                              DE.recs <- c(DE.recs, rec)
                            }
                            rm(rec.data)
                          }
                          
                         #Ligand Heatmaps:
                          
                                                                #Scaled by Celltype:
                                                                      set.seed(100)
                                                                      lig.heatmap <- plot_heatmap_srt(DE.input, 
                                                                                                      genes = ligs, 
                                                                                                      type = "bulk", 
                                                                                                      facet_by = "Lesion",
                                                                                                      scale_group = "Semi_Detailed_Celltype",
                                                                                                      cluster_by = "row", 
                                                                                                      pdf_format = "tile", 
                                                                                                      scale_by = "row",
                                                                                                      text_angle = 90,
                                                                                                      cluster_type = "kmeans", 
                                                                                                      k = round(length(ligs)^0.5, digits =0), 
                                                                                                      show_k = T,
                                                                                                      color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                                                                      log_scale_base = 10,
                                                                                                      ceiling = F,
                                                                                                      floor = F)
                                                                      
                                                                      
                                                                      #Gene order:
                                                                      #lig_genes <- c(7,4,     6,8)
                                                                      lig_genes <- c(11,10,2,8,      4,17,18)
                                                                      lig_genes <- rev(lig_genes)
                                                                      lig_genes_reorder <- c()
                                                                      
                                                                      groups.to.not.plot <-c(1,3,5,6,7,9,12,13,14,15,16)
                                                                      
                                                                      for (i in 1:length(lig_genes)) {
                                                                        int1 <- lig_genes[i]
                                                                        ind <- grep(paste0("^", int1, "$"), lig.heatmap[[2]]$cluster)
                                                                        reorder <- names(lig.heatmap[[2]]$cluster)[ind]
                                                                        lig_genes_reorder <- c(lig_genes_reorder, reorder)
                                                                      }
                                                                      
                                                                      #Heatmap Gene breaks:
                                                                      #lig_cluster_breaks <- c(6)
                                                                      lig_cluster_breaks <- c(4)
                                                                      
                                                                      lig_cluster_breaks <- rev(lig_cluster_breaks)
                                                                      
                                                                      lig_gene_breaks <- find_gene_breaks(lig.heatmap, lig_genes_reorder,lig_cluster_breaks)
                                                                      
                                                                      
                                                                      #Get full list of DE genes that I filtered out, to include in supplemental:
                                                                      genes.to.not.plot <- c()
                                                                      for (i in 1:length(groups.to.not.plot)) {
                                                                        int1 <- groups.to.not.plot[i]
                                                                        ind <- grep(paste0("^", int1, "$"), lig.heatmap[[2]]$cluster)
                                                                        reorder <- names(lig.heatmap[[2]]$cluster)[ind]
                                                                        genes.to.not.plot <- c(genes.to.not.plot, reorder)
                                                                      }         
                                                                      
                                                                      
                                                                      plot_heatmap_srt(DE.input, 
                                                                                       genes = lig_genes_reorder, 
                                                                                       gene_breaks = lig_gene_breaks,
                                                                                       gene_names = TRUE, 
                                                                                       gene_labels = lig_genes_reorder[which(str_detect(lig_genes_reorder, pattern = "^IL|^CXCL|^CCL|^TNF|^LT|^IFN|^OSM"))],
                                                                                       gene_label_side = "left",
                                                                                       gene_labels_size = 3,
                                                                                       gene_labels_nudge = 0.3,
                                                                                       type = "bulk", 
                                                                                       facet_by = "Lesion",
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
                                                                                       plot_axis = FALSE)
                                                                      
                                                                      #plot_violin_bar_srt(DE.input, gene = "LTC4S", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
                                                                      
                                                                      ligand.plot.celltype.scaled <- last_plot()
                                                                      
                                                                      plot_heatmap_srt(DE.input, 
                                                                                       genes = c(lig_genes_reorder,genes.to.not.plot),
                                                                                       gene_breaks =  c(lig_gene_breaks,lig_genes_reorder[length(lig_genes_reorder)]),
                                                                                       gene_names = TRUE, 
                                                                                       gene_labels = NULL,
                                                                                       gene_label_side = "left",
                                                                                       gene_labels_size = 3,
                                                                                       gene_labels_nudge = 0.3,
                                                                                       type = "bulk", 
                                                                                       facet_by = "Lesion",
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
                                                                                       plot_axis = FALSE)
                                                                      
                                                                      ligand.plot.celltype.scaled.full <- last_plot()
                                                                      
                                      #Not scaled by Celltype:
                                          set.seed(100)
                                          lig.heatmap <- plot_heatmap_srt(DE.input, 
                                                                          genes = DE.ligs, 
                                                                          type = "bulk", 
                                                                          facet_by = "Semi_Detailed_Celltype",
                                                                          cluster_by = "row", 
                                                                          pdf_format = "tile", 
                                                                          scale_by = "row",
                                                                          text_angle = 90,
                                                                          cluster_type = "kmeans", 
                                                                          k = round(length(DE.ligs)^0.5, digits =0), 
                                                                          show_k = T,
                                                                          color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                                          log_scale_base = 10,
                                                                          ceiling = 2,
                                                                          floor = -0.1)
                                          
                                          
                                          #Gene order:
                                          #lig_genes <- c(5,   12,8,11,   7,9,  4,  1,    3,    6,  2, 10  )
                                          lig_genes <- c(1,    5,    2,   3,   6,    4,7)
                                          
                                          lig_genes <- rev(lig_genes)
                                          lig_genes_reorder <- c()
                                          
                                          for (i in 1:length(lig_genes)) {
                                            int1 <- lig_genes[i]
                                            ind <- grep(paste0("^", int1, "$"), lig.heatmap[[2]]$cluster)
                                            reorder <- names(lig.heatmap[[2]]$cluster)[ind]
                                            lig_genes_reorder <- c(lig_genes_reorder, reorder)
                                          }
                                          
                                          #Heatmap Gene breaks:
                                          #lig_cluster_breaks <- c(12,7,4,1,3,6)
                                          lig_cluster_breaks <- c(5,2,3,6,4)
                                          
                                          lig_cluster_breaks <- rev(lig_cluster_breaks)
                                          
                                          lig_gene_breaks <- find_gene_breaks(lig.heatmap, lig_genes_reorder,lig_cluster_breaks)
                                          
                                          DE.ligs.reordered <- lig_genes_reorder
                                          DE.lig.breaks <- lig_gene_breaks
                                          
                                          plot_heatmap_srt(DE.input, 
                                                           genes = DE.ligs.reordered, 
                                                           gene_breaks = DE.lig.breaks,
                                                           gene_names = TRUE, 
                                                           gene_labels = DE.ligs.reordered[which(str_detect(DE.ligs.reordered, pattern = "^IL|^CCL|^CXCL|^IFN|^TNF|^LT|^LIF|^OSM"))],
                                                           gene_label_side = "left",
                                                           gene_labels_size = 3,
                                                           gene_labels_nudge = 0.3,
                                                           type = "bulk", 
                                                           facet_by = "Semi_Detailed_Celltype",
                                                           cluster_by = FALSE, 
                                                           pdf_format = "tile", 
                                                           scale_by = "row",
                                                           text_angle = 90,
                                                           color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                           log_scale_base = 10,
                                                           ceiling = 2,
                                                           floor = -0.1,
                                                           panel_spacing = 0.1,
                                                           plot_axis = FALSE)
                                          
                                          ligand.plot.DE <- last_plot()
                                          
                                          #plot_violin_bar_srt(DE.input, gene = "IL4", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
                                          
                                          #All Ligs:
                                          set.seed(100)
                                          lig.heatmap <- plot_heatmap_srt(DE.input, 
                                                                          genes = ligs[-which(ligs %in% DE.ligs)], 
                                                                          type = "bulk", 
                                                                          facet_by = "Semi_Detailed_Celltype",
                                                                          cluster_by = "row", 
                                                                          pdf_format = "tile", 
                                                                          scale_by = "row",
                                                                          text_angle = 90,
                                                                          cluster_type = "kmeans", 
                                                                          k = round(length(ligs[-which(ligs %in% DE.ligs)])^0.5, digits =0), 
                                                                          show_k = T,
                                                                          color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                                          log_scale_base = 10,
                                                                          ceiling = 2,
                                                                          floor = -0.1)
                                          
                                          
                                          #Gene order:
                                          lig_genes <- c(5,9,10,   7,1,14,  12,8,  13,  17,  4,11,16,   6,  15,2,3)
                                          
                                          lig_genes <- rev(lig_genes)
                                          lig_genes_reorder <- c()
                                          
                                          for (i in 1:length(lig_genes)) {
                                            int1 <- lig_genes[i]
                                            ind <- grep(paste0("^", int1, "$"), lig.heatmap[[2]]$cluster)
                                            reorder <- names(lig.heatmap[[2]]$cluster)[ind]
                                            lig_genes_reorder <- c(lig_genes_reorder, reorder)
                                          }
                                          
                                          #Heatmap Gene breaks:
                                          lig_cluster_breaks <- c(7,12,13,17,4,6,15)
                                          
                                          lig_cluster_breaks <- rev(lig_cluster_breaks)
                                          
                                          lig_gene_breaks <- find_gene_breaks(lig.heatmap, lig_genes_reorder,lig_cluster_breaks)
                                          
                                          plot_heatmap_srt(DE.input, 
                                                           genes = c(DE.ligs.reordered,lig_genes_reorder), 
                                                           gene_breaks = c(DE.lig.breaks,DE.ligs.reordered[length(DE.ligs.reordered)],lig_gene_breaks),
                                                           gene_names = TRUE, 
                                                           gene_labels =  c(DE.ligs.reordered,lig_genes_reorder)[which(str_detect(c(DE.ligs.reordered,lig_genes_reorder), pattern = "^IL|^CCL|^CXCL|^IFN|^TNF|^LT|^LIF|^OSM"))],
                                                           gene_label_side = "left",
                                                           gene_labels_size = 3,
                                                           gene_labels_nudge = 0.3,
                                                           type = "bulk", 
                                                           facet_by = "Semi_Detailed_Celltype",
                                                           cluster_by = FALSE, 
                                                           pdf_format = "tile", 
                                                           scale_by = "row",
                                                           text_angle = 90,
                                                           color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                           log_scale_base = 10,
                                                           ceiling = 2,
                                                           floor = -0.1,
                                                           panel_spacing = 0.1,
                                                           plot_axis = FALSE)
                                          
                                          ligand.plot.full <- last_plot()
                                         
                                          #Manual editing to merge DE ligs with Non-DE ligs that look like they have a trend for allergy:
                                          trend.genes <- lig_genes_reorder[c(1:which(lig_genes_reorder==lig_gene_breaks[2]))]
                                          trend.gene.breaks <- lig_gene_breaks[1:2]
                                          trend.ligs <- trend.genes
                                          
                                          plot_heatmap_srt(DE.input, 
                                                           genes = c(trend.genes, DE.ligs.reordered), 
                                                           gene_breaks = c(trend.gene.breaks,DE.lig.breaks),
                                                           gene_names = TRUE, 
                                                           gene_labels = c(trend.genes, DE.ligs.reordered)[which(str_detect(c(trend.genes, DE.ligs.reordered), pattern = "^IL|^CCL|^CXCL|^IFN|^TNF|^LT|^LIF|^OSM"))],
                                                           gene_label_side = "left",
                                                           gene_labels_size = 3,
                                                           gene_labels_nudge = 0.3,
                                                           type = "bulk", 
                                                           facet_by = "Semi_Detailed_Celltype",
                                                           cluster_by = FALSE, 
                                                           pdf_format = "tile", 
                                                           scale_by = "row",
                                                           text_angle = 90,
                                                           color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                           log_scale_base = 10,
                                                           ceiling = 2,
                                                           floor = -0.1,
                                                           panel_spacing = 0.1,
                                                           plot_axis = FALSE)
                                          
                                          ligand.plot.DE.plus.trending <- last_plot()
                                          
                                          
                    #Receptor Heatmaps:
                    
                                                                              #Scaled by Celltype:                    
                                                                                  set.seed(100)
                                                                                  rec.heatmap <- plot_heatmap_srt(DE.input, 
                                                                                                                  genes = recs, 
                                                                                                                  type = "bulk", 
                                                                                                                  facet_by = "Lesion",
                                                                                                                  scale_group = "Semi_Detailed_Celltype",
                                                                                                                  cluster_by = "row", 
                                                                                                                  pdf_format = "tile", 
                                                                                                                  scale_by = "row",
                                                                                                                  text_angle = 90,
                                                                                                                  cluster_type = "kmeans", 
                                                                                                                  k = round(length(ligs)^0.5, digits =0), 
                                                                                                                  show_k = T,
                                                                                                                  color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                                                                                  log_scale_base = 10,
                                                                                                                  ceiling = F,
                                                                                                                  floor = F)
                                                                                  
                                                                                  #Gene order:
                                                                                  #rec_genes <- c(1,7,   8)
                                                                                  rec_genes <- c(6,8,10,13,14,15,18,  2,      7,16)
                                                                                  
                                                                                  rec_genes <- rev(rec_genes) 
                                                                                  rec_genes_reorder <- c()
                                                                                  
                                                                                  groups.to.not.plot <-c(1,3,4,5,9,11,12,17)
                                                                                  
                                                                                  for (i in 1:length(rec_genes)) {
                                                                                    int1 <- rec_genes[i]
                                                                                    ind <- grep(paste0("^", int1, "$"), rec.heatmap[[2]]$cluster)
                                                                                    reorder <- names(rec.heatmap[[2]]$cluster)[ind]
                                                                                    rec_genes_reorder <- c(rec_genes_reorder, reorder)
                                                                                  }
                                                                                  
                                                                                  #Heatmap Gene breaks:
                                                                                  #rec_cluster_breaks <- c(8)
                                                                                  rec_cluster_breaks <- c(7)
                                                                                  
                                                                                  rec_cluster_breaks <- rev(rec_cluster_breaks)
                                                                                  
                                                                                  rec_gene_breaks <- find_gene_breaks(rec.heatmap, rec_genes_reorder,rec_cluster_breaks)
                                                                                  
                                                                                  
                                                                                  #Get full list of DE genes that I filtered out, to include in supplemental:
                                                                                  genes.to.not.plot <- c()
                                                                                  for (i in 1:length(groups.to.not.plot)) {
                                                                                    int1 <- groups.to.not.plot[i]
                                                                                    ind <- grep(paste0("^", int1, "$"), rec.heatmap[[2]]$cluster)
                                                                                    reorder <- names(rec.heatmap[[2]]$cluster)[ind]
                                                                                    genes.to.not.plot <- c(genes.to.not.plot, reorder)
                                                                                  }         
                                                                                  
                                                                                  
                                                                                  plot_heatmap_srt(DE.input, 
                                                                                                   genes = rec_genes_reorder, 
                                                                                                   gene_breaks = rec_gene_breaks,
                                                                                                   gene_names = TRUE, 
                                                                                                   gene_labels = rec_genes_reorder[which(str_detect(rec_genes_reorder, pattern = "^IL|^CCL|^CXCL|^IFN|^TNF|^LT|^LIF|^OSM"))],
                                                                                                   gene_label_side = "left",
                                                                                                   gene_labels_size = 3,
                                                                                                   gene_labels_nudge = 0.3,
                                                                                                   type = "bulk", 
                                                                                                   facet_by = "Lesion",
                                                                                                   scale_group = "Semi_Detailed_Celltype",
                                                                                                   cluster_by = FALSE, 
                                                                                                   pdf_format = "tile", 
                                                                                                   scale_by = "row",
                                                                                                   text_angle = 90,
                                                                                                   color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#000000","#000000","#792B45","#EC1F63"))(256),
                                                                                                   log_scale_base = 10,
                                                                                                   ceiling = F,
                                                                                                   floor = F,
                                                                                                   panel_spacing = 1,
                                                                                                   plot_axis = FALSE)
                                                                                  
                                                                                  receptor.plot.celltype.scaled <- last_plot()
                                                                                  
                                                                                  #plot_violin_bar_srt(DE.input, gene = "CXCL13", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
                                                                                  
                                                                                  
                                                                                  plot_heatmap_srt(DE.input, 
                                                                                                   genes = c(rec_genes_reorder,genes.to.not.plot),
                                                                                                   gene_breaks =  c(rec_gene_breaks,rec_genes_reorder[length(rec_genes_reorder)]),
                                                                                                   gene_names = TRUE, 
                                                                                                   gene_labels = NULL,
                                                                                                   gene_label_side = "left",
                                                                                                   gene_labels_size = 3,
                                                                                                   gene_labels_nudge = 0.3,
                                                                                                   type = "bulk", 
                                                                                                   facet_by = "Lesion",
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
                                                                                                   plot_axis = FALSE)
                                                                                  
                                                                                  receptor.plot.celltype.scaled.full <- last_plot()
                                                                                  
                                  #Not Scaled by Celltype:                    
                                        set.seed(100)
                                        rec.heatmap <- plot_heatmap_srt(DE.input, 
                                                                        genes = DE.recs, 
                                                                        type = "bulk", 
                                                                        facet_by = "Semi_Detailed_Celltype",
                                                                        cluster_by = "row", 
                                                                        pdf_format = "tile", 
                                                                        scale_by = "row",
                                                                        text_angle = 90,
                                                                        cluster_type = "kmeans", 
                                                                        k = round(length(DE.recs)^0.5, digits =0), 
                                                                        show_k = T,
                                                                        color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                                        log_scale_base = 10,
                                                                        ceiling = 2,
                                                                        floor = -0.1)
                                        
                                        #Gene order:
                                        #rec_genes <- c(12,8,9,    3,    2,6,    1,5,7,   4,10,11)
                                        rec_genes <- c(1,  7,   4,    6,3,5,2)
                                        rec_genes <- rev(rec_genes)
                                        rec_genes_reorder <- c()
                                        
                                        for (i in 1:length(rec_genes)) {
                                          int1 <- rec_genes[i]
                                          ind <- grep(paste0("^", int1, "$"), rec.heatmap[[2]]$cluster)
                                          reorder <- names(rec.heatmap[[2]]$cluster)[ind]
                                          rec_genes_reorder <- c(rec_genes_reorder, reorder)
                                        }
                                        
                                        #Heatmap Gene breaks:
                                        #rec_cluster_breaks <- c(3,2,1,4)
                                        rec_cluster_breaks <- c(7,4,6)
                                        
                                        rec_cluster_breaks <- rev(rec_cluster_breaks)
                                        
                                        rec_gene_breaks <- find_gene_breaks(rec.heatmap, rec_genes_reorder,rec_cluster_breaks)
                                        
                                        DE.recs.reordered <- rec_genes_reorder
                                        DE.rec.breaks <- rec_gene_breaks
                                        
                                        plot_heatmap_srt(DE.input, 
                                                         genes = DE.recs.reordered, 
                                                         gene_breaks = DE.rec.breaks,
                                                         gene_names = TRUE, 
                                                         gene_labels = DE.recs.reordered[which(str_detect(DE.recs.reordered, pattern = "^IL|^CCR|^CXCR|^IFN|^TNF|^LT|^LIF|^OSM"))],
                                                         gene_label_side = "left",
                                                         gene_labels_size = 3,
                                                         gene_labels_nudge = 0.3,
                                                         type = "bulk", 
                                                         facet_by = "Semi_Detailed_Celltype",
                                                         cluster_by = FALSE, 
                                                         pdf_format = "tile", 
                                                         scale_by = "row",
                                                         text_angle = 90,
                                                         color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                         log_scale_base = 10,
                                                         ceiling = 2,
                                                         floor = -0.1,
                                                         panel_spacing = 0.1,
                                                         plot_axis = FALSE)
                                        
                                        receptor.plot.DE <- last_plot()
                                        
                                      #All Receptors:
                                        set.seed(100)
                                        rec.heatmap <- plot_heatmap_srt(DE.input, 
                                                                        genes = recs[-which(recs %in% DE.recs)], 
                                                                        type = "bulk", 
                                                                        facet_by = "Semi_Detailed_Celltype",
                                                                        cluster_by = "row", 
                                                                        pdf_format = "tile", 
                                                                        scale_by = "row",
                                                                        text_angle = 90,
                                                                        cluster_type = "kmeans", 
                                                                        k = round(length(recs[-which(recs %in% DE.recs)])^0.5, digits =0), 
                                                                        show_k = T,
                                                                        color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                                        log_scale_base = 10,
                                                                        ceiling = 2,
                                                                        floor = -0.1)
                                        
                                        
                                        #Gene order:
                                        rec_genes <- c(5,14,13,     3,16,   6,   10,  7,4,  2,  8,9,   1,    12,15,   11)
                                        
                                        rec_genes <- rev(rec_genes)
                                        rec_genes_reorder <- c()
                                        
                                        for (i in 1:length(rec_genes)) {
                                          int1 <- rec_genes[i]
                                          ind <- grep(paste0("^", int1, "$"), rec.heatmap[[2]]$cluster)
                                          reorder <- names(rec.heatmap[[2]]$cluster)[ind]
                                          rec_genes_reorder <- c(rec_genes_reorder, reorder)
                                        }
                                        
                                        #Heatmap Gene breaks:
                                        rec_cluster_breaks <- c(3,6,10,7,2,8,1,12,11)
                                        
                                        rec_cluster_breaks <- rev(rec_cluster_breaks)
                                        
                                        rec_gene_breaks <- find_gene_breaks(rec.heatmap, rec_genes_reorder,rec_cluster_breaks)
                                        
                                        plot_heatmap_srt(DE.input, 
                                                         genes = c(DE.recs.reordered,rec_genes_reorder), 
                                                         gene_breaks = c(DE.rec.breaks,DE.recs.reordered[length(DE.recs.reordered)],rec_gene_breaks),
                                                         gene_names = TRUE, 
                                                         gene_labels = c(DE.recs.reordered,rec_genes_reorder)[which(str_detect(c(DE.recs.reordered,rec_genes_reorder), pattern = "^IL|^CCR|^CXCR|^IFN|^TNF|^LT|^LIF|^OSM"))],
                                                         gene_label_side = "left",
                                                         gene_labels_size = 3,
                                                         gene_labels_nudge = 0.3,
                                                         type = "bulk", 
                                                         facet_by = "Semi_Detailed_Celltype",
                                                         cluster_by = FALSE, 
                                                         pdf_format = "tile", 
                                                         scale_by = "row",
                                                         text_angle = 90,
                                                         color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                         log_scale_base = 10,
                                                         ceiling = 2,
                                                         floor = -0.1,
                                                         panel_spacing = 0.1,
                                                         plot_axis = FALSE)
                                        
                                        receptor.plot.full <- last_plot()
                                        
                                        
                                        #Manual editing to merge DE recs with Non-DE recs that look like they have a trend for allergy:
                                        trend.genes <- rec_genes_reorder[c(1:which(rec_genes_reorder==rec_gene_breaks[2]))]
                                        trend.gene.breaks <- rec_gene_breaks[1:2]
                                        trend.recs <- trend.genes
                                        
                                        plot_heatmap_srt(DE.input, 
                                                         genes = c(trend.genes, DE.recs.reordered), 
                                                         gene_breaks = c(trend.gene.breaks,DE.rec.breaks),
                                                         gene_names = TRUE, 
                                                         gene_labels = c(trend.genes, DE.recs.reordered)[which(str_detect(c(trend.genes, DE.recs.reordered), pattern = "^IL|^CCR|^CXCR|^IFN|^TNF|^LT|^LIF|^OSM"))],
                                                         gene_label_side = "left",
                                                         gene_labels_size = 3,
                                                         gene_labels_nudge = 0.3,
                                                         type = "bulk", 
                                                         facet_by = "Semi_Detailed_Celltype",
                                                         cluster_by = FALSE, 
                                                         pdf_format = "tile", 
                                                         scale_by = "row",
                                                         text_angle = 90,
                                                         color_pal = colorRampPalette(c("#000000" ,"#FFB326","#EC1F63"))(256),
                                                         log_scale_base = 10,
                                                         ceiling = 2,
                                                         floor = -0.1,
                                                         panel_spacing = 0.1,
                                                         plot_axis = FALSE)
                                        
                                        receptor.plot.DE.plus.trending <- last_plot()
                                        
                                              
                #setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/rl_network")
                
                #ggsave(ligand.plot.DE, width = 80, height = 80,units = "cm",filename = "ligand.plot.DE.pdf", limitsize = FALSE)  
                #ggsave(ligand.plot.DE.plus.trending, width = 80, height = 80,units = "cm",filename = "ligand.plot.DE.plus.trending.pdf", limitsize = FALSE)  
                #ggsave(ligand.plot.full, width = 80, height = 80,units = "cm",filename = "ligand.plot.full.pdf", limitsize = FALSE)  
                
                #ggsave(receptor.plot.DE, width = 80, height = 80,units = "cm",filename = "receptor.plot.DE.pdf", limitsize = FALSE)  
                #ggsave(receptor.plot.DE.plus.trending, width = 80, height = 80,units = "cm",filename = "receptor.plot.DE.plus.trending.pdf", limitsize = FALSE)  
                #ggsave(receptor.plot.full, width = 80, height = 80,units = "cm",filename = "receptor.plot.full.pdf", limitsize = FALSE)  
                
      #DE vs Trending Gene Stats Comparison for mean UMIs per cell and Number of cells expressing gene:
            
        RL.stats <- data.frame(Gene = c(ligs,recs),
                               Type = c(rep("Ligand", times = length(ligs)), rep("Receptor", times = length(recs))),
                               Classification = NA,
                               Positive.cell.count = NA,
                               Mean.positive.raw.UMIs = NA)
        
        RL.stats <- RL.stats %>% mutate(Classification = case_when(Gene %in% DE.ligs & Type == "Ligand" ~ "DE Ligand",
                                                                   Gene %in% trend.ligs & Type == "Ligand" ~ "Trend Ligand",
                                                                   Gene %in% DE.recs & Type == "Receptor" ~ "DE Receptor",
                                                                   Gene %in% trend.recs & Type == "Receptor" ~ "Trend Receptor",
                                                                   Gene %in% ligs & Type == "Ligand" ~ "Unclassified Ligand",
                                                                   Gene %in% recs & Type == "Receptor"~ "Unclassified Receptor"))
        
        RL.stats$Classification <- factor(RL.stats$Classification,
                                          levels = c("DE Ligand", "Trend Ligand", "Unclassified Ligand", "DE Receptor", "Trend Receptor", "Unclassified Receptor"))
      
        for (r in 1:nrow(RL.stats)) {
          
          gene = RL.stats$Gene[r]
          
          gene.val.string <- as.numeric(DE.input@assays$RNA[which(rownames(DE.input@assays$RNA) == gene),])
          
          RL.stats[r, "Positive.cell.count"] <- as.numeric(length(gene.val.string[which(gene.val.string>0)]))
          RL.stats[r, "Mean.positive.raw.UMIs"] <- as.numeric(mean(gene.val.string[which(gene.val.string>0)]))
          
        }
        
        RL.stats %>% filter(Type == "Ligand") %>% 
          ggplot(., aes(x = log10(Positive.cell.count), color = Classification)) +
          geom_density(size =2)+
          theme_classic() +
          scale_color_manual(values = c("red", "orange", "black")) +
          labs(title = "Ligand Positive Expression Cell Counts by Gene Classification", x = "Log10 (Positive Cell Count)") +
          theme(plot.title = element_text(hjust = 0.5))
        
        lig.cell.count.density <- last_plot()
        
        RL.stats %>% filter(Type == "Ligand") %>% 
          ggplot(., aes(x = log10(Mean.positive.raw.UMIs), color = Classification)) +
          geom_density(size =2)+
          theme_classic() +
          scale_color_manual(values = c("red", "orange", "black")) +
          labs(title = "Ligand Mean Positive UMI Count by Gene Classification", x = "Log10 (Mean Positive Raw UMI Count)") +
          theme(plot.title = element_text(hjust = 0.5))
        
        lig.mean.umi.density <- last_plot()
        
        lig.cell.count.density + lig.mean.umi.density
        
        lig.stat.density <- last_plot()
        
        RL.stats %>% filter(Type == "Receptor") %>% 
          ggplot(., aes(x = log10(Positive.cell.count), color = Classification)) +
          geom_density(size =2) +
          theme_classic() +
          scale_color_manual(values = c("red", "orange", "black")) +
          labs(title = "Receptor Positive Expression Cell Counts by Gene Classification", x = "Log10 (Positive Cell Count)") +
          theme(plot.title = element_text(hjust = 0.5))
        
        rec.cell.count.density <- last_plot()
        
        RL.stats %>% filter(Type == "Receptor") %>% 
          ggplot(., aes(x = log10(Mean.positive.raw.UMIs), color = Classification)) +
          geom_density(size =2)+
          theme_classic() +
          scale_color_manual(values = c("red", "orange", "black")) +
          labs(title = "Receptor Mean Positive UMI Count by Gene Classification", x = "Log10 (Mean Positive Raw UMI Count)") +
          theme(plot.title = element_text(hjust = 0.5))
        
        rec.mean.umi.density <- last_plot()
        
        rec.cell.count.density + rec.mean.umi.density
        
        rec.stat.density <- last_plot()
                
        #setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/rl_network")
        #ggsave(lig.stat.density, width = 40, height = 20,units = "cm",filename = "ligand.stat.density.plots.pdf", limitsize = FALSE)  
        #ggsave(rec.stat.density, width = 40, height = 20,units = "cm",filename = "receptor.stat.density.plots.pdf", limitsize = FALSE)  
        
        
        
      #Ligand-Receptor CellphoneDB Heatmap:
          
          CellphoneDB.network <- network %>% filter(database_source == "CellphoneDB")
          
              #Filter to DE genes:
              DE.genes <- DE.table %>% dplyr::pull(Gene) %>% unique()
          
              DE.complexes <- deconvoluted %>% filter(gene_name %in% DE.genes) %>% dplyr::pull(complex_name) %>% unique()
          
              
              DE.CellphoneDB.network <- CellphoneDB.network %>% filter(Ligand %in% c(DE.genes,DE.complexes) | Receptor %in% c(DE.genes,DE.complexes))
          
          
          #Scale by fraction of max expression:
          max.expression.vals <- DE.CellphoneDB.network %>% group_by(Ligand, Receptor) %>% summarize(max.expression = max(log10_Connection_product))
          
          CellphoneDB.network_scaled <- DE.CellphoneDB.network
          for (r in 1:nrow(max.expression.vals)) {
            
            connection.row <- max.expression.vals[r,]
            connection.max <- max.expression.vals[r,"max.expression"]
            
            connection.edges <- DE.CellphoneDB.network %>% filter(Ligand == connection.row$Ligand & Receptor == connection.row$Receptor)
            connection.edges$log10_Connection_product <- sapply(connection.edges$log10_Connection_product, FUN = function(x) as.numeric(round(10^(x)/(10^(connection.max)), digits = 3)))
            
            CellphoneDB.network_scaled[which(CellphoneDB.network_scaled$edgeID %in% connection.edges$edgeID),] <- connection.edges
          }
          
          DE.CellphoneDB.network <- CellphoneDB.network_scaled
          
          
          DE.CellphoneDB.network <- DE.CellphoneDB.network %>% pivot_wider(names_from = "Lesion", values_from = "log10_Connection_product")
          DE.CellphoneDB.network$edge_ID <- paste(DE.CellphoneDB.network$from, "_to_", DE.CellphoneDB.network$to, sep = "")
          
          cellphonedb.heat.data <- DE.CellphoneDB.network[,which(colnames(DE.CellphoneDB.network) %in% c("edge_ID", "Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))]
          cellphonedb.heat.data <- cellphonedb.heat.data[-which(duplicated(cellphonedb.heat.data$edge_ID)),]
          cellphonedb.heat.data <- column_to_rownames(cellphonedb.heat.data, var = "edge_ID")
          cellphonedb.heat.data[is.na(cellphonedb.heat.data)] <- 0
          
          
          cellphonedb.heat.data <- cellphonedb.heat.data[,c(2,1, 5,3,4)]
          
          #Scale by row:
          cellphonedb.heat.data2 <- t(apply(cellphonedb.heat.data, 1, scale))
          colnames(cellphonedb.heat.data2) <- colnames(cellphonedb.heat.data)
          cellphonedb.heat.data <- cellphonedb.heat.data2
          
                    #Add floor and ceiling:
                            #cellphonedb.heat.data <- as.matrix(cellphonedb.heat.data)
                            #cellphonedb.heat.data[which(cellphonedb.heat.data < -1)] <- -1
                            #cellphonedb.heat.data[which(cellphonedb.heat.data > 1)] <- 1
                            #cellphonedb.heat.data <- as.data.frame(cellphonedb.heat.data)
                            
          library(pheatmap)                          
          
          de.cellphonedb.heatmap <- pheatmap(cellphonedb.heat.data,
                                    cutree_cols = 5,
                                    cutree_rows = 3,
                                    color = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                    show_rownames = F,
                                    treeheight_row = 0,
                                    treeheight_col = 0,
                                    cluster_cols=FALSE)
        
          setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/rl_network")
          
          #ggsave(de.cellphonedb.heatmap, width = 80, height = 80,units = "cm",filename = "de.cellphonedb.heatmap.pdf", limitsize = FALSE)  
          
        
                
#### Old code:
#################################################        #################################################        #################################################
#################################################        #################################################        #################################################
#################################################        #################################################        #################################################
        
                  #Create node df:
                        lesion.list <- unique(network$Lesion)                                                 #Edit for sample:
                        
                        node_list <- c()
                        
                        for (i in 1:length(lesion.list)) {
                          
                          #Filter Edges:
                          lesion.edges <- filter(network, Lesion == lesion.list[i])                     #Edit for sample:
                          
                          #Create Nodes DFs:
                          
                          nodes <- unique(c(lesion.edges$from, lesion.edges$to))
                          
                          nodes.df <- data.frame(node_name = nodes,
                                                 CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                                 Lesion = lesion.list[i])
                          
                          #Assign Outputs:
                          assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
                          assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
                          
                          node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
                          
                        }
                        
                        
                        #Merge Nodes:
                        nodes.df <- purrr::reduce(mget(node_list), full_join)
                        nodes.df$Lesion <- factor(nodes.df$Lesion, 
                                                  levels = c("Nonlesional", "Irritant", "Acetone_Vehicle","Day2_Allergy", "Day4_Allergy"))
                        
                        #Create Edges.df:
                        edges.df <- network                                                                          
                        edges.df$Lesion <- factor(edges.df$Lesion, 
                                                  levels = c("Nonlesional", "Irritant", "Acetone_Vehicle" ,"Day2_Allergy", "Day4_Allergy"))
                        
                        
                        edges.df <- edges.df[,which(colnames(edges.df) %in% c("from", "to", "Lesion", "Ligand_Source", "Ligand", "Receptor_source", "Receptor","log10_Connection_product", "database_source"))]
                        
                        #Remove duplicate nichenet edges: (not sure why this happened. . . need to fix this eventually)
                            edges.df <- edges.df %>% mutate(edgeID = paste(from, "__", to, "__", Lesion, sep = ""))
                            
                            duplicates <-  edges.df$edgeID[which(duplicated(edges.df$edgeID))]         
                            duplicates <- edges.df[which(edges.df$edgeID %in% duplicates),]
                            duplicate.max.connections <- duplicates %>% group_by(edgeID) %>% summarize(max.connection = max(log10_Connection_product))
                            
                            duplicates <- duplicates[-which(duplicated(duplicates$edgeID)),]
                            duplicates <- full_join(duplicates, duplicate.max.connections)
                            duplicates <- duplicates[,-which(colnames(duplicates) == "log10_Connection_product")]
                            
                            duplicates <- duplicates %>% rename(log10_Connection_product = "max.connection")
                            
                            edges.df <- edges.df[-which(edges.df$edgeID %in% duplicates$edgeID),]
                            edges.df <- bind_rows(edges.df, duplicates)
                            
                            
                        # Scale edge with based on fold change strength relative to nonlesional connection:
                            
                            ### Instead of fold change over Nonlesional, should I try scaling by percent of max connection strength?
                            ### Why is Day4 allergy Myeloid:CXCL9-KRT:CXCR3 almost as strong as Day4 Allergy lymphocyte CXCR3 connections?
                            
                            
                            max.expression.vals <- edges.df %>% group_by(Ligand, Receptor) %>% summarize(max.expression = max(log10_Connection_product))
                            
                            scaled_edges <- edges.df
                            
                            for (r in 1:nrow(max.expression.vals)) {
                              
                              connection.row <- max.expression.vals[r,]
                              connection.max <- max.expression.vals[r,"max.expression"]
                                
                              connection.edges <- edges.df %>% filter(Ligand == connection.row$Ligand & Receptor == connection.row$Receptor)
                              connection.edges$log10_Connection_product <- sapply(connection.edges$log10_Connection_product, FUN = function(x) as.numeric(round(10^(x)/(10^(connection.max)), digits = 3)))
                              
                              scaled_edges[which(scaled_edges$edgeID %in% connection.edges$edgeID),] <- connection.edges
                            }
                            
                            edges.df <- scaled_edges
                            colnames(edges.df)[which(colnames(edges.df) == "log10_Connection_product")] <- "Scaled_Connection_Product"
                              
                            
                      #Add edge DE status:
                            edges.df$DE <- FALSE
                        
                            for (l in 1:length(unique(edges.df$Lesion))) {
                              les <- unique(edges.df$Lesion)[l]
                              
                              lesion.edges <- edges.df %>% filter(Lesion == les & Scaled_Connection_Product > 0)
                              
                              les.celltypes <- unique(lesion.edges$Ligand_Source)
                              
                                      for (lc in 1:length(les.celltypes)) {
                                
                                les.cell <- les.celltypes[lc]
                                
                                les.cell.edges <- lesion.edges %>% filter(Ligand_Source == les.cell)
                                
                                les.cell.DE <- DE.table %>% filter(Comparison == les & CellType == les.cell & Color_factor == "Significantly Up") %>% dplyr::pull(Gene)
                                
                                DE.complex.genes <- deconvoluted %>% filter(gene_name %in% les.cell.DE) %>% dplyr::pull(complex_name) %>% unique()
                                
                                les.cell.edges.DE <- les.cell.edges %>% filter(Ligand %in% unique(les.cell.DE, DE.complex.genes)) %>% dplyr::pull(edgeID)
                                
                                if (length(les.cell.edges.DE) >0) {
                                      edges.df[which(edges.df$edgeID %in% les.cell.edges.DE), "DE"] <- TRUE
                                }
                                
                          }
                            }
                                
                            
                        #Prune edges to only interconnected CellphoneDB & Nichenet edges, plus all DE CellphoneDB edges:
                          filtered.edges.df <- as.data.frame(edges.df)[0,]
                          
                          for (l in 1:length(unique(edges.df$Lesion))) {
                            
                            les <- unique(edges.df$Lesion)[l]
                          
                            lesion.edges <- edges.df %>% filter(Lesion == les & Scaled_Connection_Product > 0)
                            
                                #Prune nichenet edges:
                                lesion.nichenet.edges <- lesion.edges %>% filter(database_source == "Nichenet" & Scaled_Connection_Product > 0) 
                                
                                nichenet.edges.to.prune <- c()
                                
                                for (ne in 1:length(unique(lesion.nichenet.edges$from))) {
                                  
                                  nichenet.edge.from <- unique(lesion.nichenet.edges$from)[ne]
                                  
                                  upstream.edges <- lesion.edges %>% filter(to == nichenet.edge.from & Scaled_Connection_Product > 0 & database_source == "CellphoneDB")
                                  
                                  if (nrow(upstream.edges) == 0) {
                                    nichenet.edges.to.prune <- c(nichenet.edges.to.prune, nichenet.edge.from)
                                  }
                                  
                                }    
                                
                                lesion.nichenet.edges <- lesion.nichenet.edges %>% filter(!from %in% nichenet.edges.to.prune)
                                
                                #Prune CellphoneDB edges:
                                lesion.cellphonedb.edges <- lesion.edges %>% filter(database_source == "CellphoneDB" & Scaled_Connection_Product > 0) 
                                
                                cellphonedb.edges.to.prune <- c()
                                
                                for (ce in 1:length(unique(lesion.cellphonedb.edges$to))) {
                                  
                                  cellphonedb.edge.to <- unique(lesion.cellphonedb.edges$to)[ce]
                                  
                                  downstream.edges <- lesion.edges %>% filter(from == cellphonedb.edge.to & Scaled_Connection_Product > 0 & database_source == "Nichenet")
                                  
                                  if (nrow(downstream.edges) == 0) {
                                    cellphonedb.edges.to.prune <- c(cellphonedb.edges.to.prune, cellphonedb.edge.to)
                                  }
                                  
                                }
                                
                                    #Add back any nodes related to DE CellphoneDB edges:
                                    #cellphonedb.edges.to.not.prune <- c()
                                        
                                   # for (cep in 1:length(cellphonedb.edges.to.prune)) {
                                     # node.to <- cellphonedb.edges.to.prune[cep]
                                      
                                      #node.edges.in.all.lesions <- edges.df %>% filter(to == node.to)
                                      
                                          #if (any(node.edges.in.all.lesions$DE == TRUE)) {
                                            #cellphonedb.edges.to.not.prune <- c(cellphonedb.edges.to.not.prune, node.to)
                                          #}
                                    #}
                                
                                    #if (length(cellphonedb.edges.to.not.prune)>0) {
                                      #cellphonedb.edges.to.prune <- cellphonedb.edges.to.prune[-which(cellphonedb.edges.to.prune %in% cellphonedb.edges.to.not.prune)]
                                    #}
                                    
                                lesion.cellphonedb.edges <- lesion.cellphonedb.edges %>% filter(!to %in% cellphonedb.edges.to.prune)
                                
                            lesion.edges <- bind_rows(lesion.cellphonedb.edges, lesion.nichenet.edges)         
                          
                            filtered.edges.df <- bind_rows(filtered.edges.df, lesion.edges)
                          }
                            
                          filtered.edges.df$Lesion <- factor(filtered.edges.df$Lesion, 
                                                    levels = c("Nonlesional", "Irritant", "Acetone_Vehicle" ,"Day2_Allergy", "Day4_Allergy"))
                          
                        #Duplicate node names columns so that each is preserved after graph formation:
                          filtered.edges.df$from.named <- filtered.edges.df$from
                          filtered.edges.df$to.named <- filtered.edges.df$to
                          
                          
                        #Create filtered.nodes.df:
                          node_list <- c()
                          
                          for (i in 1:length(unique(filtered.edges.df$Lesion))) {
                            
                            les <- as.character(unique(filtered.edges.df$Lesion))[i]
                            #Filter Edges:
                            lesion.edges <- filter(filtered.edges.df, Lesion == les)                     #Edit for sample:
                            
                            #Create Nodes DFs:
                            
                            nodes <- unique(c(lesion.edges$from, lesion.edges$to))
                            
                            nodes.df <- data.frame(node_name = nodes,
                                                   CellType = str_replace_all(string = nodes, pattern = ":.*$", replacement = ""),
                                                   Lesion = lesion.list[i])
                            
                            #Assign Outputs:
                            assign(paste(lesion.list[i], "_nodes.df", sep = ""), nodes.df)
                            assign(paste(lesion.list[i], "_edges.df", sep = ""), lesion.edges)
                            
                            node_list <- c(node_list, paste(lesion.list[i], "_nodes.df", sep = ""))
                            
                          }
                          
                          
                          #Merge Nodes:
                          filtered.nodes.df <- purrr::reduce(mget(node_list), full_join)
                          filtered.nodes.df$Lesion <- factor(filtered.nodes.df$Lesion, 
                                                    levels = c("Nonlesional", "Irritant", "Acetone_Vehicle","Day2_Allergy", "Day4_Allergy"))
                          
                          
                        #Make Graph:    
                        RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
                        
                        #Correct misconnected edges:
                        graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
                        graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
                        
                        duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
                        
                        filtered.nodes.df <- rownames_to_column(filtered.nodes.df)
                        
                        for (r in 1:nrow(duplicates)) {
                          
                          les <- as.character(duplicates[r,"Lesion"])
                          
                          if (duplicates[r,"database_source"] == "CellphoneDB") {
                            from <- duplicates[r,"from.named"]
                            
                            to <- duplicates[r,"to.named"]
                            
                            duplicates[r,"from"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                            duplicates[r,"to"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                            
                          }
                          
                          if (duplicates[r,"database_source"] == "Nichenet") {
                            from <- duplicates[r,"from.named"]
                            
                            to <- duplicates[r,"to.named"]
                            
                            duplicates[r,"from"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                            duplicates[r,"to"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                            
                          }
                        }
                        
                        graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
                        
                        RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
                        
                        
                        
                        # Define Layout using merged data from all lesions:
                        nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
                        nodes <- nodes[-which(duplicated(nodes$node_name)),]  
                        
                        merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
                        
                        layout <- layout_with_stress(merged.graph)
                        colnames(layout) <- c("x", "y")
                        
                        layout <- cbind(nodes, layout)
                        
                        #save(RL_graph, file = "Data/RL_network2/savedfiles/RL_graph.Rdata")
                        #save(layout_full, file = "Data/RL_network2/savedfiles/layout_full.Rdata")
                        #save(layout, file = "Data/RL_network2/savedfiles/layout_merged.Rdata")
                        #save(filtered.nodes.df, file = "Data/RL_network2/savedfiles/filtered_nodes.df.Rdata")
                        #save(filtered.edges.df, file = "Data/RL_network2/savedfiles/filtered_edges.df.Rdata")
                        
                        #load("Data/RL_network2/savedfiles/RL_graph.Rdata")
                        #load("Data/RL_network2/savedfiles/layout_full.Rdata")
                        #load("Data/RL_network2/savedfiles/layout_merged.Rdata")
                        #load("Data/RL_network2/savedfiles/filtered_nodes.df.Rdata")
                        #load("Data/RL_network2/savedfiles/filtered_edges.df.Rdata")
                        
                            ## Ligand Edge Label mapping:
                            label.mapping.df <- filtered.edges.df
                            label.mapping.df$label_x <- NA
                            label.mapping.df$label_y <- NA
                            label.mapping.df$Edge_label <- ""
                            
                            for (l in 1:nrow(label.mapping.df)) {
                            
                                lig.from <- label.mapping.df$from.named[l]
                                lig.to <- label.mapping.df$to.named[l]
                                
                                node1_info <- layout %>% filter(node_name == lig.from)
                                node2_info <- layout %>% filter(node_name == lig.to)
                                
                                label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
                                label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
                            }
                            
                            #Only label ligands that are reasonably far away from each other:
                            ## Note to self: should i prioritize DE ligand edges?
                            max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
                            
                            for (l in 1:length(unique(label.mapping.df$Ligand))) {
                              lig <- unique(label.mapping.df$Ligand)[l]
                              
                              lig.edges <- label.mapping.df %>% filter(Ligand == lig)
                              
                              lig.edges$Edge_label <- "temp"
                              
                                  while(any(lig.edges$Edge_label =="temp")){
                                    
                                    #temporarily make a column for fraction of max distance from first ungrouped edge:
                                    grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                                    
                                    grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                                    
                                    grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.035) 
                                    
                                    #Label closest to centroid edges grouped by 5% of max layout distance:
                                    centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                                    centroid_x <- sum(centroid_x)/length(centroid_x)
                                    centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                                    centroid_y <- sum(centroid_y)/length(centroid_y)
                                    
                                    grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                                    
                                    group.edges <- grouped.lig.edges %>% arrange(centroid_dist) %>% dplyr::pull(edgeID)
                                    Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)

                                    lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                                    
                                    lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                                    
                                  }
                              
                              label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
                              
                            }
                            
                            
                        #Copy coordinates to node df containing all faceted nodes:
                        layout_full <- full_join(filtered.nodes.df, layout)
                        
                            #Set point colors for plotting:
                            GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
                            
                            ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
                            node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
                            
                        ##Plot With Facets:
                        ggraph(RL_graph, layout = "manual",
                               x = layout_full[, "x"],
                               y = layout_full[, "y"]) + 
                          geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
                          geom_node_point(aes(color = CellType), size = 0.5) + 
                          geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                          scale_edge_width(range = c(0.1, 0.5)) +
                          theme_graph() + 
                          theme(legend.position = "bottom") +
                          facet_nodes(~Lesion) +
                          scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
                          new_scale_color()+
                          geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0) +
                          scale_color_manual(values = c("#000000", "#750000"))
                          
                        ggraph(merged.graph, layout = "manual",
                               x = layout[, "x"],
                               y = layout[, "y"]) + 
                          geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
                          geom_node_point(aes(color = CellType), size = 0.5) + 
                          geom_edge_link(alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                          scale_edge_width(range = c(0.1, 0.5)) +
                          theme_graph() + 
                          theme(legend.position = "bottom") +
                          geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label), size = 2.5, alpha =1, box.padding = 0)
                          
                        
                        
                        
                        ## Find Cycles in network:
                                  FindCycles = function(g) {
                                    Cycles = NULL
                                    for(v1 in V(g)) {
                                      if(degree(g, v1, mode="in") == 0) { next }
                                      GoodNeighbors = neighbors(g, v1, mode="out")
                                      GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
                                      for(v2 in GoodNeighbors) {
                                        TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
                                        TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
                                        TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
                                        Cycles  = c(Cycles, TempCyc)
                                      }
                                    }
                                    Cycles
                        }
                                  merged.igraph <- as.igraph(merged.graph)
                                  network.cycles <- FindCycles(merged.igraph)
                                  
                                  nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
                                  
                                  for (c in 1:length(network.cycles)) {
                                    cycle <- network.cycles[[c]]
                                    
                                    named.cycle <- c()  
                                    
                                    for (n in 1:length(cycle)) {
                                      node.num <- cycle[n]
                                      node.info <- nodes[as.numeric(node.num),]
                                      named.node <- node.info$node_name
                                      
                                      named.cycle <- c(named.cycle, named.node)
                                    }
                                    
                                    if (c == 1) {
                                      named.network.cycles <-list(named.cycle)
                                    }
                                    if (c > 1) {
                                      named.network.cycles <- append(named.network.cycles, list(named.cycle))
                                    }
                        }
                                  
                                  #Edges causing cycles:
                                  
                                  for (c in 1:length(named.network.cycles)) {
                                    named.cycle <- named.network.cycles[[c]]
                                    
                                    first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])
                                    last.edges <- filtered.edges.df %>% filter(from == named.cycle[length(named.cycle)-1] & to == named.cycle[length(named.cycle)])
                                        
                                          if (c ==1) {
                                            cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                                              first.edge.from = unique(first.edges$from.named),
                                                                              first.edge.to = unique(first.edges$to.named),
                                                                              first.edges.num= nrow(first.edges),
                                                                              first.edges.database = unique(first.edges$database_source),
                                                                              last.edge.from = unique(last.edges$from.named),
                                                                              last.edge.to = unique(last.edges$to.named),
                                                                              last.edges.num = nrow(last.edges),
                                                                              last.edges.database = unique(last.edges$database_source))
                                          }
                                          
                                          if (c > 1) {
                                            cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                                           first.edge.from = unique(first.edges$from.named),
                                                                           first.edge.to = unique(first.edges$to.named),
                                                                           first.edges.num= nrow(first.edges),
                                                                           first.edges.database = unique(first.edges$database_source),
                                                                           last.edge.from = unique(last.edges$from.named),
                                                                           last.edge.to = unique(last.edges$to.named),
                                                                           last.edges.num = nrow(last.edges),
                                                                           last.edges.database = unique(last.edges$database_source))
                                            cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                                          }
                                    
                                    }
                                  cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(first.edge.from, "__", first.edge.to, sep = ""),
                                                                                      last.edge = paste(last.edge.from, "__", last.edge.to, sep = ""))
                                  cyclic.edge.info.df <- mutate(cyclic.edge.info.df, cycle_id = paste(first.edge, ".....",last.edge, sep = ""))
                                  
                                  #All cycles look to be due to weak nichenet connections at last edge:
                                  if (unique(cyclic.edge.info.df$last.edges.database) == "Nichenet") {
                                    cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$last.edge)),]
                                    cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("last.edge.from", "last.edge.to"))]
                                  }
                                  
                                  for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                                    
                                    cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$last.edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$last.edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                                  
                                    if (length(cyclic.edgeID) >0) {
                                      filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                                    }
                                    rm(cyclic.edgeID)
                                  }
                                  
                            #Plot network without cyclic edges:
                          
                                  #Make Graph:    
                                  RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
                                  
                                  #Correct misconnected edges:
                                  graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
                                  graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
                                  
                                  duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
                                  
                                  for (r in 1:nrow(duplicates)) {
                                    
                                    les <- as.character(duplicates[r,"Lesion"])
                                    
                                    if (duplicates[r,"database_source"] == "CellphoneDB") {
                                      from <- duplicates[r,"from.named"]
                                      
                                      to <- duplicates[r,"to.named"]
                                      
                                      duplicates[r,"from"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                      duplicates[r,"to"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                      
                                    }
                                    
                                    if (duplicates[r,"database_source"] == "Nichenet") {
                                      from <- duplicates[r,"from.named"]
                                      
                                      to <- duplicates[r,"to.named"]
                                      
                                      duplicates[r,"from"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                      duplicates[r,"to"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                      
                                    }
                                  }
                                  
                                  graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
                                  
                                  RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
                                  
                                  
                                  
                                  # Define Layout using merged data from all lesions:
                                  nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType"))]
                                  nodes <- nodes[-which(duplicated(nodes$node_name)),]  
                                  
                                  merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
                                  
                                  layout <- layout_with_stress(merged.graph)
                                  colnames(layout) <- c("x", "y")
                                  
                                  layout <- cbind(nodes, layout)
                                
                                  #Copy coordinates to node df containing all faceted nodes:
                                  layout_full <- full_join(filtered.nodes.df, layout)
                                  
                                  #Use igraph to find feedback edges:
                                  merged.igraph <- as.igraph(merged.graph)
                                  cycles <- feedback_arc_set(merged.igraph)
                                  cycle.edge.ids <- as_ids(cycles)
                                  

                                      #convert to list:
                                        for (e in 1:length(cycle.edge.ids)) {
                                          if (e ==1){
                                            cyclic.edges <- list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]])))
      
                                          }
                                          
                                          if (e>1){
                                            cyclic.edges <-  append(cyclic.edges, list(as.numeric(ends(merged.igraph, E(merged.igraph)[cycle.edge.ids[e]]))))
                                          }
                                        }   
                                  
                                      #Add names:
                                      for (c in 1:length(cyclic.edges)) {
                                        cycle <- cyclic.edges[[c]]
                                        
                                        named.cycle <- c()  
                                        
                                        for (n in 1:length(cycle)) {
                                          node.num <- cycle[n]
                                          node.info <- nodes[as.numeric(node.num),]
                                          named.node <- node.info$node_name
                                          
                                          named.cycle <- c(named.cycle, named.node)
                                        }
                                        
                                        if (c == 1) {
                                          named.network.cycles <-list(named.cycle)
                                        }
                                        if (c > 1) {
                                          named.network.cycles <- append(named.network.cycles, list(named.cycle))
                                        }
                                      }
                                  
                                  for (c in 1:length(named.network.cycles)) {
                                    named.cycle <- named.network.cycles[[c]]
                                    
                                    first.edges <- filtered.edges.df %>% filter(from == named.cycle[1] & to == named.cycle[2])

                                    if (c ==1) {
                                      cyclic.edge.info.df <- data.frame(cycle.num= as.numeric(c),
                                                                        edge.from = unique(first.edges$from.named),
                                                                        edge.to = unique(first.edges$to.named),
                                                                        edges.num= nrow(first.edges),
                                                                        edges.database = unique(first.edges$database_source))
                                                                        
                                    }
                                    
                                    if (c > 1) {
                                      cyclic.edge.info <- data.frame(cycle.num= as.numeric(c),
                                                                     edge.from = unique(first.edges$from.named),
                                                                     edge.to = unique(first.edges$to.named),
                                                                     edges.num= nrow(first.edges),
                                                                     edges.database = unique(first.edges$database_source))
                                                                    
                                      cyclic.edge.info.df <- bind_rows(cyclic.edge.info.df, cyclic.edge.info)
                                    }
                                    
                                  }
                                      
                                  cyclic.edge.info.df <- mutate(cyclic.edge.info.df, first.edge = paste(edge.from, "__", edge.to, sep = ""))

                                  #All cycles look to be due to weak nichenet connections at last edge:
                                  if (unique(cyclic.edge.info.df$edges.database) == "Nichenet") {
                                    cyclic_edges_to_prune = cyclic.edge.info.df[-which(duplicated(cyclic.edge.info.df$first.edge)),]
                                    cyclic_edges_to_prune = cyclic_edges_to_prune[,which(colnames(cyclic_edges_to_prune) %in% c("edge.from", "edge.to"))]
                                  }
                                  
                                  for (cycle.edge in 1:nrow(cyclic_edges_to_prune)) {
                                    
                                    cyclic.edgeID <-filtered.edges.df %>% filter(from.named == cyclic_edges_to_prune$edge.from[cycle.edge] & to.named == cyclic_edges_to_prune$edge.to[cycle.edge]) %>% dplyr::pull(edgeID)
                                    
                                    if (length(cyclic.edgeID) >0) {
                                      filtered.edges.df <- filtered.edges.df[-which(filtered.edges.df$edgeID %in% cyclic.edgeID),]
                                    }
                                    rm(cyclic.edgeID)
                                  }
                               
                                  #Plot network without cyclic edges:
                                  
                                  #Make Graph:    
                                      #add shapes to nodes based on basic celltype:
                                       filtered.nodes.df <- filtered.nodes.df %>% mutate(Basic_Celltype = case_when(CellType == "KRT" ~ "KRT",
                                                                                          CellType == "MEL" ~ "MEL",
                                                                                          CellType %in% c("Myeloid", "pDC", "LC") ~ "APC",
                                                                                          CellType %in% c("CD4", "CD8", "NK", "Treg") ~ "LYMPH"))
                                       filtered.nodes.df$Basic_Celltype <- factor(filtered.nodes.df$Basic_Celltype, levels = c("KRT", "MEL", "APC", "LYMPH"))
                                       filtered.nodes.df$CellType <- factor(filtered.nodes.df$CellType, levels = c("KRT", "MEL", "LC", "Myeloid", "pDC", "CD4", "Treg", "CD8", "NK"))

                                       
                                  RL_graph <- tbl_graph(nodes = filtered.nodes.df, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
                                  
                                  #Correct misconnected edges:
                                  graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
                                  graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
                                  
                                  duplicates <- graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),]  
                                  
                                  for (r in 1:nrow(duplicates)) {
                                    
                                    les <- as.character(duplicates[r,"Lesion"])
                                    
                                    if (duplicates[r,"database_source"] == "CellphoneDB") {
                                      from <- duplicates[r,"from.named"]
                                      
                                      to <- duplicates[r,"to.named"]
                                      
                                      duplicates[r,"from"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                      duplicates[r,"to"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                      
                                    }
                                    
                                    if (duplicates[r,"database_source"] == "Nichenet") {
                                      from <- duplicates[r,"from.named"]
                                      
                                      to <- duplicates[r,"to.named"]
                                      
                                      duplicates[r,"from"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                      duplicates[r,"to"] <- filtered.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                      
                                    }
                                  }
                                  
                                  graph.edges[which(duplicated(graph.edges$to)|duplicated(graph.edges$from)),] <- duplicates  
                                  
                                  RL_graph <- RL_graph %>% activate(edges) %>% right_join(graph.edges)
                                  
                                  # Define Layout using merged data from all lesions:
                                  nodes <- filtered.nodes.df[,which(colnames(filtered.nodes.df) %in% c("node_name", "CellType", "Basic_Celltype", "node_label"))]
                                  nodes <- nodes[-which(duplicated(nodes$node_name)),]  
                                  
                                  merged.graph <- tbl_graph(nodes = nodes, edges = filtered.edges.df, directed = TRUE, node_key = "node_name")
                                  
                                  layout <- layout_with_stress(merged.graph)
                                  colnames(layout) <- c("x", "y")
                                  
                                  layout <- cbind(nodes, layout)
                                  
                                  ## Ligand Edge Label mapping:
                                          label.mapping.df <- filtered.edges.df
                                          label.mapping.df$label_x <- NA
                                          label.mapping.df$label_y <- NA
                                          label.mapping.df$Edge_label <- ""
                                          
                                          for (l in 1:nrow(label.mapping.df)) {
                                            
                                            lig.from <- label.mapping.df$from.named[l]
                                            lig.to <- label.mapping.df$to.named[l]
                                            
                                            node1_info <- layout %>% filter(node_name == lig.from)
                                            node2_info <- layout %>% filter(node_name == lig.to)
                                            
                                            label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
                                            label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
                                          }
                                          
                                          #Only label ligands that are reasonably far away from each other:
                                          ## Note to self: should i prioritize DE ligand edges?
                                          max.layout.dist <- ((max(layout$x)-min(layout$x))^2+(max(layout$y)-min(layout$y))^2)^(0.5)
                                          
                                          for (l in 1:length(unique(label.mapping.df$Ligand))) {
                                            lig <- unique(label.mapping.df$Ligand)[l]
                                            
                                            lig.edges <- label.mapping.df %>% filter(Ligand == lig)
                                            
                                            lig.edges$Edge_label <- "temp"
                                            
                                            while(any(lig.edges$Edge_label =="temp")){
                                              
                                              #temporarily make a column for fraction of max distance from first ungrouped edge:
                                              grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                                              
                                              grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                                              
                                              grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.05) 
                                              
                                              #Label closest to centroid edges grouped by 5% of max layout distance:
                                              centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                                              centroid_x <- sum(centroid_x)/length(centroid_x)
                                              centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                                              centroid_y <- sum(centroid_y)/length(centroid_y)
                                              
                                              grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                                              
                                              group.edges <- grouped.lig.edges %>% arrange(centroid_dist) %>% dplyr::pull(edgeID)
                                              Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
                                              
                                              lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                                              
                                              lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                                              
                                            }
                                            
                                            label.mapping.df[which(label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
                                            
                                          }
                                  
                                  ## Node Labels:
                                          node.labels.df <- layout
                                          node.labels.df$node_label <- ""
                                          
                                          node.labels.df$node_gene <- sapply(node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
                                          
                                          #Only label nodes that are reasonably far away from each other:
                                          
                                          for (g in 1:length(unique(node.labels.df$node_gene))) {
                                            gene <- unique(node.labels.df$node_gene)[g]
                                            
                                            gene.nodes <- node.labels.df %>% filter(node_gene == gene)
                                            
                                            gene.nodes$node_label <- "temp"
                                            
                                            while(any(gene.nodes$node_label =="temp")){
                                              
                                              #temporarily make a column for fraction of max distance from first ungrouped edge:
                                              grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                                              
                                              grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.layout.dist), digits = 3))
                                              
                                              grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.035) 
                                              
                                              #Label closest to centroid edges grouped by 5% of max layout distance:
                                              centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                                              centroid_x <- sum(centroid_x)/length(centroid_x)
                                              centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                                              centroid_y <- sum(centroid_y)/length(centroid_y)
                                              
                                              grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                                              
                                              group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                                              nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                                              
                                              gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                                              
                                              gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                                              
                                            }
                                            
                                            node.labels.df[which(node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
                                            
                                          }
                              
                                  #Copy coordinates to node df containing all faceted nodes:
                                  layout_full <- full_join(filtered.nodes.df, node.labels.df)
                                  
                                  #Set point colors for plotting:
                                  GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
                                  
                                  ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
                                  node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)

                                  ##Plot With Facets:
                                  ggraph(RL_graph, layout = "manual",
                                         x = layout_full[, "x"],
                                         y = layout_full[, "y"]) + 
                                    geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
                                    geom_node_point(aes(color = CellType), size = 0.5) + 
                                    geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                                    scale_edge_width(range = c(0.1, 0.5)) +
                                    theme_graph() + 
                                    theme(legend.position = "bottom") +
                                    facet_nodes(~Lesion) +
                                    scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
                                    new_scale_color()+
                                    scale_color_manual(values = c("#000000", "#750000")) +
                                    geom_text_repel(data = label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0)

                                  ggraph(merged.graph, layout = "manual",
                                         x = layout[, "x"],
                                         y = layout[, "y"]) + 
                                    geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
                                    geom_node_point(aes(color = CellType), size = 2) + 
                                    geom_edge_link(alpha = 0.7,width = 0.05, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                                    theme_graph() + 
                                    theme(legend.position = "bottom")
                                  
                                  
                                  
                                  #Sugiyama Layout:
                                  sugi_layout <- layout_with_sugiyama(merged.graph, hgap = 500)
                                  sugi_layout <- as.data.frame(sugi_layout$layout)
                                  colnames(sugi_layout) <- c("x", "y")
                                 
                                  sugi_layout <- cbind(nodes, sugi_layout)
                                  
                                  ## Sugi Node Labels:
                                      sugi.node.labels.df <- sugi_layout
                                      sugi.node.labels.df$node_label <- ""
                                      
                                      sugi.node.labels.df$node_gene <- sapply(sugi.node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
                                      
                                      max.sugi.layout.dist <- ((max(sugi_layout$x)-min(sugi_layout$x))^2+(max(sugi_layout$y)-min(sugi_layout$y))^2)^(0.5)
                                      
                                      #Only label nodes that are reasonably far away from each other:
                                      for (g in 1:length(unique(sugi.node.labels.df$node_gene))) {
                                        gene <- unique(sugi.node.labels.df$node_gene)[g]
                                        
                                        gene.nodes <- sugi.node.labels.df %>% filter(node_gene == gene)
                                        
                                        gene.nodes$node_label <- "temp"
                                        
                                        while(any(gene.nodes$node_label =="temp")){
                                          
                                          #temporarily make a column for fraction of max distance from first ungrouped edge:
                                          grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                                          
                                          grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/max.sugi.layout.dist), digits = 3))
                                          
                                          grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.05) 
                                          
                                          #Label closest to centroid edges grouped by 5% of max layout distance:
                                          centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                                          centroid_x <- sum(centroid_x)/length(centroid_x)
                                          centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                                          centroid_y <- sum(centroid_y)/length(centroid_y)
                                          
                                          grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                                          
                                          group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                                          nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                                          
                                          gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                                          
                                          gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                                          
                                        }
                                        
                                        sugi.node.labels.df[which(sugi.node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
                                        
                                      }
                                  
                                  sugi_layout_full <- full_join(filtered.nodes.df, sugi.node.labels.df)
                                      
                                  
                                  #Build Custom Color Palette:
                                  flat.palette <- c("#F44336", "#C2185B", "#BA68C8", "#673AB7", "#7986CB", "#1E88E5", "#26C6DA", "#00897B", "#4CAF50",  "#AFB42B", "#FF8F00", "#FF5722", "#8D6E63", "#78909C")
                                  GetPalette = colorRampPalette(flat.palette)

                                  ColorCount = as.numeric(length(unique(nodes$CellType)))
                                  custom.colors <- GetPalette(ColorCount)
                                  
                                  
                                  ggraph(merged.graph, layout = "manual",
                                         x = sugi_layout[, "x"],
                                         y = sugi_layout[, "y"]) + 
                                    geom_node_point(data = sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
                                    geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
                                    scale_shape_manual(values = c(15, 8, 13, 16))+
                                    scale_color_manual(values = custom.colors) +
                                    geom_edge_diagonal(aes(color = database_source), alpha = 0.35, width = 0.1, label_push = TRUE, strength = 1, arrow = arrow(angle =20, length = unit(1.5, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
                                    new_scale_color() +
                                    scale_edge_colour_manual(values = c("#2E59B7", "#B72EB7")) +
                                    theme_graph() + 
                                    theme(legend.position = "bottom")+
                                    expand_limits(y = c(min(sugi.node.labels.df$y), max(sugi.node.labels.df$y)+0.15*(range(sugi.node.labels.df$y)[2]-range(sugi.node.labels.df$y)[1]))) +
                                    geom_text_repel(data = sugi.node.labels.df, aes(x = x, y = y, label = node_label), nudge_y = 0.15, size = 2.5, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5)
                                  
                                  
                                  ggraph(RL_graph, layout = "manual",
                                         x = sugi_layout_full[, "x"],
                                         y = sugi_layout_full[, "y"]) + 
                                    geom_node_point(data = sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
                                    geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
                                    scale_shape_manual(values = c(15, 8, 13,16))+
                                    geom_edge_diagonal(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
                                    scale_edge_width(range = c(0.1, 0.5)) +
                                    theme_graph() + 
                                    theme(legend.position = "bottom") +
                                    facet_nodes(~Lesion) +
                                    scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
                                    expand_limits(y = c(min(sugi.node.labels.df$y), max(sugi.node.labels.df$y)+0.05*(range(sugi.node.labels.df$y)[2]-range(sugi.node.labels.df$y)[1]))) +
                                    geom_text_repel(data = sugi.node.labels.df, aes(x = x, y = y, label = node_label), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5)
                                  
                                  
                          ### make into function for use with other data/ disease sets.        
                                  
                                  
                                  
                                  
                                       
                    #### Only Day 2 Irritant vs Day 2 allergic:
                                  
                                  d2compnetwork.nodes <- filtered.nodes.df %>% filter(Lesion %in% c("Irritant", "Day2_Allergy"))
                                  d2compnetwork.edges <- filtered.edges.df %>% filter(Lesion %in% c("Irritant", "Day2_Allergy"))
                                  
                                  Day2_graph <- tbl_graph(nodes = d2compnetwork.nodes, edges = d2compnetwork.edges, directed = TRUE, node_key = "node_name")
                                  
                                  #Correct misconnected edges:
                                  Day2_graph.nodes <- Day2_graph %>% activate(nodes) %>% as.data.frame()
                                  Day2_graph.edges <- Day2_graph %>% activate(edges) %>% as.data.frame()
                                  
                                  Day2_duplicates <- Day2_graph.edges[which(duplicated(Day2_graph.edges$to)|duplicated(Day2_graph.edges$from)),]  
                                  
                                  for (r in 1:nrow(Day2_duplicates)) {
                                    
                                    les <- as.character(Day2_duplicates[r,"Lesion"])
                                    
                                    if (Day2_duplicates[r,"database_source"] == "CellphoneDB") {
                                      from <- Day2_duplicates[r,"from.named"]
                                      
                                      to <- Day2_duplicates[r,"to.named"]
                                      
                                      Day2_duplicates[r,"from"] <- d2compnetwork.nodes %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                      Day2_duplicates[r,"to"] <- d2compnetwork.nodes %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                      
                                    }
                                    
                                    if (Day2_duplicates[r,"database_source"] == "Nichenet") {
                                      from <- Day2_duplicates[r,"from.named"]
                                      
                                      to <- Day2_duplicates[r,"to.named"]
                                      
                                      Day2_duplicates[r,"from"] <- d2compnetwork.nodes %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                      Day2_duplicates[r,"to"] <- d2compnetwork.nodes %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                      
                                    }
                                  }
                                  
                                  Day2_graph.edges[which(duplicated(Day2_graph.edges$to)|duplicated(Day2_graph.edges$from)),] <- Day2_duplicates  
                                  
                                  Day2_graph <- Day2_graph %>% activate(edges) %>% right_join(Day2_graph.edges)
                                  
                                  
                                  
                                  # Define Layout using merged data from all lesions:
                                  day2_nodes <- d2compnetwork.nodes[,which(colnames(d2compnetwork.nodes) %in% c("node_name", "CellType", "Basic_Celltype", "node_label"))]
                                  day2_nodes <- day2_nodes[-which(duplicated(day2_nodes$node_name)),]  
                                  
                                  day2_merged.graph <- tbl_graph(nodes = day2_nodes, edges = d2compnetwork.edges, directed = TRUE, node_key = "node_name")
                                  
                                  day2_layout <- layout_with_stress(day2_merged.graph)
                                  colnames(day2_layout) <- c("x", "y")
                                  
                                  day2_layout <- cbind(day2_nodes, day2_layout)
                                  
                                  ## Ligand Edge Label mapping:
                                  day2_label.mapping.df <- d2compnetwork.edges
                                  day2_label.mapping.df$label_x <- NA
                                  day2_label.mapping.df$label_y <- NA
                                  day2_label.mapping.df$Edge_label <- ""
                                  
                                  for (l in 1:nrow(day2_label.mapping.df)) {
                                    
                                    lig.from <- day2_label.mapping.df$from.named[l]
                                    lig.to <- day2_label.mapping.df$to.named[l]
                                    
                                    node1_info <- day2_layout %>% filter(node_name == lig.from)
                                    node2_info <- day2_layout %>% filter(node_name == lig.to)
                                    
                                    day2_label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
                                    day2_label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
                                  }
                                  
                                  #Only label ligands that are reasonably far away from each other:
                                  ## Note to self: should i prioritize DE ligand edges?
                                  day2_max.layout.dist <- ((max(day2_layout$x)-min(day2_layout$x))^2+(max(day2_layout$y)-min(day2_layout$y))^2)^(0.5)
                                  
                                  for (l in 1:length(unique(day2_label.mapping.df$Ligand))) {
                                    lig <- unique(day2_label.mapping.df$Ligand)[l]
                                    
                                    lig.edges <- day2_label.mapping.df %>% filter(Ligand == lig)
                                    
                                    lig.edges$Edge_label <- "temp"
                                    
                                    while(any(lig.edges$Edge_label =="temp")){
                                      
                                      #temporarily make a column for fraction of max distance from first ungrouped edge:
                                      grouped.lig.edges <- lig.edges %>% filter(Edge_label == "temp")
                                      
                                      grouped.lig.edges <- grouped.lig.edges %>% mutate(label_dist = round(((((label_x-grouped.lig.edges$label_x[1])^2+(label_y-grouped.lig.edges$label_y[1])^2)^(0.5))/day2_max.layout.dist), digits = 3))
                                      
                                      grouped.lig.edges <- grouped.lig.edges %>% filter(label_dist <= 0.05) 
                                      
                                      #Label closest to centroid edges grouped by 5% of max layout distance:
                                      centroid_x <- grouped.lig.edges %>% dplyr::pull(label_x)
                                      centroid_x <- sum(centroid_x)/length(centroid_x)
                                      centroid_y <- grouped.lig.edges %>% dplyr::pull(label_y)
                                      centroid_y <- sum(centroid_y)/length(centroid_y)
                                      
                                      grouped.lig.edges <- grouped.lig.edges %>% mutate(centroid_dist = round((((label_x-centroid_x)^2+(label_y-centroid_y)^2)^(0.5)), digits = 3))
                                      
                                      group.edges <- grouped.lig.edges %>% arrange(centroid_dist) %>% dplyr::pull(edgeID)
                                      Edges.to.label <- grouped.lig.edges %>% filter(centroid_dist == min(grouped.lig.edges$centroid_dist)) %>% dplyr::pull(edgeID)
                                      
                                      lig.edges[which(lig.edges$edgeID %in% group.edges),"Edge_label"] <- ""
                                      
                                      lig.edges[which(lig.edges$edgeID %in% Edges.to.label),"Edge_label"] <- lig
                                      
                                    }
                                    
                                    day2_label.mapping.df[which(day2_label.mapping.df$edgeID %in% lig.edges$edgeID),"Edge_label"] <- lig.edges$Edge_label
                                    
                                  }
                                  
                                  ## Node Labels:
                                  day2_node.labels.df <- day2_layout
                                  day2_node.labels.df$node_label <- ""
                                  
                                  day2_node.labels.df$node_gene <- sapply(day2_node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
                                  
                                  #Only label nodes that are reasonably far away from each other:
                                  
                                  for (g in 1:length(unique(day2_node.labels.df$node_gene))) {
                                    gene <- unique(day2_node.labels.df$node_gene)[g]
                                    
                                    gene.nodes <- day2_node.labels.df %>% filter(node_gene == gene)
                                    
                                    gene.nodes$node_label <- "temp"
                                    
                                    while(any(gene.nodes$node_label =="temp")){
                                      
                                      #temporarily make a column for fraction of max distance from first ungrouped edge:
                                      grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                                      
                                      grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/day2_max.layout.dist), digits = 3))
                                      
                                      grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.035) 
                                      
                                      #Label closest to centroid edges grouped by 5% of max layout distance:
                                      centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                                      centroid_x <- sum(centroid_x)/length(centroid_x)
                                      centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                                      centroid_y <- sum(centroid_y)/length(centroid_y)
                                      
                                      grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                                      
                                      group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                                      nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                                      
                                      gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                                      
                                      gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                                      
                                    }
                                    
                                    day2_node.labels.df[which(day2_node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
                                    
                                  }
                                  
                                  #Copy coordinates to node df containing all faceted nodes:
                                  day2_layout_full <- full_join(d2compnetwork.nodes, day2_node.labels.df)
                                  
                                  #Set point colors for plotting:
                                  GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
                                  
                                  ColorCount = as.numeric(length(unique(d2compnetwork.nodes$CellType)))
                                  node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
                                  
                                  ##Plot With Facets:
                                  ggraph(Day2_graph, layout = "manual",
                                         x = day2_layout_full[, "x"],
                                         y = day2_layout_full[, "y"]) + 
                                    geom_node_point(data = day2_layout, color = "gray", alpha = 0.2, size = 0.5) + 
                                    geom_node_point(aes(color = CellType), size = 0.5) + 
                                    geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                                    scale_edge_width(range = c(0.1, 0.5)) +
                                    theme_graph() + 
                                    theme(legend.position = "bottom") +
                                    facet_nodes(~Lesion) +
                                    scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
                                    new_scale_color()+
                                    scale_color_manual(values = c("#000000", "#750000")) +
                                    geom_text_repel(data = day2_label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0)
                                  
                                  ggraph(day2_merged.graph, layout = "manual",
                                         x = day2_layout[, "x"],
                                         y = day2_layout[, "y"]) + 
                                    geom_node_point(data = day2_layout, color = "gray", alpha = 0.2, size = 0.5) + 
                                    geom_node_point(aes(color = CellType), size = 2) + 
                                    geom_edge_link(alpha = 0.7,width = 0.05, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                                    theme_graph() + 
                                    theme(legend.position = "bottom")
                                  
                                  
                                  
                                  #Sugiyama Layout:
                                  day2_sugi_layout <- layout_with_sugiyama(day2_merged.graph, hgap = 500)
                                  day2_sugi_layout <- as.data.frame(day2_sugi_layout$layout)
                                  colnames(day2_sugi_layout) <- c("x", "y")
                                  
                                  day2_sugi_layout <- cbind(day2_nodes, day2_sugi_layout)
                                  
                                  ## Sugi Node Labels:
                                  day2_sugi.node.labels.df <- day2_sugi_layout
                                  day2_sugi.node.labels.df$node_label <- ""
                                  
                                  day2_sugi.node.labels.df$node_gene <- sapply(day2_sugi.node.labels.df$node_name, FUN  = function (x) as.character(unlist(str_split(x, pattern = ":"))[2]))
                                  
                                  day2_max.sugi.layout.dist <- ((max(day2_sugi_layout$x)-min(day2_sugi_layout$x))^2+(max(day2_sugi_layout$y)-min(day2_sugi_layout$y))^2)^(0.5)
                                  
                                  #Only label nodes that are reasonably far away from each other:
                                  for (g in 1:length(unique(day2_sugi.node.labels.df$node_gene))) {
                                    gene <- unique(day2_sugi.node.labels.df$node_gene)[g]
                                    
                                    gene.nodes <- day2_sugi.node.labels.df %>% filter(node_gene == gene)
                                    
                                    gene.nodes$node_label <- "temp"
                                    
                                    while(any(gene.nodes$node_label =="temp")){
                                      
                                      #temporarily make a column for fraction of max distance from first ungrouped edge:
                                      grouped.gene.nodes <- gene.nodes %>% filter(node_label == "temp")
                                      
                                      grouped.gene.nodes <- grouped.gene.nodes %>% mutate(label_dist = round(((((x-grouped.gene.nodes$x[1])^2+(y-grouped.gene.nodes$y[1])^2)^(0.5))/day2_max.sugi.layout.dist), digits = 3))
                                      
                                      grouped.gene.nodes <- grouped.gene.nodes %>% filter(label_dist <= 0.05) 
                                      
                                      #Label closest to centroid edges grouped by 5% of max layout distance:
                                      centroid_x <- grouped.gene.nodes %>% dplyr::pull(x)
                                      centroid_x <- sum(centroid_x)/length(centroid_x)
                                      centroid_y <- grouped.gene.nodes %>% dplyr::pull(y)
                                      centroid_y <- sum(centroid_y)/length(centroid_y)
                                      
                                      grouped.gene.nodes <- grouped.gene.nodes %>% mutate(centroid_dist = round((((x-centroid_x)^2+(y-centroid_y)^2)^(0.5)), digits = 3))
                                      
                                      group.nodes <- grouped.gene.nodes %>% arrange(centroid_dist) %>% dplyr::pull(node_name)
                                      nodes.to.label <- grouped.gene.nodes %>% filter(centroid_dist == min(grouped.gene.nodes$centroid_dist)) %>% dplyr::pull(node_name)
                                      
                                      gene.nodes[which(gene.nodes$node_name %in% group.nodes),"node_label"] <- ""
                                      
                                      gene.nodes[which(gene.nodes$node_name %in% nodes.to.label),"node_label"] <- gene
                                      
                                    }
                                    
                                    day2_sugi.node.labels.df[which(day2_sugi.node.labels.df$node_name %in% gene.nodes$node_name),"node_label"] <- gene.nodes$node_label
                                    
                                  }
                                  
                                 day2_sugi_layout_full <- full_join(d2compnetwork.nodes, day2_sugi.node.labels.df)
                                  
                                  
                                  #Build Custom Color Palette:
                                  flat.palette <- c("#F44336", "#C2185B", "#BA68C8", "#673AB7", "#7986CB", "#1E88E5", "#26C6DA", "#00897B", "#4CAF50",  "#AFB42B", "#FF8F00", "#FF5722", "#8D6E63", "#78909C")
                                  GetPalette = colorRampPalette(flat.palette)
                                  
                                  ColorCount = as.numeric(length(unique(day2_nodes$CellType)))
                                  custom.colors <- GetPalette(ColorCount)
                                  
                                  
                                  ggraph(day2_merged.graph, layout = "manual",
                                         x = day2_sugi_layout[, "x"],
                                         y = day2_sugi_layout[, "y"]) + 
                                    geom_node_point(data = day2_sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
                                    geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
                                    scale_shape_manual(values = c(15, 8, 13, 16))+
                                    scale_color_manual(values = custom.colors) +
                                    geom_edge_diagonal(aes(color = database_source), alpha = 0.35, width = 0.1, label_push = TRUE, strength = 1, arrow = arrow(angle =20, length = unit(1.5, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
                                    new_scale_color() +
                                    scale_edge_colour_manual(values = c("#2E59B7", "#B72EB7")) +
                                    theme_graph() + 
                                    theme(legend.position = "bottom")+
                                    expand_limits(y = c(min(day2_sugi.node.labels.df$y), max(day2_sugi.node.labels.df$y)+0.15*(range(day2_sugi.node.labels.df$y)[2]-range(day2_sugi.node.labels.df$y)[1]))) +
                                    geom_text_repel(data = day2_sugi.node.labels.df, aes(x = x, y = y, label = node_label), nudge_y = 0.15, size = 2.5, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5)
                                  
                                  
                                  ggraph(Day2_graph, layout = "manual",
                                         x = day2_sugi_layout_full[, "x"],
                                         y = day2_sugi_layout_full[, "y"]) + 
                                    geom_node_point(data = day2_sugi_layout, aes(shape = Basic_Celltype), color = "gray", alpha = 0.2, size = 1) + 
                                    geom_node_point(aes(color = CellType, shape = Basic_Celltype), size = 2) + 
                                    scale_shape_manual(values = c(15, 8, 13,16))+
                                    geom_edge_diagonal(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(4,'mm')) + 
                                    scale_edge_width(range = c(0.1, 0.5)) +
                                    theme_graph() + 
                                    theme(legend.position = "bottom") +
                                    facet_nodes(~Lesion) +
                                    scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
                                    expand_limits(y = c(min(day2_sugi.node.labels.df$y), max(day2_sugi.node.labels.df$y)+0.05*(range(day2_sugi.node.labels.df$y)[2]-range(day2_sugi.node.labels.df$y)[1]))) +
                                    geom_text_repel(data = day2_sugi.node.labels.df, aes(x = x, y = y, label = node_label), nudge_y = 0.15, size = 2, box.padding = 0, min.segment.length = 1.5, force = 5,  max.time = 2.5)
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                       
                                  
                                  
                                  
                                  
                                  ### IF multiple network communities:
                        
                                            ########## Community detection and Filtering to Largest Network:
                                                        
                                                        subgraphs <- igraph::clusters(merged.graph)
                                                        
                                                        #Order by size:
                                                        old_members <- subgraphs$membership
                                                        new_members <- old_members
                                                        old_sizes <- subgraphs$csize
                                                        c_order <- order(old_sizes, decreasing = T)
                                                        new_order <- seq(1:length(c_order))
                                                        for (i in 1:length(c_order)) {
                                                          int_c <- c_order[i]
                                                          int_c_nodes <- names(which(old_members == int_c))
                                                          new_members[int_c_nodes] <- new_order[i]
                                                        }
                                                        subgraphs$csize <- as.vector(table(new_members))
                                                        subgraphs$membership <- new_members
                                                        
                                                        decomposed_subgraph_old <- decompose.graph(RL_graph)
                                                        decomposed_subgraph <- vector(mode = "list", length = length(decomposed_subgraph_old))
                                                        for (i in 1:length(c_order)) {
                                                          int_decomp <- c_order[i]
                                                          decomposed_subgraph[[i]] <- decomposed_subgraph_old[[int_decomp]]
                                                        }
                                                        clusts <- as.character(as.vector(subgraphs$membership))
                                                        
                                                        merged.graph <- merged.graph %>% activate(nodes) %>% mutate(subnetwork = clusts)
                                                        
                                                        merged.nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
                                                        
                                                        layout <- full_join(layout, merged.nodes)
                                                        
                                                        RL_graph <- RL_graph %>% activate(nodes) %>% full_join(layout[,-which(colnames(layout) %in% c("x","y"))])
                                                        
                                                        graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
                                                        graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
                                                        
                                                        ggraph(RL_graph, layout = "manual",
                                                               x = layout_full[, "x"],
                                                               y = layout_full[, "y"]) + 
                                                          geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
                                                          geom_node_point(aes(color = subnetwork), size = 0.5) + 
                                                          geom_edge_link(aes(width = Scaled_Connection_Product), alpha = 0.7, arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                                                          scale_edge_width(range = c(0.1, 0.5)) +
                                                          theme_graph() + 
                                                          theme(legend.position = "bottom") +
                                                          facet_nodes(~Lesion) +
                                                          scale_edge_colour_manual(values = c("#7B7B7B", "#D20000"))
                                                        
                                                        
                                                        ## Remake Graph with Largest Network:
                                                        major.network.nodes <- graph.nodes %>% filter(subnetwork == 1) %>% as.data.frame()
                                                        
                                                        major.network.nodes.df <- filtered.nodes.df %>% filter(node_name %in% major.network.nodes$node_name)
                                                        major.network.nodes.df <- major.network.nodes.df[,-which(colnames(major.network.nodes.df) == "rowname")]
                                                        
                                                        major.network.edges.df <- filtered.edges.df %>% filter(from %in% major.network.nodes.df$node_name | to %in% major.network.nodes.df$node_name)
                                                        
                                                        RL_largest_network <- tbl_graph(nodes = major.network.nodes.df, edges = major.network.edges.df, directed = TRUE, node_key = "node_name")
                                                        
                                                        #Correct misconnected edges:
                                                        major.network.nodes <- RL_largest_network %>% activate(nodes) %>% as.data.frame()
                                                        major.network.edges <- RL_largest_network %>% activate(edges) %>% as.data.frame()
                                                        
                                                        duplicates <- major.network.edges[which(duplicated(major.network.edges$to)|duplicated(major.network.edges$from)),]  
                                
                                                        major.network.nodes.df <- rownames_to_column(major.network.nodes.df)
                                                        
                                                        for (r in 1:nrow(duplicates)) {
                                                          
                                                          les <- as.character(duplicates[r,"Lesion"])
                                                          
                                                          if (duplicates[r,"database_source"] == "CellphoneDB") {
                                                            from <- duplicates[r,"from.named"]
                                                            
                                                            to <- duplicates[r,"to.named"]
                                                            
                                                            duplicates[r,"from"] <- major.network.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                                            duplicates[r,"to"] <- major.network.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                                            
                                                          }
                                                          
                                                          if (duplicates[r,"database_source"] == "Nichenet") {
                                                            from <- duplicates[r,"from.named"]
                                                            
                                                            to <- duplicates[r,"to.named"]
                                                            
                                                            duplicates[r,"from"] <- major.network.nodes.df %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.integer()
                                                            duplicates[r,"to"] <- major.network.nodes.df %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.integer()
                                                            
                                                          }
                                                        }
                                                        
                                                        major.network.edges[which(duplicated(major.network.edges$to)|duplicated(major.network.edges$from)),] <- duplicates  
                                                        
                                                        RL_largest_network <- RL_largest_network %>% activate(edges) %>% right_join(major.network.edges)
                                                        
                                                        # Define Layout using merged data from all lesions:
                                                        merged.major.nodes <- major.network.nodes.df[,which(colnames(major.network.nodes.df) %in% c("node_name", "CellType"))]
                                                        merged.major.nodes <- merged.major.nodes[-which(duplicated(merged.major.nodes$node_name)),]  
                                                        
                                                        merged.major.graph <- tbl_graph(nodes = merged.major.nodes, edges = major.network.edges.df, directed = TRUE, node_key = "node_name")
                                                        
                                                        major.graph.layout <- layout_with_stress(merged.major.graph)
                                                        colnames(major.graph.layout) <- c("x", "y")
                                                        
                                                        major.graph.layout <- cbind(merged.major.nodes, major.graph.layout)
                                                        
                                                        
                                                        #save(RL_largest_network, file = "Data/RL_network2/temp/RL_largest_network.Rdata")
                                                        #save(major.graph.layout_full, file = "Data/RL_network2/temp/major.graph.layout_full.Rdata")
                                                        #save(major.graph.layout, file = "Data/RL_network2/temp/major.graph.layout_merged.Rdata")
                                                        #save(major.network.nodes.df, file = "Data/RL_network2/temp/major.network.nodes.df.Rdata")
                                                        #save(major.network.edges.df, file = "Data/RL_network2/temp/major.network.edges.df.Rdata")
                                                        
                                                        #load("Data/RL_network2/temp/RL_largest_network.Rdata")
                                                        #load("Data/RL_network2/temp/major.graph.layout_full.Rdata")
                                                        #load("Data/RL_network2/temp/major.graph.layout_merged.Rdata")
                                                        #load("Data/RL_network2/temp/major.network.nodes.df.Rdata")
                                                        #load("Data/RL_network2/temp/major.network.edges.df.Rdata")
                                                        
                        
                        
                                                                ## Ligand Edge Label mapping:
                                                                major.graph.label.mapping.df <- major.network.edges.df %>% filter(Edge_label != "")
                                                                major.graph.label.mapping.df$label_x <- NA
                                                                major.graph.label.mapping.df$label_y <- NA
                                                                
                                                                for (l in 1:nrow(major.graph.label.mapping.df)) {
                                                                  
                                                                  lig.from <- major.graph.label.mapping.df$from.named[l]
                                                                  lig.to <- major.graph.label.mapping.df$to.named[l]
                                                                  
                                                                  node1_info <- major.graph.layout %>% filter(node_name == lig.from)
                                                                  node2_info <- major.graph.layout %>% filter(node_name == lig.to)
                                                                  
                                                                  major.graph.label.mapping.df$label_x[l] <- (node1_info$x + node2_info$x)/2
                                                                  major.graph.label.mapping.df$label_y[l] <- (node1_info$y + node2_info$y)/2
                                                                }
                                                                
                                                                #Copy coordinates to node df containing all faceted nodes:
                                                                major.graph.layout_full <- full_join(major.network.nodes.df, major.graph.layout)
                                                             
                                                                
                                                                #Set point colors for plotting:
                                                                GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
                                                                
                                                                ColorCount = as.numeric(length(unique(major.network.nodes.df$CellType)))
                                                                node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
                                                                
                                                                
                                                                ##Plot With Facets:
                                                                ggraph(RL_largest_network, layout = "manual",
                                                                       x = major.graph.layout_full[, "x"],
                                                                       y = major.graph.layout_full[, "y"]) + 
                                                                  geom_node_point(data = major.graph.layout, color = "gray", alpha = 0.2, size = 0.5) + 
                                                                  geom_node_point(aes(color = CellType), size = 0.5) + 
                                                                  geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                                                                  scale_edge_width(range = c(0.01, 0.5)) +
                                                                  theme_graph() + 
                                                                  theme(legend.position = "bottom") +
                                                                  facet_nodes(~Lesion) +
                                                                  scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
                                                                  new_scale_color() +
                                                                  geom_text_repel(data = major.graph.label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0) +
                                                                  scale_color_manual(values = c("#000000", "#750000"))
                                                                
            
                                                                
                                                                
                                                                
              #Test new layout method:          
              
                        
                        # Define Layout using merged data from all lesions:
                        merged.major.nodes <- major.network.nodes.df[,which(colnames(major.network.nodes.df) %in% c("node_name", "CellType"))]
                        merged.major.nodes <- merged.major.nodes[-which(duplicated(merged.major.nodes$node_name)),]  
                        
                        merged.major.nodes$Node_lesions <- NA

                        for (r in 1:nrow(merged.major.nodes)) {
                          
                          node <- merged.major.nodes$node_name[r]
                          
                          node.edges <- major.network.edges.df %>% filter(from.named == node | to.named == node)
                          
                          row.node.lesions <- str_c(unique(as.character(node.edges$Lesion)), collapse = "/")
                          
                          merged.major.nodes[which(merged.major.nodes$node_name == node), "Node_lesions"] <- row.node.lesions
                          
                        }
                        
                        #Setup columns for node categorization:
                        merged.major.nodes$Node_lesion_length <- sapply(merged.major.nodes$Node_lesions, FUN = function(x) length(unlist(str_split(x, pattern = "/"))))
                        merged.major.nodes$NL_included <- str_detect(merged.major.nodes$Node_lesions, pattern = "Nonlesional") 
                        merged.major.nodes <- merged.major.nodes %>% mutate(allergy_only = case_when(str_detect(Node_lesions, pattern = "Allergy") == TRUE & str_detect(Node_lesions, pattern = "Irritant") == FALSE & str_detect(Node_lesions, pattern = "Nonlesional") == FALSE ~ TRUE))
                        merged.major.nodes$allergy_only[which(is.na(merged.major.nodes$allergy_only))] <- FALSE
                        
                       
                        merged.major.nodes <- merged.major.nodes %>% mutate(Lesion_category = case_when(Node_lesion_length == 1 ~ Node_lesions,
                                                                            Node_lesion_length > 1 & NL_included == TRUE ~ "Nonlesional_shared",
                                                                            Node_lesion_length > 1 & NL_included == FALSE & allergy_only==FALSE ~ "Inflammatory_shared",
                                                                            Node_lesion_length > 1 & NL_included == FALSE & allergy_only==TRUE ~ "Allergy_specific"))
                        
                        merged.major.nodes$Lesion_category[which(str_detect(merged.major.nodes$Lesion_category, pattern = "Allergy"))] <- "Allergy_specific"
                        
                        unique(merged.major.nodes$Lesion_category)
                        
                        merged.major.nodes$Lesion_category <- factor(merged.major.nodes$Lesion_category,
                                                               levels = c("Nonlesional", "Nonlesional_shared", "Irritant", "Inflammatory_shared","Allergy_specific"))
                        
                        
                        
                        merged.major.graph <- tbl_graph(nodes = merged.major.nodes, edges = major.network.edges.df, directed = TRUE, node_key = "node_name")
                        
                        #Try new backbone R package:
                        library(backbone)

                        ig.object <- as.igraph(merged.major.graph)
                        
                      
                        plot(ig.object, vertex.label = NA, vertex.size = 2, edge.curved=0.2, edge.arrow.size= 0.2, edge.width=1)
                        
                        bb <- disparity(ig.object, alpha = 0.5, narrative = TRUE, class = "igraph")
                        
                        plot(bb, vertex.label = NA, vertex.size = 2, edge.curved=0.2, edge.arrow.size= 0.2, edge.width=1)
                        
                        bb <- osdsm(ig.object, alpha = 0.05, narrative = TRUE, #An oSDSM backbone...
                              class = "igraph", trials = 100)
                        
                        #A weighted binary bipartite network of 20 agents & 50 artifacts; agents form two communities
                        B <- rbind(cbind(matrix(sample(0:3, 250, replace = TRUE, prob = ((1:4)^2)),10),
                                         matrix(sample(0:3, 250, replace = TRUE, prob = ((4:1)^2)),10)),
                                   cbind(matrix(sample(0:3, 250, replace = TRUE, prob = ((4:1)^2)),10),
                                         matrix(sample(0:3, 250, replace = TRUE, prob = ((1:4)^2)),10)))


                        P <- B%*%t(B) #An ordinary weighted projection...
                        
                        P
                        plot(igraph::graph_from_adjacency_matrix(P, mode = "undirected",
                                                                 weighted = TRUE, diag = FALSE)) #...is a dense hairball
                        bb <- osdsm(B, alpha = 0.05, narrative = TRUE, #An oSDSM backbone...
                                    class = "igraph", trials = 100)
                        plot(bb) #...is sparse with clear communities
                        
                        
                        
                        major.graph.layout <- layout_with_stress(merged.major.graph)
                        colnames(major.graph.layout) <- c("x", "y")
                        
                        major.graph.layout <- cbind(merged.major.nodes, major.graph.layout)
                        
                        
                        #Copy coordinates to node df containing all faceted nodes:
                        major.graph.layout_full <- full_join(major.network.nodes.df, major.graph.layout)
                        
                        #Set point colors for plotting:
                        GetPalette = colorRampPalette(c("#FF0000","#FF7000", "#1BB018", "#187BB0", "#B64FE3"))
                        
                        ColorCount = as.numeric(length(unique(filtered.nodes.df$CellType)))
                        node.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
                        
                        ##Plot With Facets:
                        ggraph(RL_largest_network, layout = "manual",
                               x = major.graph.layout_full[, "x"],
                               y = major.graph.layout_full[, "y"]) + 
                          geom_node_point(data = major.graph.layout, color = "gray", alpha = 0.2, size = 0.5) + 
                          geom_node_point(aes(color = CellType), size = 0.5) + 
                          geom_edge_link(aes(color = DE, width = Scaled_Connection_Product), alpha = 0.7, label_push = TRUE ,arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                          scale_edge_width(range = c(0.01, 0.5)) +
                          theme_graph() + 
                          theme(legend.position = "bottom") +
                          facet_nodes(~Lesion) +
                          scale_edge_colour_manual(values = c("#7B7B7B", "#D20000")) +
                          new_scale_color() +
                          geom_text_repel(data = major.graph.label.mapping.df, aes(x = label_x, y = label_y, label = Edge_label, color = DE), size = 2.5, alpha =1, box.padding = 0) +
                          scale_color_manual(values = c("#000000", "#750000"))
                        
                        #graphjs(stress.graph, bg="black", vertex.shape="sphere")
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
   #Receptor ligand heatmaps:
        ligs <- network %>% dplyr::pull(Ligand) %>% unique()
        
            complex.ligs <- deconvoluted %>% filter(complex_name %in% ligs) %>% dplyr::pull(complex_name) %>% unique()
            complex.lig.genes <- deconvoluted %>% filter(complex_name %in% ligs) %>% dplyr::pull(gene_name) %>% unique()
            
            if (length(complex.ligs)>0) {
              ligs <- ligs[-which(ligs %in% complex.ligs)]
              ligs <- unique(c(ligs, complex.lig.genes))
            }
            
        recs <- network %>% dplyr::pull(Receptor) %>% unique()
            
            complex.recs <- deconvoluted %>% filter(complex_name %in% recs) %>% dplyr::pull(complex_name) %>% unique()
            complex.rec.genes <- deconvoluted %>% filter(complex_name %in% recs) %>% dplyr::pull(gene_name) %>% unique()
            
            if (length(complex.recs)>0) {
              recs <- recs[-which(recs %in% complex.recs)]
              recs <- unique(c(recs, complex.rec.genes))
            }

        #By Lesion and Celltype:
        DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Semi_Detailed_Celltype", "Lesion"))
        
        
            #Ligands Heatmaps:
        
        
                          #Scaled by Celltype:
                          set.seed(100)
                          lig.heatmap <- plot_heatmap_srt(DE.input, 
                                                          genes = ligs, 
                                                          type = "bulk", 
                                                          facet_by = "Lesion",
                                                          scale_group = "Semi_Detailed_Celltype",
                                                          cluster_by = "row", 
                                                          pdf_format = "tile", 
                                                          scale_by = "row",
                                                          text_angle = 90,
                                                          cluster_type = "kmeans", 
                                                          k = 8, 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = F,
                                                          floor = F)
                          
                          
                          #Gene order:
                          lig_genes <- c(7,4,     6,8)
                          lig_genes <- rev(lig_genes)
                          lig_genes_reorder <- c()
                          
                          for (i in 1:length(lig_genes)) {
                            int1 <- lig_genes[i]
                            ind <- grep(paste0("^", int1, "$"), lig.heatmap[[2]]$cluster)
                            reorder <- names(lig.heatmap[[2]]$cluster)[ind]
                            lig_genes_reorder <- c(lig_genes_reorder, reorder)
                          }
                          
                          lig.heatmap[[2]]$cluster[which(names(lig.heatmap[[2]]$cluster) == "CXCL14")]
                          
                          #Heatmap Gene breaks:
                          lig_cluster_breaks <- c(6)
                          lig_cluster_breaks <- rev(lig_cluster_breaks)
                          
                          lig_gene_breaks <- find_gene_breaks(lig.heatmap, lig_genes_reorder,lig_cluster_breaks)
                          
                          #Manually add CXCL14, IFNG, and IL13 back:
                          lig_genes_reorder <- c("IFNG","IL13", lig_genes_reorder[1:which(lig_genes_reorder== "MST1")], lig_genes_reorder[c(which(lig_genes_reorder== "MST1")+1):length(lig_genes_reorder)],"CXCL14")
                          
                          ligand_labels <- lig_genes_reorder[which(str_detect(lig_genes_reorder, pattern = "^IL|^CCL|^CXCL|^LTB|^TNF|^IFN|^OSM"))]
                          
                          plot_heatmap_srt(DE.input, 
                                           genes = lig_genes_reorder, 
                                           gene_breaks = lig_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = ligand_labels,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Lesion",
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
                          
                          ligand.plot.celltype.scaled <- last_plot()
                          
                          set.seed(100)
                          ligand.plot.celltype.scaled.full <- plot_heatmap_srt(DE.input, 
                                                                genes = ligs, 
                                                                type = "bulk", 
                                                                facet_by = "Lesion",
                                                                scale_group = "Semi_Detailed_Celltype",
                                                                cluster_by = "row", 
                                                                pdf_format = "tile", 
                                                                scale_by = "row",
                                                                text_angle = 90,
                                                                show_k = T,
                                                                color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                                log_scale_base = 10,
                                                                ceiling = F,
                                                                floor = F)
                          
                          #Not scaled by Celltype:
                          set.seed(100)
                          lig.heatmap <- plot_heatmap_srt(DE.input, 
                                                          genes = ligs, 
                                                          type = "bulk", 
                                                          facet_by = "Semi_Detailed_Celltype",
                                                          cluster_by = "row", 
                                                          pdf_format = "tile", 
                                                          scale_by = "row",
                                                          text_angle = 90,
                                                          cluster_type = "kmeans", 
                                                          k = 12, 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = 2,
                                                          floor = -2)
                          
                          
                          #Gene order:
                          lig_genes <- c(5,   12,8,11,   7,9,  4,  1,    3,    6,  2, 10  )
                          lig_genes <- rev(lig_genes)
                          lig_genes_reorder <- c()
                          
                          for (i in 1:length(lig_genes)) {
                            int1 <- lig_genes[i]
                            ind <- grep(paste0("^", int1, "$"), lig.heatmap[[2]]$cluster)
                            reorder <- names(lig.heatmap[[2]]$cluster)[ind]
                            lig_genes_reorder <- c(lig_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          lig_cluster_breaks <- c(12,7,4,1,3,6)
                          lig_cluster_breaks <- rev(lig_cluster_breaks)
                          
                          lig_gene_breaks <- find_gene_breaks(lig.heatmap, lig_genes_reorder,lig_cluster_breaks)
                          
                          
                          plot_heatmap_srt(DE.input, 
                                           genes = lig_genes_reorder, 
                                           gene_breaks = lig_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = 2,
                                           floor = -2,
                                           panel_spacing = 0.1,
                                           plot_axis = TRUE)
                          
                          ligand.plot.unscaled <- last_plot()
            
                          
            #Receptor Heatmaps:
                          
                          
                          #Scaled by Celltype:                    
                          set.seed(100)
                          rec.heatmap <- plot_heatmap_srt(DE.input, 
                                                          genes = recs, 
                                                          type = "bulk", 
                                                          facet_by = "Lesion",
                                                          scale_group = "Semi_Detailed_Celltype",
                                                          cluster_by = "row", 
                                                          pdf_format = "tile", 
                                                          scale_by = "row",
                                                          text_angle = 90,
                                                          cluster_type = "kmeans", 
                                                          k = 8, 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = F,
                                                          floor = F)
                          
                          #Gene order:
                          rec_genes <- c(1,7,   8)
                          rec_genes <- rev(rec_genes)
                          rec_genes_reorder <- c()
                          
                          for (i in 1:length(rec_genes)) {
                            int1 <- rec_genes[i]
                            ind <- grep(paste0("^", int1, "$"), rec.heatmap[[2]]$cluster)
                            reorder <- names(rec.heatmap[[2]]$cluster)[ind]
                            rec_genes_reorder <- c(rec_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          rec_cluster_breaks <- c(8)
                          rec_cluster_breaks <- rev(rec_cluster_breaks)
                          
                          rec_gene_breaks <- find_gene_breaks(rec.heatmap, rec_genes_reorder,rec_cluster_breaks)

                          receptor_labels <- rec_genes_reorder[which(str_detect(rec_genes_reorder, pattern = "^IL|^CCR|^CXCR|^LTB|^TNF|^IFN|^OSM"))]

                          plot_heatmap_srt(DE.input, 
                                           genes = rec_genes_reorder, 
                                           gene_breaks = rec_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = receptor_labels,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Lesion",
                                           scale_group = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#000000","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = F,
                                           floor = F,
                                           panel_spacing = 1,
                                           plot_axis = TRUE)
                          
                          receptor.plot.celltype.scaled <- last_plot()
                          
                          set.seed(100)
                          receptor.plot.celltype.scaled.full <- plot_heatmap_srt(DE.input, 
                                                                               genes = recs, 
                                                                               type = "bulk", 
                                                                               facet_by = "Lesion",
                                                                               scale_group = "Semi_Detailed_Celltype",
                                                                               cluster_by = "row", 
                                                                               pdf_format = "tile", 
                                                                               scale_by = "row",
                                                                               text_angle = 90,
                                                                               show_k = T,
                                                                               color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                                               log_scale_base = 10,
                                                                               ceiling = F,
                                                                               floor = F)
        
                          
                          #Not Scaled by Celltype:                    
                          set.seed(100)
                          rec.heatmap <- plot_heatmap_srt(DE.input, 
                                                          genes = recs, 
                                                          type = "bulk", 
                                                          facet_by = "Semi_Detailed_Celltype",
                                                          cluster_by = "row", 
                                                          pdf_format = "tile", 
                                                          scale_by = "row",
                                                          text_angle = 90,
                                                          cluster_type = "kmeans", 
                                                          k = 12, 
                                                          show_k = T,
                                                          color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                          log_scale_base = 10,
                                                          ceiling = 2,
                                                          floor = -2)
                          
                          #Gene order:
                          rec_genes <- c(12,8,9,    3,    2,6,    1,5,7,   4,10,11)
                          rec_genes <- rev(rec_genes)
                          rec_genes_reorder <- c()
                          
                          for (i in 1:length(rec_genes)) {
                            int1 <- rec_genes[i]
                            ind <- grep(paste0("^", int1, "$"), rec.heatmap[[2]]$cluster)
                            reorder <- names(rec.heatmap[[2]]$cluster)[ind]
                            rec_genes_reorder <- c(rec_genes_reorder, reorder)
                          }
                          
                          #Heatmap Gene breaks:
                          rec_cluster_breaks <- c(3,2,1,4)
                          rec_cluster_breaks <- rev(rec_cluster_breaks)
                          
                          rec_gene_breaks <- find_gene_breaks(rec.heatmap, rec_genes_reorder,rec_cluster_breaks)
                          
                          plot_heatmap_srt(DE.input, 
                                           genes = rec_genes_reorder, 
                                           gene_breaks = rec_gene_breaks,
                                           gene_names = TRUE, 
                                           gene_labels = NULL,
                                           gene_label_side = "left",
                                           gene_labels_size = 3,
                                           gene_labels_nudge = 0.3,
                                           type = "bulk", 
                                           facet_by = "Semi_Detailed_Celltype",
                                           cluster_by = FALSE, 
                                           pdf_format = "tile", 
                                           scale_by = "row",
                                           text_angle = 90,
                                           color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#000000","#000000","#792B45","#EC1F63"))(256),
                                           log_scale_base = 10,
                                           ceiling = 2,
                                           floor = -2,
                                           panel_spacing = 0.1,
                                           plot_axis = TRUE)
                          
                          receptor.plot.unscaled <- last_plot()
        
        
        
        
        plot_heatmap_srt(DE.input,
                         genes = recs,
                         type = "bulk",
                         facet_by = c("Lesion"),
                         scale_group = "Semi_Detailed_Celltype",
                         text_angle = 90,
                         text_sizes = c(20, 10, 7, 10, 5, 5, 5),
                         color_pal = colorRampPalette(c("#005775","#003A4E","#000000","#000000","#000000","#792B45","#EC1F63"))(256),    
                         log_scale_base = 2,
                         ceiling = 2,
                         floor = -2)   
                        

      ########## Community detection and Gene/Cell-specific plots:
                        
              subgraphs <- igraph::clusters(merged.graph)
              
              #Order by size:
              old_members <- subgraphs$membership
              new_members <- old_members
              old_sizes <- subgraphs$csize
              c_order <- order(old_sizes, decreasing = T)
              new_order <- seq(1:length(c_order))
              for (i in 1:length(c_order)) {
                int_c <- c_order[i]
                int_c_nodes <- names(which(old_members == int_c))
                new_members[int_c_nodes] <- new_order[i]
              }
              subgraphs$csize <- as.vector(table(new_members))
              subgraphs$membership <- new_members
              
              decomposed_subgraph_old <- decompose.graph(RL_graph)
              decomposed_subgraph <- vector(mode = "list", length = length(decomposed_subgraph_old))
              for (i in 1:length(c_order)) {
                int_decomp <- c_order[i]
                decomposed_subgraph[[i]] <- decomposed_subgraph_old[[int_decomp]]
              }
              clusts <- as.character(as.vector(subgraphs$membership))
              
              merged.graph <- merged.graph %>% activate(nodes) %>% mutate(subnetwork = clusts)
              
              merged.nodes <- merged.graph %>% activate(nodes) %>% as.data.frame()
              
              layout <- full_join(layout, merged.nodes)
              
              RL_graph <- RL_graph %>% activate(nodes) %>% full_join(layout[,-which(colnames(layout) %in% c("x","y"))])
              
              graph.nodes <- RL_graph %>% activate(nodes) %>% as.data.frame()
              graph.edges <- RL_graph %>% activate(edges) %>% as.data.frame()
              
              ggraph(RL_graph, layout = "manual",
                     x = layout_full[, "x"],
                     y = layout_full[, "y"]) + 
                geom_node_point(data = layout, color = "gray", alpha = 0.2, size = 0.5) + 
                geom_node_point(aes(color = subnetwork), size = 0.5) + 
                geom_edge_link(aes(width = Scaled_Connection_Product), alpha = 0.7, arrow = arrow(angle =20, length = unit(1, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                scale_edge_width(range = c(0.1, 0.5)) +
                theme_graph() + 
                theme(legend.position = "bottom") +
                facet_nodes(~Lesion) +
                scale_edge_colour_manual(values = c("#7B7B7B", "#D20000"))
              
              major.network <- graph.nodes %>% filter(subnetwork == 1)
              major.network <- graph.edges %>% filter(from.named %in% major.network$node_name | to.named %in% major.network$node_name)
              
              RL_graph <- RL_graph %>% activate(nodes) %>% full_join(layout[,-which(colnames(layout) %in% c("x","y"))])
              
              
              
           
            
            
        #Filter to subnetwork communities with gene of interest:
              
              Desired_genes_or_celltypes <- c("IFNG")    
              
              Desired_subnetworks <- c()
              for (i in 1:length(Desired_genes_or_celltypes)) {
                
                associated.subnetworks <- graph.nodes[str_detect(string = graph.nodes$node_name, pattern = Desired_genes_or_celltypes[i]),] %>% select(subnetwork) %>% unlist() %>% as.character()
                
                Desired_subnetworks <- unique(c(Desired_subnetworks, associated.subnetworks))
                
              }
              
              subnet.nodes <- graph.nodes %>% filter(subnetwork %in% c(Desired_subnetworks)) %>% as.data.frame()
              
              subnet.edges <- edges.df %>% filter(from %in% subnet.nodes$node_name | to %in% subnet.nodes$node_name) %>% as.data.frame() 
              
              Filtered_graph <- tbl_graph(nodes = subnet.nodes, edges = subnet.edges, directed = TRUE, node_key = "node_name")
              
              #Correct misconnected edges:
              
              filtered.graph.nodes <- Filtered_graph %>% activate(nodes) %>% as.data.frame()
              filtered.graph.edges <- Filtered_graph %>% activate(edges) %>% as.data.frame()
              
              duplicates <- filtered.graph.edges[which(duplicated(filtered.graph.edges$to)|duplicated(filtered.graph.edges$from)),]  
              
              subnet.nodes <- rownames_to_column(subnet.nodes)
              
              for (r in 1:nrow(duplicates)) {
                
                les <- as.character(duplicates[r,"Lesion"])
                
                from <- paste(duplicates[r,"Ligand_Source"], ":", duplicates[r, "Ligand"], sep = "")
                
                to <- paste(duplicates[r,"Receptor_source"], ":", duplicates[r, "Receptor"], sep = "")
                
                duplicates[r,"from"] <- subnet.nodes %>% filter(Lesion == les) %>% filter(node_name == from) %>% select(rowname) %>% as.numeric()
                duplicates[r,"to"] <- subnet.nodes %>% filter(Lesion == les) %>% filter(node_name == to) %>% select(rowname) %>% as.numeric()
                
              }
              
              filtered.graph.edges[which(duplicated(filtered.graph.edges$to)|duplicated(filtered.graph.edges$from)),] <- duplicates  
              
              Filtered_graph <- Filtered_graph %>% activate(edges) %>% right_join(filtered.graph.edges)
              
              
              # Define Layout using merged data from all lesions:
              nodes <- filtered.graph.nodes[,which(colnames(filtered.graph.nodes) %in% c("node_name", "CellType"))]
              nodes <- nodes[-which(duplicated(nodes$node_name)),]  
              
              filtered.merged.graph <- tbl_graph(nodes = nodes, edges = subnet.edges, directed = TRUE, node_key = "node_name")
              
              filtered.layout <- layout_with_stress(filtered.merged.graph)
              colnames(filtered.layout) <- c("x", "y")
              
              filtered.layout <- cbind(nodes, filtered.layout)
              
              #Copy coordinates to node df containing all faceted nodes:
              filtered.layout_full <- full_join(filtered.graph.nodes, filtered.layout)
              
              
              ##Plot With Facets:
              ggraph(Filtered_graph, layout = "manual",
                     x = filtered.layout_full[, "x"],
                     y = filtered.layout_full[, "y"]) + 
                geom_node_point(data = filtered.layout, color = "gray", alpha = 0.1, size = 3) + 
                geom_edge_link(aes(width = Scaled_Connection_Product), alpha = 0.7, arrow = arrow(angle =20, length = unit(0.5, "mm"), type = "closed"), end_cap = circle(0.5,'mm')) + 
                geom_node_point(aes(color = CellType), size = 0.5) + 
                scale_edge_width(range = c(0.0001, 0.5)) +
                theme_graph() + 
                theme(legend.position = "bottom") +
                facet_nodes(~Lesion)
                #geom_node_label(aes(label = node_name), size = 1, repel = TRUE, point.padding = unit(0.2, "lines"), fill = "white")
              
                  
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                                              #test backbone layout:
                            
                                                    simplified.graph <- igraph::simplify(merged.graph)
                                                    simplified.graph <- igraph::as.undirected(simplified.graph)
                                                    
                                                    bb <- layout_as_backbone(simplified.graph, keep = 0.4)
                                                    
                                                    E(simplified.graph)$col <- F
                                                    E(simplified.graph)$col[bb$backbone] <- T
                                                    
                                                    ggraph(simplified.graph,
                                                           layout = "manual",
                                                           x = bb$xy[, 1],
                                                           y = bb$xy[, 2]) +
                                                      geom_edge_link0(aes(col = col), width = 0.2) +
                                                      geom_node_point(shape = 21, size = 3) +
                                                      scale_color_brewer(palette = "Set1") +
                                                      scale_fill_brewer(palette = "Set1") +
                                                      scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
                                                      theme_graph()+
                                                      theme(legend.position = "none")
                                                    
                        
  #Trim CellphoneDB interactions based on nichenet interactions:

                                
            
            
                