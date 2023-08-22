#############################################################
### Last Filter:
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

        
setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load data:
load(file = "Data/contact_derm_filtered3.Rdata")

contact_derm@meta.data$Sequence_instrument <- factor(contact_derm@meta.data$Sequence_instrument,
                                                     levels = c("NS500602", "NB501205", "VH00230", "A00439", "A00197"))

contact_derm@meta.data$Lesion <- factor(contact_derm@meta.data$Lesion,
                                                     levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))

#Last cell filter:
    lymphs.to.filter <- contact_derm@meta.data %>% filter(Basic_Celltype == "LYMPH") %>% filter(!seurat_clusters %in% c(1,8)) %>% rownames()
    KRT.to.filter <- contact_derm@meta.data %>% filter(Basic_Celltype == "KRT") %>% filter(!seurat_clusters %in% c(2,14,4,13,6,20,11,9,17,10,5,7,18,12)) %>% rownames()
    MEL.to.filter <- contact_derm@meta.data %>% filter(Basic_Celltype == "MEL") %>% filter(!seurat_clusters %in% c(16,15)) %>% rownames()
    APC.to.filter <- contact_derm@meta.data %>% filter(Basic_Celltype == "APC") %>% filter(!seurat_clusters %in% c(0,3,19,21)) %>% rownames()
    
    cell.to.filter <- c(lymphs.to.filter, KRT.to.filter, MEL.to.filter, APC.to.filter)

    contact_derm <- contact_derm[,-which(colnames(contact_derm) %in% cell.to.filter)]

    DefaultAssay(contact_derm) <- "RNA"
    
    #save(contact_derm, file = "Data/contact_derm_filtered3.Rdata")
    
#############################################################
### Final Clustering
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

plot_tsne_metadata_srt(contact_derm, color_by = "Basic_Celltype", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)

plot_tsne_gene_srt(contact_derm, gene = "TRDC")

plot_tsne_metadata_srt(contact_derm, color_by = "seurat_clusters", facet_by = c("Lesion", "Sequence_instrument"),size = 0.25, plot_labels = FALSE, plot_legend = FALSE, facet_grid = TRUE)

table(contact_derm@meta.data$Sequence_date)
              
#Look at Clusters:  
    plot_density_ridge_srt(contact_derm, val = "UMI_sum_raw", color_by = "Basic_Celltype") + geom_vline(xintercept = 300)
    
    plot_tsne_metadata_srt(contact_derm, color_by = "seurat_clusters", facet_by = c("Disease", "Lesion"), facet_grid = TRUE, plot_legend = FALSE, plot_labels = FALSE)
    plot_tsne_metadata_srt(contact_derm, color_by = "seurat_clusters", size = 0.5, plot_labels = TRUE, plot_legend = FALSE)
    plot_tsne_metadata_srt(contact_derm, color_by = "seurat_clusters", facet_by = "seurat_clusters", size = 0.5, plot_labels = FALSE, plot_legend = FALSE)
    
    
    plot_tsne_metadata_srt(DE.input, color_by = "Patient",  plot_legend = TRUE, plot_labels = FALSE)
    
    #Cluster Assignment:
    load("gene_lists/established_exclusion_genes_5.17.22.Rdata")
    
    all.markers <- fread("Data/contact_derm_filtered.markers.tsv")
    top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    
    Filtered.markers <- all.markers %>% filter(avg_log2FC > 0 & pct.2 < 0.6 & !gene %in% established_exclusion_genes)
    top10_FC_pct2 <- Filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC*pct.1)/(pct.2 + 0.00000001))
    
    Idents(contact_derm) <- "seurat_clusters"
    
    levels(contact_derm)
    
    new.cluster.ids <- c("LC", "Lymph1", "Krt1", "Mac1", "Krt2", "Krt3", "Krt4", "Krt5",
                         "Lymph2", "Krt6", "Krt7", "Krt8", "Krt9", "Krt10", "Krt11", "Mel1", "Mel2",
                         "Krt12", "Krt13", "Mac2", "Krt14", "Bcell")
    
    names(new.cluster.ids) <- levels(contact_derm)
    contact_derm <- RenameIdents(contact_derm, new.cluster.ids)
    DimPlot(contact_derm, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()    
    
    contact_derm@meta.data$Detailed_clusters <- as.character(contact_derm@active.ident)    
    
    
          #Subset KRTs from pre-integrated object to adjust their clusters:
                  Idents(contact_derm) <- "Basic_Celltype"
                  
                  KRTs <- subset(contact_derm, idents = "KRT")
                  
                  Idents(KRTs) <- "seurat_clusters"
                  DefaultAssay(KRTs) <- "integrated"
                  
                  KRTs <- RunPCA(KRTs, verbose = FALSE, npcs = 75)
                  
                  KRTs <- FindNeighbors(KRTs, dims = 1:75, k.param = 30)
                  KRTs <- FindClusters(KRTs, resolution = 0.75)
                  KRTs <- RunUMAP(KRTs, dims = 1:30, n.neighbors = 30L)
                  
                  DimPlot(KRTs, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()
                  
                  #Markers:
                  DefaultAssay(KRTs) <- "RNA"
                  
                  load("gene_lists/established_exclusion_genes_5.17.22.Rdata")
                  
                  KRT.markers <- FindAllMarkers(KRTs, min.pct = 0.25, logfc.threshold = 0.25)
                  top10 <- KRT.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                  
                  Filtered.KRT.markers <- KRT.markers %>% filter(avg_log2FC > 0 & pct.2 < 0.6 & !gene %in% established_exclusion_genes)
                  KRT_top10_FC_pct2 <- Filtered.KRT.markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC*pct.1)/(pct.2 + 0.00000001))
                
                  plot_violin_srt(KRTs, gene = "KRT14", color_by = "seurat_clusters")
                  
                  plot_tsne_gene_srt(KRTs, gene = "KRT14", log_scale = T)
      
                  plot_density_ridge_srt(KRTs, val = "nFeature_RNA", color_by = "seurat_clusters")
                  
                  plot_tsne_metadata_srt(KRTs, color_by = "seurat_clusters", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
                  plot_tsne_metadata_srt(KRTs, color_by = "seurat_clusters", facet_by = "Sequence_instrument", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
                  
                  #Rename Clusters:  
                  new.cluster.ids <- c("KRT-b2", "KRT-b1", "KRT-sp", "KRT-sp", "KRT-sp", "KRT-wr", "KRT-g", "KRT-wr", "KRT-b2", "KRT-b2", "KRT-b3", "KRT-b3", "KRT-sp", "KRT-mucl", "KRT-sp", "KRT-b3", "KRT-b2")
      
                  names(new.cluster.ids) <- levels(KRTs)
                  KRTs <- RenameIdents(KRTs, new.cluster.ids)
                  DimPlot(KRTs, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()    
                  
                  KRTs@meta.data$Detailed_clusters <- as.character(KRTs@active.ident)
                  
                  plot_tsne_metadata_srt(KRTs, color_by = "Detailed_clusters", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
                  plot_tsne_metadata_srt(KRTs, color_by = "Detailed_clusters", facet_by = "Lesion", size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                  
                  #save(KRTs, file = "Data/KRTs_subset.Rdata")
                  #write.table(KRT.markers, file = "Data/KRT.markers.tsv", row.names=FALSE, sep="\t")
                  
                  
          #Cluster Lymphs on their own for more detailed clusters:
          
                    #Subset lymphs from pre-integrated object:
                        Idents(contact_derm) <- "Basic_Celltype"
          
                        Lymphs <- subset(contact_derm, idents = "LYMPH")
                        
                        Idents(Lymphs) <- "seurat_clusters"
                        DefaultAssay(Lymphs) <- "integrated"
                        
                        
                        Lymphs <- RunPCA(Lymphs, verbose = FALSE, npcs = 75)
                        
                        Lymphs <- FindNeighbors(Lymphs, dims = 1:75, k.param = 30)
                        Lymphs <- FindClusters(Lymphs, resolution = 1.0)
                        Lymphs <- RunUMAP(Lymphs, dims = 1:30, n.neighbors = 30L)
                        
                        DimPlot(Lymphs, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()
                        
                        
                        #View integration:
                        plot_tsne_metadata_srt(Lymphs, color_by = "seurat_clusters", facet_by = c("Sequence_instrument"), size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                        plot_tsne_metadata_srt(Lymphs, color_by = "seurat_clusters", facet_by = c("Lesion", "Sequence_instrument"),size = 0.25, plot_labels = FALSE, plot_legend = FALSE, facet_grid = TRUE)
                        plot_tsne_metadata_srt(Lymphs, color_by = "seurat_clusters", facet_by = c("Lesion"),size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                        
                        #Markers:
                        DefaultAssay(Lymphs) <- "RNA"
                        
                        load("gene_lists/established_exclusion_genes_5.17.22.Rdata")
                        
                        Lymph.markers <- FindAllMarkers(Lymphs, min.pct = 0.25, logfc.threshold = 0.25)
                        top10 <- Lymph.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                        
                        Filtered.Lymph.markers <- Lymph.markers %>% filter(avg_log2FC > 0 & pct.2 < 0.6 & !gene %in% established_exclusion_genes)
                        Lymph_top10_FC_pct2 <- Filtered.Lymph.markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC*pct.1)/(pct.2 + 0.00000001))
                        
                
                        plot_violin_srt(Lymphs, gene = "RORC", color_by = "seurat_clusters")
                        
                        plot_tsne_gene_srt(Lymphs, gene = "CD8B", log_scale = T)
      
                        plot_tsne_metadata_srt(Lymphs, color_by = "seurat_clusters", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
                        
                      
                        #Rename Clusters:  
                        new.cluster.ids <- c("CD4_conv1", "CD8", "Treg", "CD4_tcm", "CD4_cd161", "CD4_trm", "CD4_conv1", "CD4_hsp", "NK_ctl", "NK_res", "CD8_CD4_active", "CD4_conv1")
                        names(new.cluster.ids) <- levels(Lymphs)
                        Lymphs <- RenameIdents(Lymphs, new.cluster.ids)
                        DimPlot(Lymphs, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()    
                        
                        Lymphs@meta.data$Detailed_clusters <- as.character(Lymphs@active.ident)
                        
                        plot_tsne_metadata_srt(Lymphs, color_by = "Detailed_clusters", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
                        plot_tsne_metadata_srt(Lymphs, color_by = "Detailed_clusters", facet_by = "Lesion", size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                        
                        #save(Lymphs, file = "Data/lymphs_subset.Rdata")
                        #write.table(Lymph.markers, file = "Data/Lymph.markers.tsv", row.names=FALSE, sep="\t")
        
                        
              #Subset APCs from pre-integrated object:
                        Idents(contact_derm) <- "Basic_Celltype"
                        
                        APCs <- subset(contact_derm, idents = "APC")
                        
                        Idents(APCs) <- "seurat_clusters"
                        DefaultAssay(APCs) <- "integrated"
                        
                        
                        APCs <- RunPCA(APCs, verbose = FALSE, npcs = 75)
                        
                        APCs <- FindNeighbors(APCs, dims = 1:75, k.param = 30)
                        APCs <- FindClusters(APCs, resolution = 1.25)
                        APCs <- RunUMAP(APCs, dims = 1:30, n.neighbors = 30L)
                        
                        DimPlot(APCs, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()
                        
                        #View integration:
                        plot_tsne_metadata_srt(APCs, color_by = "seurat_clusters", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
                        
                        plot_tsne_metadata_srt(APCs, color_by = "seurat_clusters", facet_by = c("Sequence_instrument"), size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                        plot_tsne_metadata_srt(APCs, color_by = "seurat_clusters", facet_by = c("Lesion", "Sequence_instrument"),size = 0.25, plot_labels = FALSE, plot_legend = FALSE, facet_grid = TRUE)
                        plot_tsne_metadata_srt(APCs, color_by = "seurat_clusters", facet_by = c("Lesion"),size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                        
                        #Markers:
                        DefaultAssay(APCs) <- "RNA"
                        
                        load("gene_lists/established_exclusion_genes_5.17.22.Rdata")
                        
                        APC.markers <- FindAllMarkers(APCs, min.pct = 0.25, logfc.threshold = 0.25)
                        top10 <- APC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                        
                        Filtered.APC.markers <- APC.markers %>%  filter(!gene %in% established_exclusion_genes) %>% filter(avg_log2FC > 50 | avg_log2FC > 2 & pct.2 < 0.075)
                        APC_top10_FC_pct2 <- Filtered.APC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC/(pct.2 + 0.00000001))
                        
                        plot_violin_srt(APCs, gene = "LYZ", color_by = "seurat_clusters")
                        
                        plot_tsne_gene_srt(APCs, gene = "LYZ", log_scale = T)
                        
                        plot_tsne_metadata_srt(APCs, color_by = "seurat_clusters", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
                        
                        Idents(APCs) <- "seurat_clusters"
                        
                        #Rename Clusters:  
                        new.cluster.ids <- c("LC", "Myeloid", "mixed","Myeloid", "mixed", "LC", "mixed", "Myeloid_ccl22", "Myeloid", "M2_mac", "Bcell", "LC", "cDC1")
                        names(new.cluster.ids) <- levels(APCs)
                        APCs <- RenameIdents(APCs, new.cluster.ids)
                        DimPlot(APCs, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()    
                        
                        APCs@meta.data$Detailed_clusters <- as.character(APCs@active.ident)
                      
                                #Recluster mixed clusters:
                                  Idents(APCs) <- "Detailed_clusters"
                                  
                                  APCs2 <- subset(APCs, idents = "mixed")
                                  
                                  Idents(APCs2) <- "seurat_clusters"
                                  DefaultAssay(APCs2) <- "integrated"
                                  
                                  
                                  APCs2 <- RunPCA(APCs2, verbose = FALSE, npcs = 75)
                                  
                                  APCs2 <- FindNeighbors(APCs2, dims = 1:75, k.param = 30)
                                  APCs2 <- FindClusters(APCs2, resolution = 2.0)
                                  APCs2 <- RunUMAP(APCs2, dims = 1:30, n.neighbors = 30L)
                                  
                                  DimPlot(APCs2, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()
                                
                                  
                                  plot_tsne_metadata_srt(APCs2, color_by = "seurat_clusters", facet_by = "seurat_clusters", size = 0.75, plot_labels = TRUE, plot_legend = FALSE)
                                  
                                  plot_tsne_gene_srt(APCs2, gene = "LYZ")
                                  plot_violin_srt(APCs2, gene = "CD207", color_by = "seurat_clusters")
                                  
                                  #View integration:
                                  plot_tsne_metadata_srt(APCs2, color_by = "seurat_clusters", facet_by = c("Sequence_instrument"), size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                                  plot_tsne_metadata_srt(APCs2, color_by = "seurat_clusters", facet_by = c("Lesion", "Sequence_instrument"),size = 0.25, plot_labels = FALSE, plot_legend = FALSE, facet_grid = TRUE)
                                  plot_tsne_metadata_srt(APCs2, color_by = "seurat_clusters", facet_by = c("Lesion"),size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
                                  
                                  #Markers:
                                  DefaultAssay(APCs2) <- "RNA"
                                  
                                  load("gene_lists/established_exclusion_genes_5.17.22.Rdata")
                                  
                                  APC2.markers <- FindAllMarkers(APCs2, min.pct = 0.25, logfc.threshold = 0.25)
                                  top10 <- APC2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                                  
                                  Filtered.APC2.markers <- APC2.markers %>%  filter(!gene %in% established_exclusion_genes) %>% filter(avg_log2FC > 50 | avg_log2FC > 2 & pct.2 < 0.075)
                                  APC2_top10_FC_pct2 <- Filtered.APC2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC/(pct.2 + 0.00000001))
                                  
                                  plot_heatmap_srt(APCs2, 
                                                   genes = rev(c("CD207", "FCER1A", "S100B", "CD1B", "LYZ", "PLAUR", "CD1C", "PLEK", "ALOX15")),
                                                   facet_by = "seurat_clusters",
                                                   type = "single_cell",
                                                   col_names = FALSE,
                                                   gene_names = TRUE,
                                                   cluster_by = FALSE,
                                                   ceiling = 2,
                                                   color_pal = colorRampPalette(c("#E5E5E5","#E5E5E5", "#E5E5E5","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256))
                                  
                                  
                                  #Rename Clusters:  
                                  new.cluster.ids <- c("LC", "LC", 
                                                        "Myeloid","Myeloid","Myeloid","Myeloid",
                                                       "LC", "LC","LC","LC","LC","LC")
                                  names(new.cluster.ids) <- levels(APCs2)
                                  APCs2 <- RenameIdents(APCs2, new.cluster.ids)
                                  DimPlot(APCs2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()    
                                  
                                  APCs2@meta.data$Detailed_clusters <- as.character(APCs2@active.ident)

                                  Idents(APCs2) <- "Detailed_clusters"
                                  
                                  APCs@meta.data$Detailed_clusters[which(rownames(APCs@meta.data) %in% rownames(APCs2@meta.data))] <- APCs2@meta.data$Detailed_clusters
                                  APCs@meta.data$Detailed_clusters[which(APCs@meta.data$Detailed_clusters == "Myeloid")] <- "Myeloid_1"
                                  
                                  plot_tsne_metadata_srt(APCs, color_by = "Detailed_clusters", size = 0.75, plot_labels = TRUE, plot_legend = FALSE)
                                  
                                  #save(APCs, file = "Data/APC_subset.Rdata")
                                  #write.table(APC.markers, file = "Data/APC.markers.tsv", row.names=FALSE, sep="\t")
                                  
                                  #save(APCs2, file = "Data/APC2_subset.Rdata")
                                  #write.table(APC2.markers, file = "Data/APC2.markers.tsv", row.names=FALSE, sep="\t")
                                  
            #Recluster mixed KRT populations:
                                  Idents(contact_derm) <- "Detailed_clusters"
                                  
                                  KRT_recluster <- subset(contact_derm, idents = "KRT-b3")
                                  
                                  Idents(KRT_recluster) <- "seurat_clusters"
                                  DefaultAssay(KRT_recluster) <- "integrated"
                                  
                                  
                                  KRT_recluster <- RunPCA(KRT_recluster, verbose = FALSE, npcs = 75)
                                  
                                  KRT_recluster <- FindNeighbors(KRT_recluster, dims = 1:75, k.param = 30)
                                  KRT_recluster <- FindClusters(KRT_recluster, resolution = 0.5)
                                  KRT_recluster <- RunUMAP(KRT_recluster, dims = 1:30, n.neighbors = 30L)
                                  
                                  DimPlot(KRT_recluster, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()
                                  
                                  
                                  plot_tsne_metadata_srt(KRT_recluster, color_by = "seurat_clusters", facet_by = "seurat_clusters", size = 0.75, plot_labels = TRUE, plot_legend = FALSE)
                                  
                                  plot_tsne_gene_srt(KRT_recluster, gene = "KRT77")
                                  plot_violin_srt(KRT_recluster, gene = "KRT6B", color_by = "seurat_clusters")
                                  
                                  #Markers:
                                  DefaultAssay(KRT_recluster) <- "RNA"
                                  
                                  load("gene_lists/established_exclusion_genes_5.17.22.Rdata")
                                  
                                  KRT_recluster.markers <- FindAllMarkers(KRT_recluster, min.pct = 0.25, logfc.threshold = 0.25)
                                  top10 <- KRT_recluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                                  
                                  Filtered.KRT_recluster.markers <- KRT_recluster.markers %>%  filter(!gene %in% established_exclusion_genes) %>% filter(avg_log2FC > 50 | avg_log2FC > 2 & pct.2 < 0.075)
                                  KRT_top10_FC_pct2 <- KRT_recluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC/(pct.2 + 0.00000001))

                                  plot_heatmap_srt(KRT_recluster, 
                                                   genes = rev(c("KRT77", "KRT6A","KRT6B","KRT16","KRT17", "S100A7","S100A8", "S100A9")),
                                                   facet_by = "seurat_clusters",
                                                   type = "single_cell",
                                                   col_names = FALSE,
                                                   gene_names = TRUE,
                                                   cluster_by = FALSE,
                                                   ceiling = 2,
                                                   color_pal = colorRampPalette(c("#E5E5E5","#E5E5E5", "#E5E5E5","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256))
                                  
                                  
                                  
                                  plot_heatmap_srt(KRT_recluster, 
                                                   genes = rev(D_markers_for_plot[1:21]),
                                                   facet_by = "Detailed_clusters",
                                                   type = "single_cell",
                                                   col_names = FALSE,
                                                   gene_names = TRUE,
                                                   cluster_by = FALSE,
                                                   ceiling = 2,
                                                   color_pal = colorRampPalette(c("#E5E5E5","#E5E5E5", "#E5E5E5","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256))
                                  
                                  
                                  #Rename Clusters:  
                                  new.cluster.ids <- c("KRT-b3", "KRT-wr", "KRT-wr", "KRT-wr",
                                                       "KRT-b3", "KRT-b3","KRT-b3")
                                  names(new.cluster.ids) <- levels(KRT_recluster)
                                  KRT_recluster <- RenameIdents(KRT_recluster, new.cluster.ids)
                                  DimPlot(KRT_recluster, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()    
                                  
                                  KRT_recluster@meta.data$Detailed_clusters <- as.character(KRT_recluster@active.ident)
                                  
                                  Idents(KRT_recluster) <- "seurat_clusters"
                                  
                                  plot_tsne_metadata_srt(KRT_recluster, color_by = "Detailed_clusters", size = 0.75, plot_labels = TRUE, plot_legend = FALSE)
                                  
                                  #save(KRT_recluster, file = "Data/KRT_mini_subset.Rdata")
                                  
                                  
        #Detailed_clusters assignment:
        contact_derm@meta.data$Detailed_clusters[which(rownames(contact_derm@meta.data) %in% rownames(Lymphs@meta.data))] <- Lymphs@meta.data$Detailed_clusters
        contact_derm@meta.data$Detailed_clusters[which(rownames(contact_derm@meta.data) %in% rownames(KRTs@meta.data))] <- KRTs@meta.data$Detailed_clusters
        
        contact_derm@meta.data$Detailed_clusters[which(rownames(contact_derm@meta.data) %in% rownames(APCs@meta.data))] <- APCs@meta.data$Detailed_clusters
        contact_derm@meta.data$Detailed_clusters[which(rownames(contact_derm@meta.data) %in% rownames(APCs2@meta.data))] <- APCs2@meta.data$Detailed_clusters
        
        contact_derm@meta.data$Detailed_clusters[which(contact_derm@meta.data$Detailed_clusters == "Myeloid")] <- "Myeloid_1"
        
        contact_derm@meta.data$Detailed_clusters[which(rownames(contact_derm@meta.data) %in% rownames(KRT_recluster@meta.data))] <- KRT_recluster@meta.data$Detailed_clusters
        
        #Rename clusters to celltype:
        colnames(contact_derm@meta.data)[which(colnames(contact_derm@meta.data) == "Detailed_clusters")] <- "Detailed_Celltype"

        Idents(contact_derm) <- "Detailed_Celltype"
        
        #Semi_Detailed_Celltype assignment:
        contact_derm@meta.data$Semi_Detailed_Celltype <- contact_derm@meta.data$Detailed_Celltype
        
        contact_derm@meta.data$Semi_Detailed_Celltype[which(str_detect(contact_derm@meta.data$Semi_Detailed_Celltype, pattern = "KRT"))] <- "KRT"
        contact_derm@meta.data$Semi_Detailed_Celltype[which(str_detect(contact_derm@meta.data$Semi_Detailed_Celltype, pattern = "Mel"))] <- "MEL"
        contact_derm@meta.data$Semi_Detailed_Celltype[which(str_detect(contact_derm@meta.data$Semi_Detailed_Celltype, pattern = "NK"))] <- "NK"
        contact_derm@meta.data$Semi_Detailed_Celltype[which(str_detect(contact_derm@meta.data$Semi_Detailed_Celltype, pattern = "CD8"))] <- "CD8"
        contact_derm@meta.data$Semi_Detailed_Celltype[which(str_detect(contact_derm@meta.data$Semi_Detailed_Celltype, pattern = "CD4"))] <- "CD4"
        contact_derm@meta.data$Semi_Detailed_Celltype[which(str_detect(contact_derm@meta.data$Semi_Detailed_Celltype, pattern = "Myeloid|M2_mac|cDC1"))] <- "Myeloid"
        
        #Rename KRT-b3:
        contact_derm@meta.data$Detailed_Celltype[which(str_detect(contact_derm@meta.data$Detailed_Celltype, pattern = "KRT-b3"))] <- "KRT-77"
        
        #Rename CD4_trm to marker-based name:
        contact_derm@meta.data$Detailed_Celltype[which(str_detect(contact_derm@meta.data$Detailed_Celltype, pattern = "CD4_trm"))] <- "CD4_crem"

        #Rename NK_res to marker-based name:
        contact_derm@meta.data$Detailed_Celltype[which(str_detect(contact_derm@meta.data$Detailed_Celltype, pattern = "NK_res"))] <- "NK_areg"
        
        #Rename Bcell to pDC:
        contact_derm@meta.data$Detailed_Celltype[which(str_detect(contact_derm@meta.data$Detailed_Celltype, pattern = "Bcell"))] <- "pDC"
        contact_derm@meta.data$Semi_Detailed_Celltype[which(str_detect(contact_derm@meta.data$Semi_Detailed_Celltype, pattern = "Bcell"))] <- "pDC"
        
            plot_violin_srt(contact_derm, gene = "JCHAIN", color_by = "Detailed_Celltype")
            plot_violin_srt(contact_derm, gene = "IGKC", color_by = "Detailed_Celltype")
            plot_violin_srt(contact_derm, gene = "IRF8", color_by = "Detailed_Celltype")
            plot_violin_srt(contact_derm, gene = "PLAC8", color_by = "Detailed_Celltype")
            
            plot_violin_srt(contact_derm, gene = "IL3RA", color_by = "Detailed_Celltype")
            plot_violin_srt(contact_derm, gene = "NRP1", color_by = "Detailed_Celltype")
            plot_violin_srt(contact_derm, gene = "CLEC4C", color_by = "Detailed_Celltype")
            plot_violin_srt(contact_derm, gene = "MS4A1", color_by = "Detailed_Celltype")
            plot_violin_srt(contact_derm, gene = "CD19", color_by = "Detailed_Celltype")
            
        
        #Visualize:
            #Build Custom Color Palette:
            flat.palette <- c("#F44336", "#C2185B", "#BA68C8", "#673AB7", "#7986CB", "#1E88E5", "#26C6DA", "#00897B", "#4CAF50",  "#AFB42B", "#FF8F00", "#FF5722", "#8D6E63", "#78909C")
            GetPalette = colorRampPalette(flat.palette)
            
            ColorCount = as.numeric(length(unique(contact_derm@meta.data$Detailed_Celltype)))
            discrete.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
            plot_tsne_metadata_srt(contact_derm, color_by = "Detailed_Celltype", size = 0.25, plot_labels = TRUE, plot_legend = TRUE, colors = discrete.colors1)
            
            plot_tsne_metadata_srt(contact_derm, color_by = "Semi_Detailed_Celltype", size = 0.25, plot_labels = TRUE, plot_legend = TRUE)
            plot_tsne_metadata_srt(contact_derm, color_by = "Basic_Celltype", size = 0.25, plot_labels = TRUE, plot_legend = TRUE)
            
            
    #Convert metagenes to most popular/familiar gene name:
        meta_backto_symbols <- fread("gene_lists/meta_backto_symbols.txt")
        
        contact_derm <- RenameGenesSeurat_frisoli(contact_derm, meta_backto_symbols)
    
        VlnPlot(contact_derm, features = "FOXP3")
        FeaturePlot(contact_derm, features = "FOXP3")
        
    
    #Add Lesion.score to metadata:
        library(lubridate)
        library(readxl)
        lesion.scores <- read_xlsx(path = "Data/master_SADBE_sample_info.xlsx")
        
        lesion.scores$Blister_Date <- as_date(lesion.scores$Blister_Date)
        
        lesion.scores <- lesion.scores %>% filter(Sequenced == "yes")
        lesion.scores <- lesion.scores[,-which(colnames(lesion.scores) %in% c("Paired Microneedle Samples?", "Patient Note", "Freezer Location", "Bioanalyzer", "Lesion Photo", "R1 Sequence Index", "DolphinNext Report", "Sample Cell Concentration per mL", "Total cells", "Sequenced"))]
        
        
        #Format before merge:
            colnames(lesion.scores)[which(colnames(lesion.scores) == "Sequence Date")] <- "Sequence_date"
            colnames(lesion.scores)[which(colnames(lesion.scores) == "Blister_Date")] <- "Blister_date"
            colnames(lesion.scores)[which(colnames(lesion.scores) == "Skin Sample Site location")] <- "Anatomic_site"
            colnames(lesion.scores)[which(colnames(lesion.scores) == "Lesion Visual Score")] <- "Lesion_visual_score"
            colnames(lesion.scores)[which(colnames(lesion.scores) == "SADBE Sensitization Volume (drops)")] <- "SADBE_sensitization_drop_volume"
            colnames(lesion.scores)[which(colnames(lesion.scores) == "SADBE Elicitation Dose (%wt/%wt)")] <- "SADBE_elicitation_dose_percentage"
            
            lesion.scores[which(lesion.scores$Lesion == "NL"), "Lesion"] <- "Nonlesional"
            lesion.scores[which(lesion.scores$Lesion == "Day 2 SLS"), "Lesion"] <- "Irritant"
            lesion.scores[which(lesion.scores$Lesion == "Day 2 Acetone"), "Lesion"] <- "Acetone_Vehicle"
            lesion.scores[which(lesion.scores$Lesion == "Day 2 SADBE"), "Lesion"] <- "Day2_Allergy"
            lesion.scores[which(lesion.scores$Lesion == "Day 4 SADBE"), "Lesion"] <- "Day4_Allergy"
            
            lesion.scores$Lesion <- factor(lesion.scores$Lesion,
                                           levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
            
            #Merge:
            metadata <- contact_derm@meta.data
            metadata <- rownames_to_column(metadata)
        
            if (any(colnames(metadata) %in% "Lesion_visual_score")) {
              metadata <- metadata[,-which(colnames(metadata) %in% "Lesion_visual_score")]
            }
            
            test.meta <- full_join(metadata, lesion.scores)
            test.meta <- test.meta[-which(is.na(test.meta$rowname)),]
        
            test.meta <- column_to_rownames(test.meta, var = "rowname")
            
            contact_derm@meta.data <- test.meta

          #Separate Reactive from Non-reactive Acetone Vehicle:
            contact_derm@meta.data$Lesion_type <- as.character(contact_derm@meta.data$Lesion)
            
            contact_derm@meta.data$Lesion_type[which(contact_derm@meta.data$Lesion_type == "Acetone_Vehicle")] <- "Acetone_Vehicle_Nonreactive"
            
            reactive_acetone_cellIDs <- contact_derm@meta.data %>% filter(Patient %in% c("HC8", "HC12") & 
                                               Lesion == "Acetone_Vehicle" & 
                                               Blister_date %in% as_date(c("2021-11-12", "2021-09-24"))) %>% rownames()
            
            contact_derm@meta.data[which(rownames(contact_derm@meta.data) %in% reactive_acetone_cellIDs),"Lesion_type"] <- "Acetone_Vehicle_Reactive"
          
            contact_derm@meta.data$Lesion_type <- factor(contact_derm@meta.data$Lesion_type,
                                           levels = c("Nonlesional", "Irritant", "Acetone_Vehicle_Nonreactive", "Acetone_Vehicle_Reactive","Day2_Allergy", "Day4_Allergy"))

            contact_derm@meta.data <- contact_derm@meta.data %>% mutate(sex = case_when(Patient %in% c("HC1", "HC8", "HC12","H13", "HC15", "HC16", "HC23",    "CB19","CB22", "CB27", "CB31", "CB48") ~ "Male",
                                                              Patient %in% c("HC3", "HC7", "HC11", "HC14", "HC17", "HC18", "HC20", "HC21", "HC22",    "CB17", "CB20","CB21", "CB23","CB43", "CB44", "CB46") ~ "Female"))
            
        ### Final edits:
            # Change visual score to match standardized patch test scoring system:
            table(contact_derm@meta.data$Lesion_visual_score)
            contact_derm@meta.data <- contact_derm@meta.data %>% mutate(Lesion_severity_score = case_when(Lesion_visual_score == 2 ~ "++",
                                                                                       Lesion_visual_score %in% c(0.5, 1) ~ "+",
                                                                                       Lesion_visual_score %in% c(0.25) ~ "?/+",
                                                                                       Lesion_visual_score == 0 ~ "-"))
            contact_derm@meta.data$Lesion_severity_score <- factor(contact_derm@meta.data$Lesion_severity_score,
                                                                levels = c("-", "?/+", "+", "++"))
            
            # Rename cell types in seurat object to match what I want final names to be:
            contact_derm@meta.data$Detailed_Celltype <- str_replace_all(contact_derm@meta.data$Detailed_Celltype, pattern = "_", replacement = "-")
            
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "KRT-77")] <- "KRT-krt77"
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "KRT-mucl")] <- "KRT-mucl1"
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "Mel1")] <- "MEL-1"
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "Mel2")] <- "MEL-2"
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "M2-mac")] <- "MAC-fcn1"
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "CD4-hsp")] <- "CD4-hspa6"
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "CD4-cd161")] <- "CD4-klrb1"
            contact_derm@meta.data$Detailed_Celltype[which(contact_derm@meta.data$Detailed_Celltype == "CD4-tcm")] <- "CD4-il7r"


            #save(contact_derm, file = "Data/contact_derm_filtered3.Rdata")

            #Rename Patient IDs to exclude dropout patients:
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC1")] <- "P1"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC3")] <- "P2"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC7")] <- "P3"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC8")] <- "P4"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC12")] <- "P5"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC14")] <- "P6"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC16")] <- "P7"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC17")] <- "P8"
            contact_derm@meta.data$Patient[which(contact_derm@meta.data$Patient == "HC18")] <- "P9"
            
            #save(contact_derm, file = "Data/contact_derm_filtered_final.Rdata")
            
            
#################################################################
### Cell Count, Proportion Analysis, and cluster marker heatmaps:
#################################################################
library(Seurat)
library(SeuratObject)
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
library(ggnewscale) 

        
    setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")
    
    #Load data:
    load(file = "Data/contact_derm_filtered_final.Rdata")
    
    #Set Celltypes as factors:
          contact_derm@meta.data$Basic_Celltype <- factor(contact_derm@meta.data$Basic_Celltype,
                                                          levels = c("KRT", "MEL", "APC", "LYMPH"))
          
          contact_derm@meta.data$Semi_Detailed_Celltype <- factor(contact_derm@meta.data$Semi_Detailed_Celltype, 
                                                                  levels = c("KRT", "MEL", "LC","Myeloid","pDC", "CD4", "Treg", "CD8", "NK"))
          
          contact_derm@meta.data$Detailed_Celltype <- factor(contact_derm@meta.data$Detailed_Celltype, 
                                                             levels = c("KRT-b1", "KRT-b2", "KRT-sp", "KRT-g", "KRT-krt77","KRT-mucl1","KRT-wr", "MEL-1", "MEL-2", "LC", "Myeloid-1", "Myeloid-ccl22", "MAC-fcn1", "cDC1","pDC", 
                                                                        "CD4-conv1", "CD4-hspa6", "CD4-klrb1", "CD4-il7r", "CD4-crem", "Treg", "CD8-CD4-active","CD8", "NK-ctl", "NK-areg"))
          
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
          
          #load(file = "Data/APC_subset.Rdata")
          #APCs@meta.data$Detailed_clusters <- factor(APCs@meta.data$Detailed_clusters, levels = c("LC", "Myeloid_1", "M2_mac", "Myeloid_ccl22", "cDC1", "Bcell"))
          #plot_tsne_metadata_srt(APCs, color_by = "Detailed_clusters", size = 0.25, plot_labels = FALSE, plot_legend = TRUE, colors = c("#BF581D","#77AF39","#EC9E0D", "#BA27D1", "#FF9C33","#67008B"))
          
    plot_tsne_metadata_srt(contact_derm, color_by = "Basic_Celltype", size = 0.25, plot_labels = FALSE, plot_legend = TRUE)
    
    plot_tsne_metadata_srt(contact_derm, color_by = "Semi_Detailed_Celltype", size = 0.25, plot_labels = FALSE, plot_legend = TRUE)
    
    plot_tsne_metadata_srt(contact_derm, color_by = "Detailed_Celltype", size = 0.25, plot_labels = TRUE, plot_legend = FALSE)
    
    plot_tsne_metadata_srt(DE.input, color_by = "Detailed_Celltype", facet_by = "Lesion", size = 0.25, plot_labels = FALSE, plot_legend = FALSE)
    
    
    #Build Custom Color Palette:
    flat.palette <- c("#F44336", "#C2185B", "#BA68C8", "#673AB7", "#7986CB", "#1E88E5", "#26C6DA", "#00897B", "#4CAF50",  "#AFB42B", "#FF8F00", "#FF5722", "#8D6E63", "#78909C")
    GetPalette = colorRampPalette(flat.palette)
    
    ColorCount = as.numeric(length(unique(DE.input@meta.data$Detailed_Celltype)))
    discrete.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
    
    plot_tsne_metadata_srt(DE.input, color_by = "Detailed_Celltype", size = 0.1, colors = discrete.colors1, plot_legend = TRUE, plot_labels = FALSE)
    
    detailed.umap <- last_plot()
    
    #Build Custom Color Palette:
    flat.palette <- c("#F44336", "#C2185B", "#BA68C8", "#673AB7", "#7986CB", "#1E88E5", "#26C6DA", "#00897B", "#4CAF50",  "#AFB42B", "#FF8F00", "#FF5722", "#8D6E63", "#78909C")
    GetPalette = colorRampPalette(flat.palette)
    
    ColorCount = as.numeric(length(unique(DE.input@meta.data$Patient)))
    discrete.colors2 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
    
    plot_tsne_metadata_srt(DE.input, color_by = "Patient", size = 0.1, colors = discrete.colors2, plot_legend = TRUE, plot_labels = FALSE)
    
    
    load("Data/APC_subset.Rdata")
    APCs@meta.data$Detailed_clusters[which(APCs@meta.data$Detailed_clusters == "Bcell")] <- "pDC"
    APCs@meta.data$Detailed_clusters <- factor(APCs@meta.data$Detailed_clusters,
                                               levels = c("LC", "Myeloid_1", "Myeloid_ccl22", "M2_mac", "cDC1","pDC"))
    unique(DE.input@meta.data$Detailed_Celltype)
    discrete.colors2 <- discrete.colors1[10:15]
    
    plot_tsne_metadata_srt(APCs, color_by = "Detailed_clusters",colors = discrete.colors2, size = 0.25, plot_legend = TRUE, plot_labels = FALSE)
    detailed.umap.apcs <- last_plot()
    
    #setwd("Plots/umap")
    #ggsave(detailed.umap, width = 40, height = 40,units = "cm",filename = "detailed.umap.png", limitsize = FALSE)  
    #ggsave(detailed.umap.apcs, width = 25, height = 25,units = "cm",filename = "detailed.umap.apcs.png", limitsize = FALSE)  
    
    flat.palette <- c("#F44336", "#C2185B", "#BA68C8", "#673AB7", "#7986CB", "#1E88E5", "#26C6DA", "#00897B", "#4CAF50",  "#AFB42B", "#FF8F00", "#FF5722", "#8D6E63", "#78909C")
    GetPalette = colorRampPalette(flat.palette)
    
    ColorCount = as.numeric(length(unique(DE.input@meta.data$Semi_Detailed_Celltype)))
    discrete.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
    
    plot_tsne_metadata_srt(DE.input, color_by = "Semi_Detailed_Celltype", size = 0.1, colors = discrete.colors1, plot_legend = TRUE, plot_labels = FALSE)
    
    
#Doublet Analysis:
    load(file = "Data/unfiltered.seurat.Rdata")
    
    #Doublets filtered by my method:
    master_data@meta.data$Final.cell <- rownames(master_data@meta.data) %in% colnames(contact_derm)

    plot_tsne_metadata_srt(master_data, color_by = "Final.cell", size = 0.5, colors = c("#F65937", "#5A5A5A"), plot_labels = FALSE, plot_legend = TRUE)
    
    #Comparison to DoubletFinder:
    library(DoubletFinder)

    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list <- paramSweep_v3(master_data, PCs = 1:10, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats) 
    
    ggplot(bcmvn, aes(x = pK, y = BCmetric)) +
      geom_point() +
      geom_line(group = 1, color = "gray", linetype = "dashed") +
      theme_classic()
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- master_data@meta.data$seurat_clusters
    
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.08*nrow(master_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset--> I chose 0.08 based on the plot above
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    master_data <- doubletFinder_v3(master_data, PCs = 1:10, pN = 0.25, pK = 0.08, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    master_data <- doubletFinder_v3(master_data, PCs = 1:10, pN = 0.25, pK = 0.08, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = TRUE)
    
    plot_tsne_metadata_srt(master_data, color_by = "DF.classifications_0.25_0.08_4938", size = 0.5, colors = c("#F65937", "#5A5A5A"), plot_labels = FALSE, plot_legend = TRUE)
    
    # Cell Count Analysis Relative to Published Vitiligo Paper:
    
       #By Patient:
              Vit.pub.cellstats <- fread("Data/vit_published_data/vitiligo.pub.cellstats.tsv")
              Vit.pub.cellstats <- Vit.pub.cellstats %>% group_by_at("Patient") %>% summarize(cellcount = sum(cellcount))
              Vit.pub.cellstats$Data_group <- "Published"
              Vit.pub.cellstats$Disease <- NA
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Patient, pattern = "^CB"))] <- "Published Healthy"
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Patient, pattern = "^VB"))] <- "Vitiligo"
              Vit.pub.cellstats$Sequence_date = NA
              
              cell.stats <- contact_derm@meta.data %>% group_by_at(c("Patient", "Sequence_date")) %>% summarize(cellcount = n())
              cell.stats$Data_group <- "Data from previous TRAC"
              cell.stats$Data_group[which(cell.stats$Sequence_date %in% c("12_3_21", "4_28_22"))] <- "New Data"
              
              cell.stats$Disease <- NA
              cell.stats$Disease[which(str_detect(cell.stats$Patient, pattern = "^CB"))] <- "Healthy"
              cell.stats$Disease[which(str_detect(cell.stats$Patient, pattern = "^HC"))] <- "Contact Dermatitis"
              
              cell.stats <- bind_rows(cell.stats, Vit.pub.cellstats)
              cell.stats$Disease <- factor(cell.stats$Disease, 
                                           levels = c("Published Healthy", "Vitiligo", "Healthy", "Contact Dermatitis"))
              
              
              new.patients <- cell.stats %>% filter(Data_group == "New Data") %>% select(Patient) %>% unlist()
              
              Contact.derm.samples <- cell.stats %>% filter(Disease == "Contact Dermatitis") 
              
              Contact.derm.samples$Patient <- factor(Contact.derm.samples$Patient,
                                                     levels = c("HC1", "HC3", "HC7", "HC8","HC12", "HC14", "HC16", "HC17", "HC18"))
              Contact.derm.samples <- arrange(Contact.derm.samples, Sequence_date, Patient)
              Contact.derm.samples$Blister_date = c("10-11-19", "11-22-19", "10-2-20", "11-20-20", "8-5-21", "8-13-21", "9-17-21", "9-24-21", "11-12-21", "2-11-22", "11-19-21")

              cell.stats <- cell.stats %>% group_by_at(c("Patient", "Disease")) %>% summarize(cellcount = sum(cellcount),
                                                                                            Data_group = str_c(Data_group, collapse = "__"))

              cell.stats$Data_group[which(str_detect(cell.stats$Data_group, pattern = "Published"))] <- "Published"
              cell.stats$Data_group[which(cell.stats$Data_group == "New Data")] <- "New patient data"
              cell.stats$Data_group[which(str_detect(cell.stats$Data_group, pattern = "__"))] <- "Repeat patient with new data"
              cell.stats$Data_group <- factor(cell.stats$Data_group,
                                              levels = c("Published", "Data from previous TRAC", "New patient data", "Repeat patient with new data"))
              
              
              ggplot(cell.stats, aes(x = Disease, y = log10(cellcount), color = Data_group)) + 
                geom_jitter(size = 3, width = 0.25) + 
                theme_classic() +
                scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1), expand = c(0,0))+
                theme(text = element_text(size = 20),
                      axis.text.x = element_text(size = 15),
                      axis.text.y = element_text(size = 12),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 12),
                      axis.title = element_text(size = 17)) +
                labs(y = "log 10 Cell Count per Patient", fill = "Data group", color = "Data group", title = "Cell count by Patient") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                scale_color_manual(values = c("black", "grey38", "red", "red"))
              
              ggplot(Contact.derm.samples, aes(x = Blister_date, y = log10(cellcount), color = Patient)) + 
                geom_point(size = 4) + 
                theme_classic() +
                scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1), expand = c(0,0))+
                theme(text = element_text(size = 20),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 10),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(y = "log10 (Cell Count)", color = "Patient", title = "Contact Derm cell counts per blister visit") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
              
              
              ggplot(cell.stats, aes(x = Disease, y = cellcount, color = Data_group, fill = Data_group)) + 
                geom_bar(stat = "summary")  + 
                theme_classic() +
                scale_y_continuous(limits = c(0, 40000), breaks = seq(0, 40000, by = 5000), expand = c(0,0))+
                theme(text = element_text(size = 20),
                      axis.text = element_text(size = 10),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(y = "Total Cell Count", fill = "Data group", color = "Data group", title = "Total Cells by Disease Category") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
              
              
        #By Patient and Lesion:
              Vit.pub.cellstats <- fread("Data/vit_published_data/vitiligo.pub.cellstats.tsv")
              colnames(Vit.pub.cellstats)[2] <- "Lesion"
              Vit.pub.cellstats$Data_group <- "Published"
            
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Patient, pattern = "^CB"))] <- "Published Healthy"
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "NonLesional"))] <- "Vitiligo Nonlesional"
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Disease"))] <- "Vitiligo Lesional"
              
          
              cell.stats <- contact_derm@meta.data %>% group_by_at(c("Patient","Lesion")) %>% summarize(cellcount = n())
              cell.stats[,1:2] <- sapply(cell.stats[,1:2], FUN = as.character)
              cell.stats$Data_group <- "Manuscript Data"
              
              
              cell.stats$Lesion[which(str_detect(cell.stats$Patient, pattern = "^CB"))] <- "Healthy"
              
              cell.stats <- bind_rows(cell.stats, Vit.pub.cellstats)
              cell.stats$Lesion <- str_replace_all(cell.stats$Lesion, pattern = "_", replacement = " ")
              cell.stats$Lesion <- factor(cell.stats$Lesion, 
                                          levels = c("Published Healthy", "Vitiligo Nonlesional", "Vitiligo Lesional", "Healthy", "Nonlesional", "Irritant", "Acetone Vehicle","Day2 Allergy", "Day4 Allergy"))
              
              ggplot(cell.stats, aes(x = Lesion, y = log10(cellcount), color = Data_group)) + 
                geom_jitter(size = 2.5, width = 0.15) + 
                theme_classic() +
                scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1), expand = c(0,0)) +
                theme(text = element_text(size = 20),
                      axis.text = element_text(size = 10),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(y = "Log10 Cell Count", fill = "Data group", color = "Data group", title = "Cell count by Patient and Lesion") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                scale_color_manual(values = c("red", "black"))
              
      
              ggplot(cell.stats, aes(x = Lesion, y = log10(cellcount), color = Data_group)) + 
                geom_jitter(size = 2.5, width = 0.15) + 
                theme_classic() +
                scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1), expand = c(0,0)) +
                theme(text = element_text(size = 20),
                      axis.text = element_text(size = 10),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(y = "Log10 Cell Count", fill = "Data group", color = "Data group", title = "Cell count by Patient and Lesion") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                geom_hline(yintercept = 2, color = "gray38", linetype = "dashed") +
                scale_color_manual(values = c("red", "black"))
              
              
        #By Disease:
              Vit.pub.cellstats <- fread("Data/vit_published_data/vitiligo.pub.cellstats.tsv")
              Vit.pub.cellstats$Disease <- NA
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Patient, pattern = "^CB"))] <- "Published Healthy"
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Patient, pattern = "^VB"))] <- "Published Vitiligo"
              
              Vit.pub.cellstats <- Vit.pub.cellstats %>% group_by_at("Disease") %>% summarize(cellcount = sum(cellcount))
              Vit.pub.cellstats$Data_group <- "Published"
              
              cell.stats <- contact_derm@meta.data %>% group_by_at(c("Lesion", "Disease", "Sequence_date")) %>% summarize(cellcount = n())
              cell.stats[,1:2] <- sapply(cell.stats[,1:2], FUN = as.character)
              
              cell.stats$Disease[which(cell.stats$Disease == "Allergy")] <- "Contact Dermatitis"
              cell.stats$Data_group <- "Data from previous TRAC"
              cell.stats$Data_group[which(cell.stats$Sequence_date %in% c("12_3_21", "4_28_22"))] <- "New Data"
              cell.stats <- cell.stats %>% group_by_at(c("Disease", "Data_group")) %>% summarize(cellcount = sum(cellcount))
              
              cell.stats <- bind_rows(cell.stats, Vit.pub.cellstats)
              cell.stats$Disease <- factor(cell.stats$Disease, 
                                          levels = c("Published Healthy", "Published Vitiligo", "Healthy", "Contact Dermatitis"))
              
              cell.stats$Data_group <- factor(cell.stats$Data_group,
                                              levels = rev(c("Published", "Data from previous TRAC", "New Data")))
              
              
              ggplot(cell.stats, aes(x = Disease, y = cellcount, color = Data_group, fill = Data_group)) + 
                geom_bar(stat = "summary")  + 
                theme_classic() +
                scale_y_continuous(limits = c(0, 40000), breaks = seq(0, 40000, by = 5000), expand = c(0,0))+
                theme(text = element_text(size = 20),
                      axis.text.x = element_text(size = 15),
                      axis.text.y = element_text(size = 12),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 12),
                      axis.title = element_text(size = 17)) +
                labs(y = "Total Cell Count", fill = "Data group", color = "Data group", title = "Total Cells by Disease Category") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                scale_color_manual(values = rev(c("black", "grey38", "red"))) +
                scale_fill_manual(values = rev(c("black", "grey38", "red")))
                
              
        #By Lesion:
              Vit.pub.cellstats <- fread("Data/vit_published_data/vitiligo.pub.cellstats.tsv")
              colnames(Vit.pub.cellstats)[2] <- "Lesion"
              Vit.pub.cellstats <- Vit.pub.cellstats %>% group_by_at("Lesion") %>% summarize(cellcount = sum(cellcount))
              Vit.pub.cellstats$Data_group <- "Published"
              
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Healthy"))] <- "Published Healthy"
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "NonLesional"))] <- "Vitiligo Nonlesional"
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Disease"))] <- "Vitiligo Lesional"
              
              Vit.pub.cellstats$Disease <- NA
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Healthy"))] <- "Published Healthy"
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Vitiligo"))] <- "Published Vitiligo"
              
              
              cell.stats <- contact_derm@meta.data %>% group_by_at(c("Disease", "Lesion")) %>% summarize(cellcount = n())
              cell.stats[,1:2] <- sapply(cell.stats[,1:2], FUN = as.character)
              cell.stats$Lesion[which(str_detect(cell.stats$Disease, pattern = "Healthy"))] <- "Healthy"
              cell.stats$Data_group <- "Manuscript Data"
              
              cell.stats <- bind_rows(cell.stats, Vit.pub.cellstats)
              cell.stats$Lesion <- str_replace_all(cell.stats$Lesion, pattern = "_", replacement = " ")
              cell.stats$Lesion <- factor(cell.stats$Lesion, 
                                          levels = c("Published Healthy", "Vitiligo Nonlesional", "Vitiligo Lesional", "Healthy", "Nonlesional", "Irritant","Acetone Vehicle", "Day2 Allergy", "Day4 Allergy"))
              
              ggplot(cell.stats, aes(x = Lesion, y = log10(cellcount), color = Data_group, fill = Data_group)) + 
                geom_bar(stat = "summary")  + 
                theme_classic() +
                scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1), expand = c(0,0)) +
                theme(text = element_text(size = 20),
                      axis.text = element_text(size = 10),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(y = "Total Log10 Cell Count", fill = "Data group", color = "Data group", title = "Total cell count by lesion") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                scale_fill_manual(values = c("red", "black")) +
                scale_color_manual(values = c("red", "black"))
              
              
        #By Lesion and CellType:
              Vit.pub.cellstats <- fread("Data/vit_published_data/vitiligo.pub.cellstats.bycelltype.tsv")
              colnames(Vit.pub.cellstats)[2] <- "Lesion"
              Vit.pub.cellstats$CellType[which(Vit.pub.cellstats$CellType %in% c("MAC", "DC"))] <- "APC"
              Vit.pub.cellstats$CellType[which(Vit.pub.cellstats$CellType == "TC")] <- "LYMPH"
              
              Vit.pub.cellstats <- Vit.pub.cellstats %>% group_by_at(c("Lesion", "CellType")) %>% summarize(cellcount = sum(cellcount))
              Vit.pub.cellstats$Data_group <- "Published"
              
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Healthy"))] <- "Published Healthy"
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "NonLesional"))] <- "Vitiligo Nonlesional"
              Vit.pub.cellstats$Lesion[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Disease"))] <- "Vitiligo Lesional"
              
              Vit.pub.cellstats$Disease <- NA
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Healthy"))] <- "Published Healthy"
              Vit.pub.cellstats$Disease[which(str_detect(Vit.pub.cellstats$Lesion, pattern = "Vitiligo"))] <- "Published Vitiligo"
              colnames(Vit.pub.cellstats)[2] <- "Basic_Celltype"
              
              cell.stats <- contact_derm@meta.data %>% group_by_at(c("Disease", "Lesion", "Basic_Celltype")) %>% summarize(cellcount = n())
              cell.stats[,1:2] <- sapply(cell.stats[,1:2], FUN = as.character)
              cell.stats$Lesion[which(str_detect(cell.stats$Disease, pattern = "Healthy"))] <- "Healthy"
              cell.stats$Data_group <- "Manuscript Data"
              
              cell.stats <- bind_rows(cell.stats, Vit.pub.cellstats)
              cell.stats$Lesion <- str_replace_all(cell.stats$Lesion, pattern = "_", replacement = " ")
              cell.stats$Lesion <- factor(cell.stats$Lesion, 
                                          levels = c("Published Healthy", "Vitiligo Nonlesional", "Vitiligo Lesional", "Healthy", "Nonlesional", "Irritant","Acetone Vehicle", "Day2 Allergy", "Day4 Allergy"))
            
              ggplot(cell.stats, aes(x = Lesion, y = log10(cellcount), color = Data_group, fill = Data_group)) + 
                geom_bar(stat = "summary")  + 
                theme_classic() +
                facet_wrap(.~Basic_Celltype, scales = "free") +
                scale_y_continuous(expand = c(0,0)) +
                geom_text(aes(y=log10(cellcount) * 1.25, label="")) +
                theme(text = element_text(size = 20),
                      axis.text = element_text(size = 10),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(y = "Total Cell Count", fill = "Data group", color = "Data group", title = "Cell count by Celltype and Lesion") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                scale_fill_manual(values = c("red", "black")) +
                scale_color_manual(values = c("red", "black"))

# Mean UMI count for each patient by basic celltype:
              
    UMI.data_per_patient_celltype <- function(input,celltype, patient = "Patient") {
                
                bulks <- input@meta.data %>% group_by_at(c(celltype,patient)) %>% summarize(cell.count = n()) %>% as.data.frame()
                
                bulks$mean.UMI.sum = NA
                bulks$mean.gene.sum = NA
                
                for (row in 1:nrow(bulks)) {
                  
                  pat = bulks[row, patient]
                  cell = bulks[row, celltype]
                  
                  group.cells <- input@meta.data %>% filter(!!as.symbol(patient) == pat & !!as.symbol(celltype) == cell) %>% rownames()
                  
                  group.gene.data <- input@assays$RNA@counts[,which(colnames(input@assays$RNA@counts) %in% group.cells)]
                  
                  group.gene.data2 <- group.gene.data
                  group.gene.data2[which(group.gene.data2 > 0)] <- 1
                  
                  if (length(group.cells) >1) {
                    gSums <- apply(group.gene.data2, 2, sum)
                    UMI.sum <- apply(group.gene.data, 2, sum)   
                  }
                  
                  if (length(group.cells) == 1) {
                    gSums <- sum(group.gene.data2)
                    UMI.sum <- sum(group.gene.data)   
                  }
                  
                  bulks$mean.UMI.sum[row] <- round(mean(UMI.sum), digits = 2)
                  bulks$mean.gene.sum[row] <- round(mean(gSums), digits = 2)
                  
                }
                return(bulks)
              }
              
      Vit.pub.umi.data <- fread("Data/vit_published_data/vit_umi.cell.data.tsv")
      Vit.pub.umi.data$Disease <- NA
      Vit.pub.umi.data$Disease[which(str_detect(Vit.pub.umi.data$Patient, pattern = "CB"))] <- "Published Healthy"
      Vit.pub.umi.data$Disease[which(str_detect(Vit.pub.umi.data$Patient, pattern = "VB"))] <- "Vitiligo"
      

      UMI_gene_counts <- UMI.data_per_patient_celltype(contact_derm, celltype = "Basic_Celltype")
         
      UMI_gene_counts$Disease <- NA
      UMI_gene_counts$Disease[which(str_detect(UMI_gene_counts$Patient, pattern = "CB"))] <- "Healthy"
      UMI_gene_counts$Disease[which(str_detect(UMI_gene_counts$Patient, pattern = "HC"))] <- "Contact Dermatitis"
      
      UMI_gene_counts <- bind_rows(UMI_gene_counts, Vit.pub.umi.data)
      
      
      UMI_gene_counts$Basic_Celltype <- factor(UMI_gene_counts$Basic_Celltype,
                                               levels = c("LYMPH", "APC", "MEL", "KRT"))
      UMI_gene_counts$Disease <- factor(UMI_gene_counts$Disease,
                                               levels = c("Published Healthy","Vitiligo","Healthy", "Contact Dermatitis"))
      
        ggplot(UMI_gene_counts, aes(x = Basic_Celltype, y = log10(mean.UMI.sum), color = Basic_Celltype)) + 
          geom_jitter(size = 2, width = 0.25) + 
          facet_grid(~Disease)+
          theme_classic() +
          scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, by = 1), expand = c(0,0)) +
          theme(text = element_text(size = 20),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                plot.title = element_text(hjust = 0.5),
                legend.text=element_text(size=12),
                legend.title = element_text(size = 12),
                axis.title = element_text(size = 17)) +
          labs(y = "log10 mean cell UMI sum per Patient", color = "Basic Celltype", title = "UMI averages") +
          scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

           
              
# Cell Cluster Marker Analysis:
     
          #Manual Selected Markers:
                MEL_specific_markers <- c("TYRP1", "MLANA", "PMEL", "MITF", "DCT", "KIT")
                KRT_specific_markers <- c("KRT10", "KRT14", "KRT17", "DMKN", "KRTDAP", "MUCL1")
                LYMPH_specific_markers <- c("TRAC", "TRDC", "NCAM1")
                APC_specific_markers <- c("CD207", "LYZ", "HLA.DR", "CD86", "CD83", "CLEC4A", "CLEC4C", "MRC1", "NRP1", "JCHAIN")
                RBC_specific_markers <- c("HBB", "HBA2")
                
                Basic_cell_specific_markers <- c(MEL_specific_markers, KRT_specific_markers, LYMPH_specific_markers, APC_specific_markers, RBC_specific_markers)
                additional.markers.from.vit.pub = c("CFTR", "KRT2", "KRT1", "KRT5", "CD207", "PLAUR", "TYR")
                more.lymph.markers = c("CD4", "CD8A", "FOXP3", "IL2RA", "JCHAIN", "IGKC")
      
                plot_heatmap_srt(contact_derm,
                                 genes = unique(c(Basic_cell_specific_markers, additional.markers.from.vit.pub, more.lymph.markers)),
                                 facet_by = "Semi_Detailed_Celltype",
                                 type = "single_cell",
                                 group_names = FALSE,
                                 cluster_by = "both",
                                 ceiling = 2,
                                 color_pal = colorRampPalette(c("#1f1137","#000004FF", "#000004FF","#ff6666","#FF7B25","#FCFDBFFF"))(256))
          
                color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
                  scale = (length(lut)-1)/(max-min)
                  
                  dev.new(width=1.75, height=5)
                  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
                  axis(2, ticks, las=1)
                  for (i in 1:(length(lut)-1)) {
                    y = (i-1)/scale + min
                    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
                  }
                }
                
                color.bar(colorRampPalette(c("#5d00ba", "#000004FF","#000004FF","#ff6666","#FCFDBFFF"))(256), -1)
          
              
                
                
    #Recalculate markers:
          downsampled.CD <- contact_derm[,sample(colnames(contact_derm), size = 20000, replace = FALSE)]
          
          plot_tsne_metadata_srt(downsampled.CD, color_by = "Detailed_Celltype", size = 0.2)
          
          load("gene_lists/established_exclusion_genes_5.17.22.Rdata")
          
              #somehow 3 ribosomal, mitochondrial genes didn't get excluded:
              mito.rpl.genes <- rownames(contact_derm)[which(str_detect(rownames(contact_derm), pattern = "^RPL|^RPS|MT-"))]
              established_exclusion_genes <- unique(c(established_exclusion_genes, mito.rpl.genes[-which(mito.rpl.genes %in% established_exclusion_genes)]))
              
          load("Data/final_celltype_markers/basic_celltype_markers.Rdata")
          load("Data/final_celltype_markers/semi_detailed_celltype_markers.Rdata")
          load("Data/final_celltype_markers/detailed_celltype_markers.Rdata")
          
          #Detailed Markers:
          #Detailed.markers <- FindAllMarkers(downsampled.CD, min.pct = 0.25, logfc.threshold = 0.25)
          Filtered.Detailed.markers <- Detailed.markers %>% filter(!gene %in% established_exclusion_genes) %>% filter(avg_log2FC > 50 | avg_log2FC > 2 & pct.2 < 0.075)
          D_top10_FC_pct2 <- Filtered.Detailed.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC/(pct.2 + 0.00000001)) %>% ungroup() %>% select(gene) %>% unlist() %>% as.character() %>% unique()
          
                #Make sure some of my manually selected markers are included:
                old.D.plot.genes <- rev(c("KRT15", "KRT14", "KRT5", "DMKN", "KRT10", "KRT1", "KRTDAP", "S100A14", "DSC1", "KRT2", "KRT6A", "KRT17", "S100A9", "S100A8", "S100A7", "S100P", "MUCL1", "CEACAM6",
                                          "SCGB1D2", "LSAMP", "AKAP12","SRPX", "KIT", "GPM6B", "CAPN3", "MLANA", "PMEL", "EDNRB", "MFSD12", "APOD", "FCER1A", "CD207", "S100B", "F13A1", "ALOX15", "CD1B", "MMP12",
                                          "PLEK", "CD1C", "LYZ", "IL1B", "PLAUR", "NEAT1", "NRP2", "GPR157", "IDO1", "PLAAT3", "CCL22", "TRAF1", "TCF4", "TSPAN13", "IGKC", "JCHAIN", "IRF8", "CD3G", "SERPINH1",
                                          "HSPA6", "LTB", "KLRB1", "IL7R", "KLF2", "CXCR4", "CREM", "FOXP3", "CTLA4", "IL2RA", "TNFRSF4", "MIR155HG", "MCM6", "GINS2", "GZMK", "CCL5", "CD8A", "NKG7", "GZMA", "KLRD1",
                                          "GNLY", "XCL1", "XCL2", "CTSW", "SPINK2"))
                
                D.genes.to.plot <- unique(c(D_top10_FC_pct2, old.D.plot.genes))
                
          #Semi-Detailed Markers:
          #Idents(downsampled.CD) <- "Semi_Detailed_Celltype"
          #Semi_Detailed.markers <- FindAllMarkers(downsampled.CD, min.pct = 0.25, logfc.threshold = 0.25)
          Filtered.Semi_Detailed.markers <- Semi_Detailed.markers %>% filter(!gene %in% established_exclusion_genes) %>% filter(avg_log2FC > 50 | avg_log2FC > 2 & pct.2 < 0.075)
          SD_top10_FC_pct2 <- Filtered.Semi_Detailed.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC/(pct.2 + 0.00000001)) %>% ungroup() %>% select(gene) %>% unlist() %>% as.character() %>% unique()
          
          old.SD.plot.genes <- rev(c("S100A16", "S100A14", "DMKN", "KRT10", "KRT1", "KRTDAP", "KRT5", "KRT14", "KRT15",
                                     "PMEL", "TYRP1", "EDNRB", "TYR", "DCT", "PCSK2", "GPM6B", "KIT", "OCA2", "LSAMP",
                                     "FCER1A", "CD207", "S100B", "CD1B", "MMP12", "PLEK", "CD1C", "ALOX15", "LYZ", "PLAUR", 
                                     "NEAT1", "CD83", "IL1B", "TCF4", "TSPAN13", "IGKC", "JCHAIN", "IRF8", "CD2", "CD3D", "CD3G", 
                                     "LAT", "TRBC1", "CXCR4", "BCL11B", "FOXP3", "TNFRSF4", "CTLA4", "CCND2", "IL2RA", "GZMK", "CCL5",
                                     "NKG7", "CD8B", "CD8A", "AREG", "XCL2", "XCL1", "CTSW", "GNLY"))
          
          SD.genes.to.plot <- unique(c(SD_top10_FC_pct2, old.SD.plot.genes))
          
          #Basic Markers:
          #Idents(downsampled.CD) <- "Basic_Celltype"

          #Basic.markers <- FindAllMarkers(downsampled.CD, min.pct = 0.25, logfc.threshold = 0.25)
          Filtered.Basic.markers <- Basic.markers %>% filter(!gene %in% established_exclusion_genes) %>% filter(avg_log2FC > 50 | avg_log2FC > 2 & pct.2 < 0.075)
          Basic_top10_FC_pct2 <- Filtered.Basic.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC/(pct.2 + 0.00000001)) %>% ungroup() %>% select(gene) %>% unlist() %>% as.character() %>% unique()
          
          #save(Basic.markers, file = "Data/final_celltype_markers/basic_celltype_markers.Rdata")
          #save(Semi_Detailed.markers, file = "Data/final_celltype_markers/semi_detailed_celltype_markers.Rdata")
          #save(Detailed.markers, file = "Data/final_celltype_markers/detailed_celltype_markers.Rdata")
          
          
          #Basic Marker Heatmap:
              set.seed(200)
              Basic_markers <- plot_heatmap_srt(contact_derm,
                                             genes = Basic_top10_FC_pct2,
                                             facet_by = "Basic_Celltype",
                                             type = "single_cell",
                                             col_names = FALSE,
                                             gene_names = FALSE,
                                             cluster_by = "row",
                                             ceiling = 2,
                                             color_pal = colorRampPalette(c("#97A0CD","#E5E5E5", "#E5E5E5","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256),
                                             cluster_type = "kmeans",
                                             k = 4,
                                             show_k = TRUE)
              
              int_Basic_markers <- c(3,1,2,4)
              int_Basic_markers <- rev(int_Basic_markers)
              Basic_markers_reorder <- c()
              
              for (i in 1:length(int_Basic_markers)) {
                int1 <- int_Basic_markers[i]
                ind <- grep(paste0("^", int1, "$"), Basic_markers[[2]]$cluster)
                reorder <- names(Basic_markers[[2]]$cluster)[ind]
                Basic_markers_reorder <- c(Basic_markers_reorder, reorder)
              }
              
              #remove markers that don't look pure:
              test <- Filtered.Basic.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC/(pct.2 + 0.00000001))
              Basic_markers_reorder <- Basic_markers_reorder[-which(Basic_markers_reorder %in% c("ITGB2", "FAM107B", "SLC2A3", "S100A2"))]
              
             plot_heatmap_srt(contact_derm, 
                               genes = Basic_markers_reorder,
                               facet_by = "Basic_Celltype",
                               type = "single_cell",
                               col_names = FALSE,
                               gene_names = TRUE,
                               cluster_by = FALSE,
                               ceiling = 2,
                               color_pal = colorRampPalette(c("#F7FAEE","#F7FAEE", "#F7FAEE","#F7FAEE","#F7FAEE", "#FC8FE1","#B000A8"))(256))
              
             basic_marker.heatmap <- last_plot()

          #Semi_Detailed Heatmap:
                set.seed(200)
                SD_markers <- plot_heatmap_srt(contact_derm,
                                               genes = SD.genes.to.plot,
                                               facet_by = "Semi_Detailed_Celltype",
                                               type = "single_cell",
                                               col_names = FALSE,
                                               gene_names = FALSE,
                                               cluster_by = "row",
                                               ceiling = 2,
                                               color_pal = colorRampPalette(c("#97A0CD","#E5E5E5", "#E5E5E5","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256),
                                               cluster_type = "kmeans",
                                               k = 12,
                                               show_k = TRUE)

                int_SD_markers <- c(11,8,4,   12,3,9,1,   6,5,7,10,2)
                int_SD_markers <- rev(int_SD_markers)
                SD_markers_reorder <- c()
                
                for (i in 1:length(int_SD_markers)) {
                  int1 <- int_SD_markers[i]
                  ind <- grep(paste0("^", int1, "$"), SD_markers[[2]]$cluster)
                  reorder <- names(SD_markers[[2]]$cluster)[ind]
                  SD_markers_reorder <- c(SD_markers_reorder, reorder)
                }
                
                
                #remove markers that don't look pure:
                test <- Semi_Detailed.markers %>% filter(gene %in% SD_markers_reorder)
                SD_markers_reorder <- SD_markers_reorder[-which(SD_markers_reorder %in% c("DNAJA4", "MPEG1", "UGCG", "CD48", "FYB1", "IL2RG", "FAM107B", 
                                                                                          "TAPBPL", "LTB", "WIPF1", "TNFRSF1B", "BATF", "TNFRSF18", "IL2RB", 
                                                                                          "APOBEC3G", "ITM2A", "GZMB", "CD7", "CD96", "CD68", "SLC2A3", "S100A2", "GLRX",
                                                                                          "APOBEC3C", "CLEC4A", "S100B", "LGMN", "CORO1B", "IKZF1", "BCL11B", "CDC42SE2", "XIST",
                                                                                          "TSPAN13", "CD74", "CD1C", "FNBP1", "FRZB", "LSAMP", "FN1", "SLC6A17", "TCF4", "CD2", "CD3D",
                                                                                          "PCSK2", "GPM6B", "CXCR4", "MMP12", "NEAT1", "ALOX15", "SERPINF1","SPOCK2", "LAT", "CSF2RB"))]
              
                
                #Adjust order:
                SD.genes.to.plot <- rev(SD_markers_reorder)
                    
                    SD.genes.to.plot <- SD.genes.to.plot[-which(SD.genes.to.plot %in% c("CD207"))]
                    SD.genes.to.plot <- c(SD.genes.to.plot[1:(which(SD.genes.to.plot == "FCER1A")-1)], c("CD207"), SD.genes.to.plot[(which(SD.genes.to.plot == "FCER1A"):length(SD.genes.to.plot))])
                    
                    
                    SD.genes.to.plot <- SD.genes.to.plot[-which(SD.genes.to.plot %in% c("NKG7"))]
                    SD.genes.to.plot <- c(SD.genes.to.plot[1:(which(SD.genes.to.plot == "XCL2")-1)], c("NKG7"), SD.genes.to.plot[(which(SD.genes.to.plot == "XCL2"):length(SD.genes.to.plot))])
                    
                    
                    SD.genes.to.plot <- SD.genes.to.plot[-which(SD.genes.to.plot %in% c("CTSW"))]
                    SD.genes.to.plot <- c(SD.genes.to.plot[1:(which(SD.genes.to.plot == "XCL2")-1)], c("CTSW"), SD.genes.to.plot[(which(SD.genes.to.plot == "XCL2"):length(SD.genes.to.plot))])
                    
                    SD.genes.to.plot.withCD4 <- c(SD.genes.to.plot[1:(which(SD.genes.to.plot == "FOXP3")-1)], c("CD4"), SD.genes.to.plot[(which(SD.genes.to.plot == "FOXP3"):length(SD.genes.to.plot))])
                    
                    
                    
                SD_markers_reorder <- rev(SD.genes.to.plot)
                
                plot_heatmap_srt(contact_derm, 
                                 genes = SD_markers_reorder,
                                 facet_by = "Semi_Detailed_Celltype",
                                 type = "single_cell",
                                 col_names = FALSE,
                                 gene_names = TRUE,
                                 cluster_by = FALSE,
                                 ceiling = 2,
                                 color_pal = colorRampPalette(c("#F7FAEE","#F7FAEE", "#F7FAEE","#F7FAEE","#F7FAEE", "#FC8FE1","#B000A8"))(256))
                
                SD_marker.heatmap <- last_plot()
                
                
         #Detailed Heatmap:
                  set.seed(100)
                  D_markers <- plot_heatmap_srt(contact_derm,
                                                 genes = D.genes.to.plot,
                                                 facet_by = "Detailed_Celltype",
                                                 type = "single_cell",
                                                 col_names = FALSE,
                                                 gene_names = FALSE,
                                                 cluster_by = "row",
                                                 ceiling = 2,
                                                 color_pal = colorRampPalette(c("#97A0CD","#E5E5E5", "#E5E5E5","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256),
                                                 cluster_type = "kmeans",
                                                 k = 25,
                                                 show_k = TRUE)

                  int_D_markers <- c(23,9,18,25,21,16,
                                     15,13,
                                     12,8,3,1,4,6,
                                     2,17,22,10,24,19,5,7,11,20,14)
            
                  int_D_markers <- rev(int_D_markers)
                  D_markers_reorder <- c()
                  
                  for (i in 1:length(int_D_markers)) {
                    int1 <- int_D_markers[i]
                    ind <- grep(paste0("^", int1, "$"), D_markers[[2]]$cluster)
                    reorder <- names(D_markers[[2]]$cluster)[ind]
                    D_markers_reorder <- c(D_markers_reorder, reorder)
                  }
                  
                  #remove markers that don't look pure:
                  test <- Detailed.markers %>% filter(gene %in% D_markers_reorder)
                  D_markers_reorder <- D_markers_reorder[-which(D_markers_reorder %in% c("KRT6B", "AL137127.1", "SDC1", "SERPINA12", "SULT2B1","ADIRF","MAL","CD68","DNAJA4","MPEG1", "UGCG",
                                                                                         "LAT","FYB1","CD3G", "CD3D", "CD2", "SPOCK2", "TRBC1", "UGCG","CD40LG", "TAPBPL", "FLT3LG", "FNBP1",
                                                                                         "BATF", "SUN2", "TNFRSF18","TMUB1", "CD69", "SLC2A3","CREM", "FOSL2","TNFRSF1B", "CDC45","SH2D2A","IL2RB",
                                                                                         "SHMT2","ITGB7", "SLURP1", "S100A2","DUSP5", "MAP4K4","CCR7","GADD45A","WIPF1","CD7", "CD96","TNFAIP3",
                                                                                         "REL","DUSP4","LTA","CNTRL","SLC29A1","CD8B", "AREG","MCOLN3", "TMEM47", "LPL", "GZMB", "NR4A2", "CD3E", "DEGS1", "BBOX1",
                                                                                         "GLRX", "SERPINF1", "S100B", "CD1C", "CCDC88A", "F13A1", "CD40","AOC1", "IDO1","TCF4", "FAU", "NACA", "ZNF395",
                                                                                         "ITM2A", "CDC42SE2", "SAMSN1", "CYLD", "LPXN", "CLEC2D", "APOBEC3G"))]
                
                  
                  #slightly adjust gene order:
                  D_markers_for_plot <- rev(D_markers_reorder)
                        
                        D_markers_for_plot <- c("KRT15", D_markers_for_plot[-which(D_markers_for_plot == "KRT15")])
                  
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("KRT2","KRT77", "MUCL1"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "SCGB1D2")-1)], c("KRT2","KRT77", "MUCL1"), D_markers_for_plot[(which(D_markers_for_plot == "SCGB1D2"):length(D_markers_for_plot))])
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("SERPINH1", "HSPA6"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "LTB")-1)], c("SERPINH1", "HSPA6"), D_markers_for_plot[(which(D_markers_for_plot == "LTB"):length(D_markers_for_plot))])
                  
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("KLRD1", "GNLY"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "XCL1")-1)], c("KLRD1", "GNLY"), D_markers_for_plot[(which(D_markers_for_plot == "XCL1"):length(D_markers_for_plot))])
                        
                        #Add CD3G:
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "SERPINH1")-1)], c("CD3G"), D_markers_for_plot[(which(D_markers_for_plot == "SERPINH1"):length(D_markers_for_plot))])
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("CXCR4"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "FOXP3")-1)], c("CXCR4", "CREM"), D_markers_for_plot[(which(D_markers_for_plot == "FOXP3"):length(D_markers_for_plot))])
                        
                        all.mel.markers <- c("KIT", "GPM6B", "CAPN3", "MLANA", "TYRP1", "DCT", "TYR", "PMEL", "MFSD12", "EDNRB","APOD","SGCD", "LSAMP", "AKAP12", "SRPX", "FN1","PTPRM")
                        new.mel.markers <- c("TYR", "DCT", "KIT", "GPM6B", "CAPN3", "MLANA", "PMEL", "EDNRB", "MFSD12", "APOD")
                          
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% all.mel.markers)]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "FCER1A")-1)], new.mel.markers, D_markers_for_plot[(which(D_markers_for_plot == "FCER1A"):length(D_markers_for_plot))])
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% "CD207")]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "FCER1A")-1)], "CD207", D_markers_for_plot[(which(D_markers_for_plot == "FCER1A"):length(D_markers_for_plot))])
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("CD83"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "PLAAT3")-1)], c("CD83"), D_markers_for_plot[(which(D_markers_for_plot == "PLAAT3"):length(D_markers_for_plot))])
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("IL1B"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "PLEK")-1)], c("IL1B"), D_markers_for_plot[(which(D_markers_for_plot == "PLEK"):length(D_markers_for_plot))])
                        
                        #Add BATF3 and CLEC9A for cDC1 markers:
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "IRF8")-1)], c("BATF3"), D_markers_for_plot[(which(D_markers_for_plot == "IRF8"):length(D_markers_for_plot))])
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("IRF8", "BATF3"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "TSPAN13")-1)], c("BATF3","IRF8", "CLEC9A"), D_markers_for_plot[(which(D_markers_for_plot == "TSPAN13"):length(D_markers_for_plot))])
                        

                        #Add more Mac markers:
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "BATF3")-1)], c("MRC1", "FCN1"), D_markers_for_plot[(which(D_markers_for_plot == "BATF3"):length(D_markers_for_plot))])
                        
                        #Reorder Lymph genes:
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("CD3G"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "BCL11B")-1)], c("CD3G"), D_markers_for_plot[(which(D_markers_for_plot == "BCL11B"):length(D_markers_for_plot))])
                        
                            #Add beta chain:
                            D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "IL7R")-1)], c("TRBC"), D_markers_for_plot[(which(D_markers_for_plot == "IL7R"):length(D_markers_for_plot))])
                         
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("TRAC"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "TRBC")-1)], c("TRAC"), D_markers_for_plot[(which(D_markers_for_plot == "TRBC"):length(D_markers_for_plot))])
                        
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("IL7R"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "KLF2")-1)], c("IL7R"), D_markers_for_plot[(which(D_markers_for_plot == "KLF2"):length(D_markers_for_plot))])
                        
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("XCL1"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "XCL2")-1)], c("XCL1"), D_markers_for_plot[(which(D_markers_for_plot == "XCL2"):length(D_markers_for_plot))])
                        
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("CTSW"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "XCL1")-1)], c("CTSW"), D_markers_for_plot[(which(D_markers_for_plot == "XCL1"):length(D_markers_for_plot))])
                        
                        D_markers_for_plot <- c(D_markers_for_plot, "AREG")
                        
                        D_markers_for_plot <- D_markers_for_plot[-which(D_markers_for_plot %in% c("KRT86"))]
                        D_markers_for_plot <- c(D_markers_for_plot[1:(which(D_markers_for_plot == "AREG")-1)], c("KRT86"), D_markers_for_plot[(which(D_markers_for_plot == "AREG"))])
                  
                  D_markers_reorder <- rev(D_markers_for_plot)
                  
                  
                  plot_heatmap_srt(DE.input, 
                                   genes = D_markers_reorder,
                                   facet_by = "Detailed_Celltype",
                                   type = "single_cell",
                                   col_names = FALSE,
                                   gene_names = TRUE,
                                   cluster_by = FALSE,
                                   ceiling = 2,
                                   color_pal = colorRampPalette(c("#F7FAEE","#F7FAEE", "#F7FAEE","#F7FAEE","#F7FAEE", "#FC8FE1","#B000A8"))(256),
                                   scfacet_width_proportioal_to_cellnum = FALSE)
                  
                  
                  D_marker.heatmap <- last_plot()
                  
                  setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/marker_heatmaps")
                  #ggsave(basic_marker.heatmap, width = 25, height = 20,units = "cm",filename = "basic_marker.heatmap.pdf", limitsize = FALSE)  
                  #ggsave(SD_marker.heatmap, width = 25, height = 20,units = "cm",filename = "SD_marker.heatmap.pdf", limitsize = FALSE)  
                  #ggsave(D_marker.heatmap, width = 30, height = 25,units = "cm",filename = "D_marker.heatmap.equalfacetwidth.png", limitsize = FALSE)  
                  #ggsave(D_marker.heatmap, width = 30, height = 40,units = "cm",filename = "D_marker.heatmap.equalfacetwidth.tall.png", limitsize = FALSE)  
                  
                  
                  
                  
                  
                  
                  
                                
      ## Initial analysis by Basic_CellType:
              sample.cell.composition <- contact_derm@meta.data %>% group_by_at(c("Lesion", "Patient", "Basic_Celltype", "Sequence_date", "Lesion_visual_score")) %>%
                                                                  summarise(n()) %>%
                                                                  rename_at("n()", ~ "cellcount") %>% as.data.frame()
              
              sample.cell.composition$Lesion_visual_score[which(is.na(sample.cell.composition$Lesion_visual_score))] <- 0
              
                    ##Only analyze samples with greater than 100 cells:
                    sample.cell.count <- contact_derm@meta.data %>% group_by_at(c("Lesion", "Patient", "Sequence_date")) %>%
                                                                  summarise(n()) %>%
                                                                  rename_at("n()", ~ "cellcount") %>% as.data.frame()
                                  
                    sample.cell.count$inclusion.factor <- NA
                    sample.cell.count <- sample.cell.count %>% mutate(inclusion.factor = ifelse(cellcount >= 100, "yes", "no"))
                    sample.cell.composition <- full_join(sample.cell.composition, sample.cell.count[,-which(colnames(sample.cell.count) == "cellcount")])
                    sample.cell.composition <- sample.cell.composition %>% filter(inclusion.factor == "yes")
                    sample.cell.composition <- sample.cell.composition[,-which(colnames(sample.cell.composition) == "inclusion.factor")]
              
              #Add percentage column:
              sample.cell.composition <- spread(sample.cell.composition, key = Basic_Celltype, value = cellcount)
              sample.cell.composition[is.na(sample.cell.composition)] <- 0
              sample.cell.composition$Total <- apply(sample.cell.composition[,-c(1:4)], MARGIN = 1, FUN = sum)
              sample.cell.composition <- gather(sample.cell.composition, key = "Basic_Celltype", value = "cellcount", -c(Lesion, Patient, Sequence_date, Lesion_visual_score, Total))
              sample.cell.composition$percentage <- round((sample.cell.composition$cellcount/sample.cell.composition$Total)*100, digits = 1)
              sample.cell.composition$Basic_Celltype <- factor(sample.cell.composition$Basic_Celltype, 
                                                               levels = levels(contact_derm@meta.data$Basic_Celltype))
              
              #Reformat for graph:
              sample.cell.composition$Lesion <- as.character(sample.cell.composition$Lesion)
              sample.cell.composition$Lesion[which(str_detect(sample.cell.composition$Patient, pattern = "^CB"))] <- "Healthy"
              sample.cell.composition$Lesion <- str_replace_all(sample.cell.composition$Lesion, pattern = "_", replacement = " ")
              sample.cell.composition$Lesion <- factor(sample.cell.composition$Lesion,
                                                       levels = c("Healthy", str_replace_all(levels(contact_derm@meta.data$Lesion), pattern = "_", replacement = " ")))
              
              sample.cell.composition$Lesion_visual_score <- factor(sample.cell.composition$Lesion_visual_score,
                                                                      levels = c(0,0.25,0.5,1,2))
              
              
              ggplot(sample.cell.composition, aes(x = Basic_Celltype, y = percentage, fill = Lesion)) +
                geom_boxplot() +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) + expand_limits(y=100) +
                labs(title = "Basic Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text = element_text(size = 15),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15),
                      axis.text.x = element_text(vjust = -0.5)) +
                labs(x = "Basic Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Basic Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
              
              
          
              ggplot(sample.cell.composition, aes(x = Basic_Celltype, y = percentage, color = Lesion)) +
                geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) + expand_limits(y=100) +
                labs(title = "Basic Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text = element_text(size = 15),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15),
                      axis.text.x = element_text(vjust = -0.5)) +
                labs(x = "Basic Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Basic Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
              
              #Build Custom Color Palette:
              color.palette <- c("#080808","#F50000")
              GetPalette = colorRampPalette(color.palette)
              
              ColorCount = as.numeric(length(unique(sample.cell.composition$Lesion_visual_score)))
              custom.colors <- GetPalette(ColorCount)
              
              ggplot(sample.cell.composition, aes(x = Lesion, y = percentage)) +
                geom_violin(width = 1.25) +
                geom_jitter(width = 0.1, aes(color = Lesion_visual_score)) +
                scale_color_manual(values = custom.colors) +
                theme_classic() +
                facet_wrap(~Basic_Celltype, scales = "free") +
                scale_y_continuous(expand = c(0,0)) +
                labs(title = "Basic Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text = element_text(size = 12),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(x = "Basic Celltype", y = "Percent of Total Cells", color = "Lesion visual score", title = "Basic Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                geom_text(aes(y=percentage * 1.1, label=""))

      ## Same analysis by Semi_Detailed_Celltype:
              
              sample.cell.composition <- DE.input@meta.data %>% group_by_at(c("Lesion", "Patient", "Semi_Detailed_Celltype", "Sequence_date", "Lesion_visual_score")) %>%
                                                                  summarise(n()) %>%
                                                                  rename_at("n()", ~ "cellcount") %>% as.data.frame()
              
              sample.cell.composition$Lesion_visual_score[which(is.na(sample.cell.composition$Lesion_visual_score))] <- 0
              
              
                                ##Only analyze samples with greater than 100 cells:
                                sample.cell.count <- contact_derm@meta.data %>% group_by_at(c("Lesion", "Patient", "Sequence_date")) %>%
                                  summarise(n()) %>%
                                  rename_at("n()", ~ "cellcount") %>% as.data.frame()
                                
                                sample.cell.count$inclusion.factor <- NA
                                sample.cell.count <- sample.cell.count %>% mutate(inclusion.factor = ifelse(cellcount >= 100, "yes", "no"))
                                sample.cell.composition <- full_join(sample.cell.composition, sample.cell.count[,-which(colnames(sample.cell.count) == "cellcount")])
                                sample.cell.composition <- sample.cell.composition %>% filter(inclusion.factor == "yes")
                                sample.cell.composition <- sample.cell.composition[,-which(colnames(sample.cell.composition) == "inclusion.factor")]
                                
              #Add percentage column:
              sample.cell.composition <- spread(sample.cell.composition, key = Semi_Detailed_Celltype, value = cellcount)
              sample.cell.composition[is.na(sample.cell.composition)] <- 0
              sample.cell.composition$Total <- apply(sample.cell.composition[,-c(1:4)], MARGIN = 1, FUN = sum)
              sample.cell.composition <- gather(sample.cell.composition, key = "Semi_Detailed_Celltype", value = "cellcount", -c(Lesion, Patient, Sequence_date,Lesion_visual_score, Total))
              sample.cell.composition$percentage <- round((sample.cell.composition$cellcount/sample.cell.composition$Total)*100, digits = 1)
              
              unique(sample.cell.composition$Semi_Detailed_Celltype)
              sample.cell.composition$Semi_Detailed_Celltype <- factor(sample.cell.composition$Semi_Detailed_Celltype, 
                                                               levels = levels(contact_derm@meta.data$Semi_Detailed_Celltype))
              #Reformat for graph:
              sample.cell.composition$Lesion <- as.character(sample.cell.composition$Lesion)
              sample.cell.composition$Lesion[which(str_detect(sample.cell.composition$Patient, pattern = "^CB"))] <- "Healthy"
              sample.cell.composition$Lesion <- str_replace_all(sample.cell.composition$Lesion, pattern = "_", replacement = " ")
              sample.cell.composition$Lesion <- factor(sample.cell.composition$Lesion,
                                                       levels = c("Healthy", str_replace_all(levels(contact_derm@meta.data$Lesion), pattern = "_", replacement = " ")))
              
              sample.cell.composition$Lesion_visual_score <- factor(sample.cell.composition$Lesion_visual_score,
                                                                    levels = c(0,0.25,0.5,1,2))
              
              
              
              ggplot(sample.cell.composition, aes(x = Semi_Detailed_Celltype, y = percentage, fill = Lesion)) +
                geom_boxplot() +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) + expand_limits(y=100) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text = element_text(size = 15),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15),
                      axis.text.x = element_text(vjust = -0.5)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
              
              ggplot(sample.cell.composition, aes(x = Semi_Detailed_Celltype, y = percentage, color = Lesion)) +
                geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width= 0.75)) +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) + expand_limits(y=100) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text = element_text(size = 15),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15),
                      axis.text.x = element_text(vjust = -0.5)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
              
              #Build Custom Color Palette:
              color.palette <- c("#080808","#F50000")
              GetPalette = colorRampPalette(color.palette)
              
              ColorCount = as.numeric(length(unique(sample.cell.composition$Lesion_visual_score)))
              custom.colors <- GetPalette(ColorCount)
             
              ggplot(sample.cell.composition, aes(x = Lesion, y = percentage)) +
                geom_violin(width = 1)+
                geom_jitter(width = 0.15, aes(color = Lesion_visual_score)) +
                scale_color_manual(values = custom.colors) +
                theme_classic() +
                facet_wrap(~Semi_Detailed_Celltype, scales = "free") +
                scale_y_continuous(expand = c(0,0)) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text.y = element_text(size = 10),
                      axis.text.x = element_text(size = 8),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion visual score", title = "Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                geom_text(aes(y=percentage * 1.1, label=""))

              
              total.cell.counts <- DE.input@meta.data %>% group_by(Semi_Detailed_Celltype, Lesion) %>% summarize(cellcount = n())
              
              Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
              
              ggplot(total.cell.counts, aes(x = Lesion, y = log2(cellcount), fill = Lesion, color = Lesion)) +
                geom_col()+
                facet_grid(Semi_Detailed_Celltype~., scales="free_y")+
                theme_classic()+
                scale_color_manual(values = Lesion.colors)+
                scale_fill_manual(values = Lesion.colors)
                
              semi_detailed_log2cellcount<- last_plot()
              ggsave(semi_detailed_log2cellcount, width = 25, height = 25,units = "cm",filename = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_heatmaps/semi_detailed_log2cellcount.pdf", limitsize = FALSE) 
              
              
    ## Same analysis by Detailed_Celltype:
              sample.cell.composition <- DE.input@meta.data %>% group_by_at(c("Lesion", "Patient", "Detailed_Celltype", "Sequence_date", "Lesion_severity_score")) %>%
                summarise(n()) %>%
                rename_at("n()", ~ "cellcount") %>% as.data.frame()
              
              sample.cell.composition$Lesion_severity_score[which(is.na(sample.cell.composition$Lesion_severity_score))] <- "-"
              
              
                    ##Only analyze samples with greater than 100 cells:
                    sample.cell.count <- DE.input@meta.data %>% group_by_at(c("Lesion", "Patient", "Sequence_date")) %>%
                      summarise(n()) %>%
                      rename_at("n()", ~ "cellcount") %>% as.data.frame()
                    
                    sample.cell.count$inclusion.factor <- NA
                    sample.cell.count <- sample.cell.count %>% mutate(inclusion.factor = ifelse(cellcount >= 100, "yes", "no"))
                    sample.cell.composition <- full_join(sample.cell.composition, sample.cell.count[,-which(colnames(sample.cell.count) == "cellcount")])
                    sample.cell.composition <- sample.cell.composition %>% filter(inclusion.factor == "yes")
                    sample.cell.composition <- sample.cell.composition[,-which(colnames(sample.cell.composition) == "inclusion.factor")]
                    
              #Add percentage column:
              sample.cell.composition <- spread(sample.cell.composition, key = Detailed_Celltype, value = cellcount)
              sample.cell.composition[is.na(sample.cell.composition)] <- 0
              sample.cell.composition$Total <- apply(sample.cell.composition[,-c(1:4)], MARGIN = 1, FUN = sum)
              sample.cell.composition <- gather(sample.cell.composition, key = "Detailed_Celltype", value = "cellcount", -c(Lesion, Patient, Sequence_date,Lesion_severity_score, Total))
              sample.cell.composition$percentage <- round((sample.cell.composition$cellcount/sample.cell.composition$Total)*100, digits = 1)
              
              sample.cell.composition$Detailed_Celltype <- factor(sample.cell.composition$Detailed_Celltype, 
                                                                       levels = levels(contact_derm@meta.data$Detailed_Celltype))
              
              #Reformat for graph:
              sample.cell.composition$Lesion <- as.character(sample.cell.composition$Lesion)
              sample.cell.composition$Lesion[which(str_detect(sample.cell.composition$Patient, pattern = "^CB"))] <- "Healthy"
              sample.cell.composition$Lesion <- str_replace_all(sample.cell.composition$Lesion, pattern = "_", replacement = " ")
              sample.cell.composition$Lesion <- factor(sample.cell.composition$Lesion,
                                                       levels = c("Healthy", str_replace_all(levels(contact_derm@meta.data$Lesion), pattern = "_", replacement = " ")))
              
              sample.cell.composition$Lesion_severity_score <- factor(sample.cell.composition$Lesion_severity_score,
                                                                    levels = levels(DE.input@meta.data$Lesion_severity_score))
              
              
              # Statistical T tests between Day2 Allergy vs Day 2 Irritant:
              library(rstatix)
              comparisons <- data.frame(reference = c("Irritant"),
                                        comparison = c("Day2 Allergy"))
              
              
              for (c in 1:nrow(comparisons)) {
                
                ref = comparisons$reference[c]
                comp = comparisons$comparison[c]
                
                #Paired T tests for significance:
                comparison.ttest <- sample.cell.composition %>% filter(Lesion %in% c(ref, comp)) %>% group_by(Detailed_Celltype) %>% rstatix::pairwise_t_test(
                                                                    percentage ~ Lesion, paired = FALSE, pool.sd = FALSE, var.equal = FALSE, 
                                                                    p.adjust.method = "BH") 
              
                
                #Fold Change:
                comparison.means <- sample.cell.composition %>% filter(Lesion %in% c(ref, comp)) %>% group_by(Detailed_Celltype, Lesion) %>% summarize(mean = mean(percentage))
                comparison.means <- comparison.means %>% pivot_wider(names_from = Lesion, values_from = mean)
                comparison.means <- comparison.means %>% mutate(foldchange = eval(as.symbol(comp))/eval(as.symbol(ref)))
                
                comparison.ttest <- cbind(comparison.ttest, comparison.means$foldchange[match(comparison.ttest$Detailed_Celltype, comparison.means$Detailed_Celltype)])
                colnames(comparison.ttest)[12] <- "Fold_change"
                
                comparison.ttest <- comparison.ttest %>% filter(!(is.nan(p.adj)))
                comparison.ttest$p.adj.signif <- factor(comparison.ttest$p.adj.signif, levels = c("ns", "*", "**", "***", "****"))
                
                comparison.ttest <- comparison.ttest %>% arrange(p.adj)
                
                #Label comparison conditions:
                comparison.ttest$Reference <- ref
                comparison.ttest$Comparison <- comp
                
                if (c == 1) {
                  cell.composition.stats <- comparison.ttest
                }
                
                if (c > 1) {
                  cell.composition.stats <- bind_rows(cell.composition.stats, comparison.ttest)
                }
                
              }

              to.plot <- cell.composition.stats$Detailed_Celltype[1:4]
              high.fc.cells <- cell.composition.stats %>% filter(Fold_change > 2) %>% dplyr::select(Detailed_Celltype) %>% unlist() %>% as.character()
              
              
              ## Linear Modeling Cell proportions between Day2 Allergy vs Day2 Irritant with Lesion Severity as a covariate:
              
              
              for (c in 1:length(unique(sample.cell.composition$Detailed_Celltype))) {
              
                  cell.data <- sample.cell.composition %>% filter(Detailed_Celltype == unique(sample.cell.composition$Detailed_Celltype)[c] & Lesion %in% c("Irritant", "Day2 Allergy"))

                  lm.model <- lm(percentage ~ Lesion + as.numeric(Lesion_severity_score), data = cell.data)
                  lm.coefficients <- as.data.frame(broom::tidy(lm.model))
                  lm.coefficients$Celltype <- unique(sample.cell.composition$Detailed_Celltype)[c]
              
                  lm.coefficients$term <- str_replace(lm.coefficients$term, pattern = "^Lesion", replacement = "Lesion_")
                  lm.coefficients$term[which(lm.coefficients$term == "as.numeric(Lesion_severity_score)")] <- "Lesion_severity_score"

                  colnames(lm.coefficients)[which(colnames(lm.coefficients) == "estimate")] <- "coefficient"
                  
                  if (c==1) {
                    cell.lm.stats <- lm.coefficients
                  }
                  if (c>1) {
                    cell.lm.stats <- bind_rows(cell.lm.stats, lm.coefficients)
                  }
                
              }
              
             allergy.lm.stats <- cell.lm.stats %>% filter(term == "Lesion_Day2 Allergy")
             allergy.lm.stats$lower_95_CI <- allergy.lm.stats$coefficient-1.96*allergy.lm.stats$std.error
             allergy.lm.stats$upper_95_CI <- allergy.lm.stats$coefficient+1.96*allergy.lm.stats$std.error
             
             #Compare to Lesion Severity Score Coefficient Estimates:
             lesion.score.lm.stats <- cell.lm.stats %>% filter(term == "Lesion_severity_score")
             lesion.score.lm.stats$lower_95_CI <- lesion.score.lm.stats$coefficient-1.96*lesion.score.lm.stats$std.error
             lesion.score.lm.stats$upper_95_CI <- lesion.score.lm.stats$coefficient+1.96*lesion.score.lm.stats$std.error
             
                 plot.min.lim <- min(allergy.lm.stats$lower_95_CI,lesion.score.lm.stats$lower_95_CI)
                 plot.max.lim <- max(allergy.lm.stats$upper_95_CI,lesion.score.lm.stats$upper_95_CI)
                 
                 
             
             
             # Plot Coefficient Estimates:
             
                   p1 <-ggplot(allergy.lm.stats, aes(x=Celltype, y=coefficient)) +
                     geom_errorbar(aes(ymin=lower_95_CI, ymax=upper_95_CI), 
                                   width = 0.2,size  = 0.5,
                                   position = "dodge") +
                     geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = "dashed") +
                     geom_point() + coord_flip() + theme_classic() +
                     scale_x_discrete(limits = rev(unique(allergy.lm.stats$Celltype))) +
                     scale_y_continuous(limits = c(-30,30)) +
                     labs(y = "Allergy vs Irritant Coefficient Estimate with 95% Confidence Interval", title = "Lesion type coefficient estimate in linear models of cell proportion for each celltype") +
                     theme(plot.title = element_text(hjust = 0.5))
             
                   
                   p2 <- ggplot(lesion.score.lm.stats, aes(x=Celltype, y=coefficient)) +
                     geom_errorbar(aes(ymin=lower_95_CI, ymax=upper_95_CI), 
                                   width = 0.2,size  = 0.5,
                                   position = "dodge") +
                     geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = "dashed") +
                     geom_point() + coord_flip() + theme_classic() +
                     scale_x_discrete(limits = rev(unique(allergy.lm.stats$Celltype))) +
                     scale_y_continuous(limits = c(-30,30)) +
                     labs(y = "Lesion Severity Score Coefficient Estimate with 95% Confidence Interval", title = "Lesion Severity Score coefficient estimate in linear models of cell proportion for each celltype") +
                     theme(plot.title = element_text(hjust = 0.5))
                   
                   p1+p2
                   linear.model.coefficient.plot <- last_plot()
                   #ggsave(linear.model.coefficient.plot, height = 25, width = 40, units = "cm", filename = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/cell.proprtion.lm.estimates.pdf")
                   
             # Statistical T tests between Day2 Allergy / Day 2 Irritant and Nonlesional to find common cell proportion changes:
                     library(rstatix)
                     comparisons <- data.frame(reference = c("Nonlesional", "Nonlesional"),
                                               comparison = c("Day2 Allergy", "Irritant"))
                     
                     
                     for (c in 1:nrow(comparisons)) {
                       
                       ref = comparisons$reference[c]
                       comp = comparisons$comparison[c]
                       
                       #Paired T tests for significance:
                       comparison.ttest <- sample.cell.composition %>% filter(Lesion %in% c(ref, comp)) %>% group_by(Detailed_Celltype) %>% rstatix::pairwise_t_test(
                         percentage ~ Lesion, paired = FALSE, pool.sd = FALSE, var.equal = FALSE, 
                         p.adjust.method = "BH") 
                       
                       
                       #Fold Change:
                       comparison.means <- sample.cell.composition %>% filter(Lesion %in% c(ref, comp)) %>% group_by(Detailed_Celltype, Lesion) %>% summarize(mean = mean(percentage))
                       comparison.means <- comparison.means %>% pivot_wider(names_from = Lesion, values_from = mean)
                       comparison.means <- comparison.means %>% mutate(foldchange = eval(as.symbol(comp))/eval(as.symbol(ref)))
                       
                       comparison.ttest <- cbind(comparison.ttest, comparison.means$foldchange[match(comparison.ttest$Detailed_Celltype, comparison.means$Detailed_Celltype)])
                       colnames(comparison.ttest)[12] <- "Fold_change"
                       
                       comparison.ttest <- comparison.ttest %>% filter(!(is.nan(p.adj)))
                       comparison.ttest$p.adj.signif <- factor(comparison.ttest$p.adj.signif, levels = c("ns", "*", "**", "***", "****"))
                       
                       comparison.ttest <- comparison.ttest %>% arrange(p.adj)
                       
                       #Label comparison conditions:
                       comparison.ttest$Reference <- ref
                       comparison.ttest$Comparison <- comp
                       
                       if (c == 1) {
                         cell.composition.stats <- comparison.ttest
                       }
                       
                       if (c > 1) {
                         cell.composition.stats <- bind_rows(cell.composition.stats, comparison.ttest)
                       }
                       
                     }
                     
                     high.fc.cells.irr <- cell.composition.stats %>% filter(group2 == "Irritant") %>% filter(Fold_change > 2 | Fold_change < 0.5) %>% dplyr::select(Detailed_Celltype) %>% unlist() %>% as.character()
                     high.fc.cells.all <- cell.composition.stats %>% filter(group2 == "Day2 Allergy") %>% filter(Fold_change > 2 | Fold_change < 0.5) %>% dplyr::select(Detailed_Celltype) %>% unlist() %>% as.character()
                     
                     high.fc.common.cells <- high.fc.cells.irr[which(high.fc.cells.irr %in% high.fc.cells.all)]
                     
                     
                     #filter high FC cells to make sure they are statistically significant from Nonlesional/Healthy in at least 1 condition
                     high.fc.common.cells <- cell.composition.stats %>% filter(Detailed_Celltype %in% high.fc.common.cells) %>% filter(p.adj <= 0.1) %>% dplyr::pull(Detailed_Celltype) %>% as.character()
                     #add KRT-wr because they are so close:
                     high.fc.common.cells <- c(high.fc.common.cells, "KRT-wr")
                     
                     
                     #filter out cells that were also significantly different between Irritant and Allergic:
                     high.fc.common.cells <- high.fc.common.cells[-which(high.fc.common.cells %in% high.fc.cells)]
                     common.cell.changes.df <- sample.cell.composition %>% filter(Detailed_Celltype %in% high.fc.common.cells)
                     
                     cell.composition.stats.to.plot <- cell.composition.stats %>% filter(Detailed_Celltype %in% high.fc.common.cells)
                     
                     common.cell.changes.df <- sample.cell.composition %>% filter(Detailed_Celltype %in% high.fc.common.cells)
                     
                     #Build Custom Color Palette:
                     color.palette <- c("#080808","#F50000")
                     GetPalette = colorRampPalette(color.palette)
                     
                     ColorCount = as.numeric(length(unique(common.cell.changes.df$Lesion_severity_score)))
                     custom.colors <- GetPalette(ColorCount)
                     
                     
                     ggplot(common.cell.changes.df, aes(x = Lesion, y = percentage)) +
                       geom_boxplot(color = "#C2C2C2", outlier.alpha = 0)+
                       geom_jitter(width = 0.2, height = 0, aes(color = Lesion_severity_score)) +
                       scale_color_manual(values = custom.colors) +
                       theme_classic() +
                       facet_wrap(~Detailed_Celltype, scales = "free") +
                       scale_y_continuous(expand = c(0.05,0.05)) +
                       labs(title = "Cell Compositions") +
                       theme(text = element_text(size = 18),
                             axis.text.y = element_text(size = 10),
                             axis.text.x = element_text(size = 8),
                             plot.title = element_text(hjust = 0.5),
                             legend.text=element_text(size=12),
                             legend.title = element_text(size = 15)) +
                       labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions") +
                       scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                       geom_text(aes(y=percentage * 1.1, label=""))
                     
                     common.CD.cell.proportion.plot <- last_plot()    

                     #ggsave(common.CD.cell.proportion.plot, height = 30, width = 30, units = "cm", filename = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/common.CD.cell.proportion.plot.pdf")
                     
              
              #Only plot celltypes with trend:
              #to.plot <- c("KRT-wr", "LC", "Myeloid1", "M2_mac", "cDC1", "CD4_conv1", "CD4_hsp", "CD4_cd161", "CD4_tcm", "Treg", "CD8_CD4_active", "CD8", "NK_ctl")
              #to.plot <- c("Mel1", "Mel2", "LC", "Myeloid_1", "Myeloid_ccl22", "M2_mac", "cDC1", "CD4_conv1", "CD4_hsp", "CD4_cd161", "CD4_tcm", "Treg", "CD8_CD4_active", "CD8", "NK_ctl")
              #to.plot <- c("KRT-b1", "KRT-b2", "KRT-sp", "KRT-g", "KRT-77", "KRT-mucl")
              #to.plot <- c("cDC1", "CD4-hspa6", "NK-ctl", "CD8", "CD8-CD4-active")
              
             
              sample.cell.composition.to.plot <- sample.cell.composition %>% filter(Detailed_Celltype %in% to.plot)
              
              ggplot(sample.cell.composition.to.plot, aes(x = Detailed_Celltype, y = percentage, fill = Lesion)) +
                geom_boxplot() +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) + expand_limits(y=50) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text = element_text(size = 10),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15),
                      axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions")
              
              ggplot(sample.cell.composition.to.plot, aes(x = Detailed_Celltype, y = percentage, color = Lesion)) +
                geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width= 0.75)) +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) + expand_limits(y=50) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text = element_text(size = 15),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15),
                      axis.text.x = element_text(vjust = -0.5)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5))
              
              ggplot(sample.cell.composition.to.plot, aes(x = Lesion, y = percentage, color = Lesion)) +
                geom_boxplot(color = "#C2C2C2") +
                geom_jitter(width = 0.15) +
                theme_classic() +
                facet_wrap(~Detailed_Celltype, scales = "free") +
                scale_y_continuous(expand = c(0,0)) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text.y = element_text(size = 10),
                      axis.text.x = element_text(size = 8),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                geom_text(aes(y=percentage * 1.1, label=""))
              
              cell.proportion.plot <- last_plot()
              
              #ggsave(cell.proportion.plot, height = 20, width = 50, units = "cm", filename = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/cell.proprtion.plot.pdf")
              
              #Build Custom Color Palette:
              color.palette <- c("#080808","#F50000")
              GetPalette = colorRampPalette(color.palette)
              
              ColorCount = as.numeric(length(unique(sample.cell.composition.to.plot$Lesion_severity_score)))
              custom.colors <- GetPalette(ColorCount)
              
              ggplot(sample.cell.composition.to.plot, aes(x = Lesion, y = percentage)) +
                geom_boxplot(color = "#C2C2C2", outlier.alpha = 0)+
                geom_jitter(width = 0.2, height = 0, aes(color = Lesion_severity_score)) +
                scale_color_manual(values = custom.colors) +
                theme_classic() +
                facet_wrap(~Detailed_Celltype, scales = "free") +
                scale_y_continuous(expand = c(0.05,0.05)) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text.y = element_text(size = 10),
                      axis.text.x = element_text(size = 8),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                geom_text(aes(y=percentage * 1.1, label=""))
              
              cell.proportion.plot <- last_plot()
              #ggsave(cell.proportion.plot, height = 30, width = 20, units = "cm", filename = "/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/cell.proprtion.plot2.pdf")

              ggplot(sample.cell.composition, aes(x = Lesion, y = percentage)) +
                geom_jitter(width = 0.15, aes(color = Lesion_severity_score)) +
                scale_color_manual(values = custom.colors) +
                theme_classic() +
                facet_grid(Detailed_Celltype~., scales = "free") +
                scale_y_continuous(expand = c(0,0)) +
                labs(title = "Cell Compositions") +
                theme(text = element_text(size = 18),
                      axis.text.y = element_text(size = 10),
                      axis.text.x = element_text(size = 8),
                      plot.title = element_text(hjust = 0.5),
                      legend.text=element_text(size=12),
                      legend.title = element_text(size = 15)) +
                labs(x = "Celltype", y = "Percent of Total Cells", color = "Lesion", title = "Cell Compositions") +
                scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
                geom_text(aes(y=percentage * 1.1, label=""))+
                geom_hline(yintercept = 0)
              
              
## Composition analysis with propeller:  (Log transformation of proportions didn't work well for my data. . . I did not use this)
library(speckle)
library(limma)
              
              #Define cell compositions and filter to samples with 100 or greater cells:
              sample.cell.composition <- DE.input@meta.data %>% group_by_at(c("Lesion", "Patient", "Detailed_Celltype", "Sequence_date", "Lesion_visual_score")) %>%
                summarise(n()) %>%
                rename_at("n()", ~ "cellcount") %>% as.data.frame()
              
              sample.cell.composition$Lesion_visual_score[which(is.na(sample.cell.composition$Lesion_visual_score))] <- 0
              
              
              ##Only analyze samples with greater than 100 cells:
              sample.cell.count <- DE.input@meta.data %>% group_by_at(c("Lesion", "Patient", "Sequence_date")) %>%
                summarise(n()) %>%
                rename_at("n()", ~ "cellcount") %>% as.data.frame()
              
              sample.cell.count$inclusion.factor <- NA
              sample.cell.count <- sample.cell.count %>% mutate(inclusion.factor = ifelse(cellcount >= 100, "yes", "no"))
              sample.cell.composition <- full_join(sample.cell.composition, sample.cell.count[,-which(colnames(sample.cell.count) == "cellcount")])
              sample.cell.composition <- sample.cell.composition %>% filter(inclusion.factor == "yes")
              sample.cell.composition <- sample.cell.composition[,-which(colnames(sample.cell.composition) == "inclusion.factor")]
              
              sample.cell.composition$samples <- paste(sample.cell.composition$Patient, sample.cell.composition$Sequence_date,  sep = "__")

              
              #Setup for Propeller:
              cell_info_for_propeller <- DE.input@meta.data[which(colnames(DE.input@meta.data) %in% c("Lesion", "Patient", "Detailed_Celltype", "Sequence_date"))]
              cell_info_for_propeller$sample <- paste(cell_info_for_propeller$Patient, cell_info_for_propeller$Sequence_date, sep = "__")
              cell_info_for_propeller$sample_lesion <- paste(cell_info_for_propeller$Patient, cell_info_for_propeller$Sequence_date,cell_info_for_propeller$Lesion, sep = "__")
              
              cell_info_for_propeller <- cell_info_for_propeller %>% filter(sample %in% sample.cell.composition$samples)
              
              
              # Analyze by Propeller:
              prop.list <- getTransformedProps(clusters = cell_info_for_propeller$Detailed_Celltype, 
                                           sample = cell_info_for_propeller$sample_lesion, 
                                           transform="logit")
              

              grp <- data.frame(cbind(as.character(cell_info_for_propeller$Lesion), cell_info_for_propeller$sample_lesion)) %>% distinct()
              
              rownames(grp) <- grp$X2
              grp <- grp[colnames(prop.list$Proportions),"X1"]
              design <- model.matrix(~0+grp)
              colnames(design) <- str_replace(colnames(design), pattern = "^grp", replacement = "")
              design <- design[,c("Nonlesional", "Irritant","Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy")]
              
              contrasts <- c(0,1, 0, -1,0)
              ttest.stats <- propeller.ttest(prop.list = prop.list, design = design, contrasts = contrasts, trend = TRUE, sort = TRUE, robust=TRUE)
              
              #Plot Log-transformed proportions:
             transformed.props <- as.data.frame(t(prop.list$TransformedProps))
             transformed.props$Lesion <- sapply(transformed.props$sample, FUN = function(x) unlist(str_split(x, pattern = "__"))[3])
             transformed.props$Lesion <- factor(transformed.props$Lesion,
                                                levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy"))
             
             ggplot(transformed.props, aes(x = Lesion, y = Freq)) +
               geom_jitter(width = 0.15, height = 0, size = 0.5) +
               scale_color_manual(values = custom.colors) +
               theme_classic() +
               facet_wrap(~clusters, scales = "free") +
               labs(title = "Logit Transformed Cell Frequencies") +
               theme(text = element_text(size = 18),
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 8),
                     plot.title = element_text(hjust = 0.5),
                     legend.text=element_text(size=12),
                     legend.title = element_text(size = 15)) +
               labs(x = "Celltype", y = "Freq") +
               scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
               geom_text(aes(y=Freq * 1.1, label=""))
             
              plotCellTypeMeanVar(props$Counts)
              plotCellTypePropsMeanVar(props$Counts)
              

              propeller.anova(prop.list=props, design=design, 
                              robust=TRUE, trend=FALSE, sort=TRUE)
              
                                                #Bulk Aggregate DE Analysis:
                                                    #Only Run DE on samples with visual reaction scores of 1 or greater:
    
                                                    contact_derm@meta.data <- contact_derm@meta.data %>% mutate(DE_inclusion_factor = case_when(Lesion == "Nonlesional" ~ "yes",
                                                                                                                                              Lesion == "Acetone_Vehicle" & !Patient %in% c("HC8", "HC12") ~ "yes",
                                                                                                                                              Lesion == "Irritant" & Lesion_visual_score >= 1 ~ "yes",
                                                                                                                                              Lesion == "Day2_Allergy" & Lesion_visual_score >= 1 ~ "yes",
                                                                                                                                              Lesion == "Day4_Allergy" & Lesion_visual_score >= 1 ~ "yes",
                                                                                                                                              TRUE ~ "no"))
                                                          
                                                    
                                                    DE.input <- contact_derm
                                                    Idents(DE.input) <- "DE_inclusion_factor"
                                                    DE.input <- subset(DE.input, idents = "yes")
                                                    Idents(DE.input) <- "Semi_Detailed_Celltype"
                                                    
                                                    #Metadata formating:
                                                        DE.input@meta.data$Blister_date <- as.character(DE.input@meta.data$Blister_date)
                                                        DE.input@meta.data$Blister_date[which(is.na(DE.input@meta.data$Blister_date))] <- "NA"
                                                        
                                                        DE.input@meta.data$Sequence_instrument <- as.character(DE.input@meta.data$Sequence_instrument)

                                                        DE.input@meta.data <- DE.input@meta.data %>% mutate(BatchDE = case_when(Sequence_instrument %in% c("NS500602", "NB501205", "A00197") ~ "Batch1",
                                                                                                                                Sequence_instrument %in% c("VH00230", "A00439") ~ "Batch2"))
                                                        DE.input@meta.data$BatchDE <- factor(DE.input@meta.data$BatchDE,
                                                                                             levels = c("Batch1", "Batch2"))
                                                        
                                                        DE.input@meta.data$Lesion <- as.character(DE.input@meta.data$Lesion)
                                                        DE.input@meta.data <- DE.input@meta.data %>% mutate(Lesion_combos = case_when(Lesion %in% c("Nonlesional", "Irritant", "Acetone_Vehicle") ~ "Control",
                                                                                                                                      Lesion == "Day2_Allergy" ~ "Allergy",
                                                                                                                                      Lesion == "Day4_Allergy" ~ "Allergy"))
                                                        DE.input@meta.data$Lesion_combos <- factor(DE.input@meta.data$Lesion_combos,
                                                                                                   levels = c("Control", "Allergy"))
                                                  
                                                    #Aggregate as pseudobulks:
                                                        DE.input <- calc_agg_bulk_srt_new(DE.input, aggregate_by = c("Lesion", "Patient", "Blister_date", "Sequence_instrument"))
                                                        bulk_data <- DE.input@assays$RNA@misc$aggregate.bulk
                                                        
                                                          #code for old method:
                                                                  #bulk_data <- as.data.frame(t(bulk_data$RNA))
                                                                  
                                                                  #bulk_data <- rownames_to_column(bulk_data, var = "rowname")
                                                        
                                                                  #bulk_data$Lesion <- sapply(bulk_data$rowname, FUN = function(x) unlist(str_split(x, pattern = "__"))[1])
                                                                  #bulk_data$Patient <-sapply(bulk_data$rowname, FUN = function(x) unlist(str_split(x, pattern = "__"))[2])
                                                                  #bulk_data$Blister_date <- sapply(bulk_data$rowname, FUN = function(x) unlist(str_split(x, pattern = "__"))[3])
                                                                  
                                                                  #bulk_data <- bulk_data[,-which(colnames(bulk_data) == "rowname")]
                                                        
                                                          sample.lesion.info <- DE.input@meta.data %>% group_by_at(c("Lesion", "Patient", "Blister_date", "Lesion_visual_score", "Lesion_combos", "Sequence_instrument")) %>%
                                                                                                          summarise(n()) %>%
                                                                                                          rename_at("n()", ~ "cellcount") %>% as.data.frame()
                                                          sample.lesion.info$Blister_date <- as.character(sample.lesion.info$Blister_date)
                                                          
                                                        bulk_data <- full_join(bulk_data, sample.lesion.info)
                                                        bulk_data$Lesion_visual_score[which(is.na(bulk_data$Lesion_visual_score))] <- 0

                                                      #Batch Correct:
                                                        library(Harman)
                                                        metadata <- bulk_data[,-which(colnames(bulk_data) %in% rownames(DE.input))]
                                                        bulk_counts <- bulk_data[,which(colnames(bulk_data) %in% rownames(DE.input))]
                                                        
                                                        
                                                        harman.object <- harman(datamatrix = t(bulk_counts), 
                                                                                     expt = as.vector(metadata$Lesion), 
                                                                                     batch = as.vector(metadata$Sequence_instrument), 
                                                                                     limit = 0.95)
                                                        plot(harman.object)
                                                        
                                                        harman.corrected.df <- bulk_data
                                                        harman.corrected.df[,which(colnames(harman.corrected.df) %in% rownames(DE.input))] <- t(reconstructData(harman.object))
                                                        
                                                        
                                                        
                                                      #Bulk DE:
                                                        library(edgeR)
                                                        
                                                        bulk_counts <- t(bulk_data[,-which(colnames(bulk_data) %in% colnames(sample.lesion.info))])
                                                        name.vector <- paste(bulk_data$Lesion, bulk_data$Patient, bulk_data$Blister_date, sep = "_")
                                                        name.vector <- str_replace(name.vector, pattern = "_Vehicle", replacement = "")
                                                        colnames(bulk_counts) <- name.vector
                                                        
                                                        bulk.meta.data <- bulk_data[which(colnames(bulk_data) %in% colnames(sample.lesion.info))]
                                                        
                                                        d <- DGEList(counts = bulk_counts,
                                                                           group = factor(bulk.meta.data$Lesion,
                                                                                          levels = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy")))
                                                        
                                                        dim(d.full)
                                                        dim(d)
                                                        #Filter:
                                                              d.full <- d # keep the old one in case we mess up
                                                              keep <- rowSums(cpm(d)>5) >= 2
                                                              d <- d[keep,]
                                                              
                                                        #Normalize & Scale:
                                                        d <- calcNormFactors(d)
                                                      
                                                        #Estimate Dispersion and Generalized Linear Model:  
                                                        d1 <- estimateCommonDisp(d, verbose=T)
                                                        
                                                        design.mat <- model.matrix(~ 0 + d$samples$group)
                                                        colnames(design.mat) <- levels(d$samples$group)
                                                        d2 <- estimateGLMCommonDisp(d,design.mat)
                                                        d2 <- estimateGLMTrendedDisp(d2,design.mat)
                                                        # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
                                                        # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
                                                        d2 <- estimateGLMTagwiseDisp(d2,design.mat)
                                                        
                                                        fit <- glmFit(d2, design.mat)
                                                        colnames(design.mat)
                                                        
                                                        #Run DE:
                                                        comparisons <- as.character(sort(unique(DE.input@meta.data$Lesion)))
                                                        comparisons <- expand.grid(comparisons, comparisons, stringsAsFactors = FALSE)
                                                        colnames(comparisons) <- c("reference", "comparison")
                                                        comparisons <- comparisons[c(5,10,15,20,9,14,6,11),]
                                                        
                                                        bulk.de <- data.frame(matrix(ncol = 8, nrow = 0))
                                                        colnames(bulk.de) <- c("Gene","logFC", "logCPM", "LR", "PValue", "FDR", "Reference", "Comparison")
                                                        bulk.de[,2:6] <-lapply(bulk.de[,1:5], MARGIN = 2, FUN = as.numeric)
                                                        bulk.de[,c(1,7,8)] <-lapply(bulk.de[,c(1,7,8)], MARGIN = 2, FUN = as.character)
                                                        
                                                        for (r in 1:nrow(comparisons)){
                                                          ref = comparisons$reference[r]
                                                          comp = comparisons$comparison[r]
                                                         
                                                          contrast = c(0,0,0,0,0)
                                                          names(contrast) = c("Nonlesional", "Irritant", "Acetone_Vehicle", "Day2_Allergy", "Day4_Allergy")
                                                        
                                                          contrast <- contrast - as.numeric(names(contrast) %in% c(ref))
                                                          contrast <- contrast + as.numeric(names(contrast) %in% c(comp))
                                                          
                                                          lrt <- glmLRT(fit, contrast=contrast)
                                                          de <- as.data.frame(edgeR::topTags(lrt, p.value = 1, n = Inf, sort.by = "logFC"))
                                                          
                                                          de <- rownames_to_column(de, var = "Gene")
                                                          de$Reference <- ref
                                                          de$Comparison <- comp
                                                          
                                                          bulk.de <- bind_rows(bulk.de, de)
                                                          
                                                        }
                                                        
                                                        bulk.de <- bulk.de %>% mutate(ref_comparison = paste(Reference, Comparison, sep = "__"))
                                                        
                                                        genes.to.exclude <- rownames(contact_derm)[which(str_detect(rownames(contact_derm), pattern = "^MT-|^MT.|^RPL|^RPS|^AC[:digit:]"))]
                                                        
                                                        significant.bulk.de.genes <- bulk.de %>% filter(!Gene %in% genes.to.exclude) %>% filter(logFC > 1 & FDR < 0.5) %>%
                                                                                                filter(!ref_comparison %in% c("Nonlesional__Acetone_Vehicle", "Nonlesional__Irritant")) %>%
                                                                                                group_by(ref_comparison) %>%
                                                                                                top_n(n = 5000, wt = (logFC)) %>%
                                                                                                dplyr::pull(Gene) %>% unique()
                                                  
                                                      #PCA & UMAP:  
                                                      bulk_pca_object <- prcomp(harman.corrected.df[,which(colnames(harman.corrected.df) %in% significant.bulk.de.genes)])
                                                      
                                                      bulk_pca <- cbind(harman.corrected.df[,which(colnames(harman.corrected.df) %in% colnames(sample.lesion.info))], bulk_pca_object$x)

                                                      set.seed(100)
                                                      bulk.umap <- umap(bulk_pca[,which(colnames(bulk_pca) %in% c(paste("PC", seq(1:20), sep = "")))])
                                                      colnames(bulk.umap$layout) <- c("UMAP_x", "UMAP_y")
                                                      
                                                      bulk_pca_umap <- cbind(bulk_pca, bulk.umap$layout, bulk_data[,which(colnames(bulk_data) %in% significant.bulk.de.genes)])
                                                      
                                                      bulk_pca_umap$Lesion_visual_score <- factor(bulk_pca_umap$Lesion_visual_score,
                                                                                                  levels = c(0,0.25,0.5,1,2))
                                                      
                                                      bulk_pca_umap$Lesion <- factor(bulk_pca_umap$Lesion,
                                                                                levels = levels(contact_derm@meta.data$Lesion))
                                                      
                                                      bulk.pca.long <- bulk_pca %>% pivot_longer(cols = which(str_detect(colnames(bulk_pca), pattern = "^PC")), names_to = "PC", values_to = "pc_val")
                                                      bulk.pca.long$PC <- factor(bulk.pca.long$PC,
                                                                                 levels = paste("PC", seq(1:length(unique(bulk.pca.long$PC))), sep = ""))
                                                      bulk.pca.long$Lesion <- factor(bulk.pca.long$Lesion,
                                                                                     levels = levels(contact_derm@meta.data$Lesion))
                                                      
                                                      bulk.pca.long %>% filter(PC %in% paste("PC", seq(1:30),sep = "")) %>%
                                                        ggplot(., aes(x = PC, y = pc_val, color = Lesion)) +
                                                        geom_point(position = position_dodge(width = 0.85))+
                                                        theme_classic()
                                                      
                                                      bulk.pca.long %>% filter(PC %in% paste("PC", seq(1:30),sep = "")) %>%
                                                        ggplot(., aes(x = PC, y = pc_val, color = Lesion)) +
                                                        geom_boxplot()+
                                                        theme_classic()
                                                    
                                                      ggplot(bulk_pca_umap, aes(x = PC1, y = PC2, color = Lesion)) +
                                                        geom_point(size = 2) +
                                                        theme_classic()
                                                      
                                                      ggplot(bulk_pca_umap, aes(x = UMAP_x, y = UMAP_y, color = Lesion)) +
                                                        geom_point(size = 2) +
                                                        theme_classic()
                                                      
                                                      #Build Custom Color Palette:
                                                      color.palette <- c("#080808","#F50000")
                                                      GetPalette = colorRampPalette(color.palette)
                                                      
                                                      ColorCount = as.numeric(length(unique(bulk_pca_umap$Lesion_visual_score)))
                                                      custom.colors <- GetPalette(ColorCount)
                                                      
                                                      ggplot(bulk_pca_umap, aes(x = UMAP_x, y = UMAP_y, color = Lesion_visual_score)) +
                                                                    geom_point(size = 2) +
                                                                    theme_classic()+
                                                                    scale_color_manual(values = custom.colors) +
                                                                    facet_wrap(~Lesion)
                                                      
                                                      ggplot(bulk_pca_umap, aes(x = UMAP_x, y = UMAP_y, color = Lesion_visual_score, shape = Lesion)) +
                                                        geom_point(size = 4) +
                                                        theme_classic()+
                                                        scale_color_manual(values = custom.colors) +
                                                        scale_shape_manual(values = c(16, 8, 0, 7, 12))
                                                      
                                                      
          
                
#Bulk Cytokine Barplots:   
      contact_derm <- calc_agg_bulk_srt(contact_derm, aggregate_by = c("Lesion", "Semi_Detailed_Celltype"), group_by = "Lesion")
      
      plot_violin_srt(contact_derm, gene = "IFNG", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
      plot_violin_srt(contact_derm, gene = "IL4", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
      plot_violin_srt(contact_derm, gene = "IL5", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
      plot_violin_srt(contact_derm, gene = "IL13", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
      
      
      plot_barplot_srt(contact_derm, gene = "IFNG", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      IFG.weighted.barplot <- last_plot()
      plot_barplot_srt(contact_derm, gene = "IL4", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      IL4.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "IL13", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      IL13.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "IL5", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      IL5.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "CXCL11", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      CXCL11.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "CXCL10", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      CXCL10.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "CXCL9", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      CXCL9.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "TNFRSF11B", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      TNFRSF11B.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "CCL2", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = TRUE)
      CCL2.weighted.barplot <- last_plot()
      
      plot_barplot_srt(contact_derm, gene = "TNF", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion", cell_proportion_weighted = FALSE)
      
      
      #ggsave(IFG.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/IFNG.weighted.barplot.pdf")
      #ggsave(IL4.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/IL4.weighted.barplot.pdf")
      #ggsave(IL13.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/IL13.weighted.barplot.pdf")
      #ggsave(IL5.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/IL5.weighted.barplot.pdf")
      
      #ggsave(CXCL9.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/CXCL9.weighted.barplot.pdf")
      #ggsave(CXCL10.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/CXCL10.weighted.barplot.pdf")
      #ggsave(CXCL11.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/CXCL11.weighted.barplot.pdf")
      #ggsave(CCL2.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/CCL2.weighted.barplot.pdf")
      #ggsave(TNFRSF11B.weighted.barplot, width = 40, height = 20, units = "cm", file = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Plots/TNFRSF11B.weighted.barplot.pdf")
      
      
      
# Supplemental Cytokine Violin Bar Plots:
      
      Chemokines.in.data <- rownames(DE.input@assays$SCT)[which(str_detect(rownames(DE.input@assays$SCT), pattern = "^CCL|^CXCL"))]
      Chemokines.in.data <- sort(Chemokines.in.data)
      Other.cytokines.in.data <- rownames(DE.input@assays$SCT)[which(str_detect(rownames(DE.input@assays$SCT), pattern = "^IFN|^IL"))]
      Other.cytokines.in.data <- Other.cytokines.in.data[-which(str_detect(Other.cytokines.in.data, pattern = "R"))]
      Other.cytokines.in.data <- sort(Other.cytokines.in.data)
      
      
      Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
      
      for (c in 1:length(Chemokines.in.data)) {
        g <- Chemokines.in.data[c]
        
        plot <- plot_violin_bar_srt(DE.input, gene = g, color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
        
        assign(paste(g, "_plot", sep = ""), plot)
        
        if (c == 1) {
          plot.list <- list(get(paste(g, "_plot", sep = "")))
        }
        if (c > 1) {
          plot.list[[c]] <- get(paste(g, "_plot", sep = ""))
        }
        
        
      }
      
      Chemokines.finished.plot <- ggpubr::ggarrange(plotlist = plot.list, common.legend = TRUE)
      
      
      for (c in 1:length(Other.cytokines.in.data)) {
        g <- Other.cytokines.in.data[c]
        
        plot <- plot_violin_bar_srt(DE.input, gene = g, color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
        
        assign(paste(g, "_plot", sep = ""), plot)
        
        if (c == 1) {
          plot.list2 <- list(get(paste(g, "_plot", sep = "")))
        }
        if (c > 1) {
          plot.list2[[c]] <- get(paste(g, "_plot", sep = ""))
        }
        
        
      }
      
      Other.Cytokines.finished.plot <- ggpubr::ggarrange(plotlist = plot.list2, common.legend = TRUE)
      
      Other.genes.for.violin.barplot <- c("KIT", "L1CAM")
      
      for (c in 1:length(Other.genes.for.violin.barplot)) {
        g <- Other.genes.for.violin.barplot[c]
        
        plot <- plot_violin_bar_srt(DE.input, gene = g, color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
        
        assign(paste(g, "_plot", sep = ""), plot)
        
        if (c == 1) {
          plot.list3 <- list(get(paste(g, "_plot", sep = "")))
        }
        if (c > 1) {
          plot.list3[[c]] <- get(paste(g, "_plot", sep = ""))
        }
        
        
      }
      
      Other.Genes.finished.plot <- ggpubr::ggarrange(plotlist = plot.list3, common.legend = TRUE)
      
      #ggsave(Chemokines.finished.plot, width = 100, height = 100,units = "cm",filename = "Plots/barplots/All.Chemokines.violin.barplot.png", limitsize = FALSE)  
      #ggsave(Other.Cytokines.finished.plot, width = 100, height = 100,units = "cm",filename = "Plots/barplots/All.other.cytokines.violin.barplot.png", limitsize = FALSE)  
      #ggsave(Other.Genes.finished.plot, width = 30, height = 15,units = "cm",filename = "Plots/barplots/Other.genes.of.interest.violin.barplot.png", limitsize = FALSE)  
      
      plot_violin_bar_srt(DE.input, gene = "CXCL14", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
      plot_violin_srt(DE.input, gene = "IL16", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", colors = Lesion.colors)
      
      DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Lesion"))
      plot_barplot_srt(DE.input, gene = "CXCL14", color_by = "Lesion", number_labels = F, colors = Lesion.colors, cell_proportion_weighted = FALSE)
      CXCL14.bulk.barplot <- last_plot()
      
      DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Semi_Detailed_Celltype", "Lesion"))
      
      
      plot_violin_bar_srt(DE.input, gene = "CD207", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
      CD207.violin.barplot <- last_plot()
      
      
      plot_violin_bar_srt(DE.input, gene = "CXCL14", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
      CXCL14.violin.barplot <- last_plot()
      
      
      plot_violin_bar_srt(DE.input, gene = "CXCL8", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
      
      CXCL8.violin.barplot <- last_plot()
      
      plot_violin_bar_srt(DE.input, gene = "CCL17", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
      CCL17.violin.barplot <- last_plot()
      
      plot_violin_bar_srt(DE.input, gene = "HLA.G", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", plot_patient_pos_frac = FALSE, colors = Lesion.colors)
      KRT17.violin.barplot <- last_plot()
      
      
      #ggsave(CXCL14.bulk.barplot, width = 10, height = 22,units = "cm",filename = "Plots/barplots/CXCL14.bulk.barplot.pdf", limitsize = FALSE)  
      #ggsave(CD207.violin.barplot, width = 50, height = 35,units = "cm",filename = "Plots/barplots/CD207.violin.barplot.png", limitsize = FALSE)  
      #ggsave(CXCL14.violin.barplot, width = 50, height = 45,units = "cm",filename = "Plots/barplots/CXCL14.violin.barplot.png", limitsize = FALSE)  
      

      
# Single-cell Cytokine Heatmaps:        
      cytokine_genes <- unlist(fread("gene_lists/cytokine_genes.tsv"))
              
      set.seed(100)
      cytokines <- plot_heatmap_srt(contact_derm,
                                    genes = cytokine_genes,
                                    facet_by = c("Lesion", "Semi_Detailed_Celltype"),
                                    type = "single_cell",
                                    log_scale_base = 2,
                                    ceiling = 5,
                                    col_names = FALSE,
                                    gene_names = TRUE,
                                    cluster_by = "row",
                                    text_sizes = c(15, 10, 10, 10, 5, 5,5),
                                    color_pal = colorRampPalette(c("#E5E5E5","#E5E5E5", "#E5E5E5","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256))
                                    
# Bulk Cytokine Heatmaps:
              
    
    #By Lesion and Celltype:
        cytokine_genes <- unlist(fread("gene_lists/cytokine_genes.tsv"))

        contact_derm <- calc_agg_bulk_srt(contact_derm, aggregate_by = c("Lesion", "Semi_Detailed_Celltype"))
        
        plot_heatmap_srt(contact_derm,
                         genes = cytokine_genes,
                         type = "bulk",
                         facet_by = c("Lesion"),
                         scale_group = "Semi_Detailed_Celltype",
                         text_angle = 90,
                         text_sizes = c(20, 10, 7, 10, 5, 5, 5),
                         color_pal = colorRampPalette(c("#97A0CD","#DBDDE8","#E5E5E5","#E5E5E5", "#FFC900","#EC1F1F"))(256),
                         log_scale_base = 2,
                         ceiling = 2,
                         floor = -2)
        
        #With Gene Order:
        plot_heatmap_srt(contact_derm, 
                         genes = cytokines_reorder, 
                         facet_by = "Semi_Detailed_Celltype",
                         type = "bulk", 
                         cluster_by = FALSE, 
                         pdf_format = "tile", 
                         scale_by = "row", 
                         gene_names = T, 
                         col_names = TRUE, 
                         title = "Cytokines",
                         text_angle = 90,
                         text_sizes = c(20, 10, 7, 10, 5, 5, 5),
                         color_pal = colorRampPalette(c("#97A0CD", "#FFFFFF","#FFFFFF", "#FFC900","#EC1F1F"))(256),
                         log_scale_base = 2,
                         ceiling = 3)
        
        plot_violin_srt(contact_derm, gene = "CXCL13", facet_by = "Semi_Detailed_Celltype", color_by = "Lesion")
        
        
        
# Cluster GSEA Analysis:
library(clusterProfiler)

      #Treg gene set:
            
            Treg.enriched.genes.ref1 <- c("SH3BGRL2", "PRNP", "ITGAE1", "KLRG1", "CTLA4", "RGS1", "SLC22A2", "TNFRSF4", "TNFRSF9", "IL2RB", "IL2RA", "NRP1", "GPR83", "NT5E",
                                     "RGS16", "SOCS2", "TIAM1", "MAP3K8", "PLAGL1", "IKZF4", "IRF4", "IKZF2", "MDFIC", "FOXP3", "DUSP4")
            
            Treg.enriched.genes.ref2 <- c("FOXP3", "DUSP4", "SDC4", "NINJ2", "PTTG1", "TIAF1", "TRIB1", "S100A10", "GBP2", "GATA3", "IL2RA", "BHLHB2", "CEB1", "CTLA4", 
                                          "TFRC", "HLA-DMA", "AKAP2", "TNFRSF1B","CCR5", "GPR2", "IL2RB", "SHMt2", "HLA-DRB1", "HLA-DRB3", "TP53INP1", "GBP5", "EPSTI1", 
                                          "LGALS3", "SLAMF1", "TRAF1", "LGALS1", "S100A4", "G1P2")
            
            Treg.reduced.genes.ref1 <- c("IGFBP4", "ITGB3", "VIPR1", "KLRD1", "IL1RL2", "APP", "RGMB", "TGFBR3", "PTGER2", "SEMA4F", "SGK1", "MYO10", "APPL2", "ARHGAP29",
                                        "HS3ST3B1", "RAPGEF4", "ST3GAL6", "PDE3B", "ENC1", "POLE2")
            Treg.reduced.genes.ref2 <- c("SATB1", "PIM1", "ACTN1", "STAT4", "ID2", "NELL2", "SLC40A1", "IL1RL2","DGKA", "ITGB2", "STAT6", "GZMA", "MYC", "FHIT", "TCF7", 
                                         "IL7R", "CCF7", "PITPNC1", "RBMS1", "XBP1", "GZMK", "TNFSF5", "TRGV9", "CD81", "CNOT2", "CCL5", "NOSIP", "IFITM1", "PECAM1", "TNFRSF10B")
              
            Treg.enriched.genes <- unique(c(Treg.enriched.genes.ref1, Treg.enriched.genes.ref2))
            Treg.reduced.genes <- unique(c(Treg.reduced.genes.ref1, Treg.reduced.genes.ref2))

                #DE:
                  #DE_Treg_vs_otherCD4s <- FindMarkers(contact_derm, ident.1 = "Treg", ident.2 = c("CD4_conv1", "CD4_hsp", "CD4_cd161","CD4_tcm", "CD4_trm"), min.pct = 0, logfc.threshold = 0)
                  #DE_Treg_vs_otherCD4s <- rownames_to_column(DE_Treg_vs_otherCD4s, var = "Gene")
                  #DE_Treg_vs_otherCD4s <- DE_Treg_vs_otherCD4s %>% arrange(desc(avg_log2FC))
                  #save(DE_Treg_vs_otherCD4s, file = "Data/DE_Treg_vs_otherCD4s.Rdata")
                  
                #GSEA:
                  load("Data/DE_Treg_vs_otherCD4s.Rdata")
                  temp.genelist <- DE_Treg_vs_otherCD4s$avg_log2FC
                  names(temp.genelist) <- DE_Treg_vs_otherCD4s$Gene
                  
                  universe_id <- bitr(names(temp.genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                  universe_id$ENTREZID <- as.numeric(universe_id$ENTREZID)
                  universe_id$avg_log2FC <- temp.genelist[match(universe_id$SYMBOL, names(temp.genelist))]

                  gsea_set <- as.vector(universe_id$avg_log2FC)
                  names(gsea_set) <- universe_id$ENTREZID
                  
                  gsea_set <- sort(gsea_set, decreasing = TRUE)
                  
                  #Term2gene dataframe for custom Gene lists:
                          Treg.enriched.entrezid <- bitr(Treg.enriched.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                          Treg.enriched.entrezid$ENTREZID <- as.numeric(Treg.enriched.entrezid$ENTREZID)
                          Treg.enriched.entrezid$term <- "Treg_enriched"
                          
                          Treg.reduced.entrezid <- bitr(Treg.reduced.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                          Treg.reduced.entrezid$ENTREZID <- as.numeric(Treg.reduced.entrezid$ENTREZID)
                          Treg.reduced.entrezid$term <- "Treg_reduced"
                          
                          Term2gene.df <- bind_rows(Treg.enriched.entrezid, Treg.reduced.entrezid)[,-1]
                          Term2gene.df <- data.frame(term = Term2gene.df$term,
                                                     gene = Term2gene.df$ENTREZID)                          
                  
                  treg.gsea <- GSEA(geneList = gsea_set,
                                    TERM2GENE = Term2gene.df,
                                    by = "fgsea",
                                    nPerm = 10000,
                                    pvalueCutoff = 1,
                                    minGSSize = 0,
                                    maxGSSize = 25000,
                                    
                                    exponent = 1,
                                    eps = 1e-10,
                                    pAdjustMethod = "BH", 
                                    TERM2NAME = NA,
                                    verbose = TRUE,
                                    seed = FALSE)
            
                  treg.gsea.result <- treg.gsea@result
                  treg.gsea.result <- treg.gsea.result %>% filter(ID == "Treg_enriched")

library(enrichplot)
                  
                  RES1 <- gseaplot2(treg.gsea, geneSetID = 1, color = c("black", "gray"), subplots = 1) + 
                                  geom_hline(yintercept = 0, width = 0.5) +
                                  scale_x_continuous(breaks = NULL) +
                                  scale_y_continuous(limits = c(-1,1)) +
                                  ggtitle("Treg Enriched Gene Set") +
                                  theme(plot.title = element_text(hjust = 0.5))
                  
                  RES2 <- gseaplot2(treg.gsea, geneSetID = 2, color = c("black", "gray"), subplots = 1) + 
                                  geom_hline(yintercept = 0, width = 0.5) +
                                  scale_x_continuous(breaks = NULL) +
                                  scale_y_continuous(limits = c(-1,1))  +
                                  ggtitle("Treg Reduced Gene Set") +
                                  theme(plot.title = element_text(hjust = 0.5))
                        
                  
                  RES <- RES1+RES2
                  
                  NES <- ggplot(treg.gsea.result, aes(x = ID, y = NES)) +
                                geom_col(width = 0.7) +
                                theme_classic() +
                                geom_hline(yintercept = 0, width = 0.5) +
                                scale_y_continuous(limits = c(-1.5,1.5)) +
                                labs(y = "Normalized Enrichment Score")
                  
                  #ggsave(RES, width = 40, height = 40,units = "cm",filename = "Plots/marker_heatmaps/Treg_GSEA_RES.pdf", limitsize = FALSE)  
                  #ggsave(NES, width = 40, height = 40,units = "cm",filename = "Plots/marker_heatmaps/Treg_GSEA_NES.pdf", limitsize = FALSE)  
                  


#Other plots:
  # KRT IFNG response GSEA:
              
              #Load KRT in vitro IFNG response DE list:
                KRT_invitro_IFNG_DE <- fread("gene_lists/KRT.IFNg.DE.tsv")  
                KRT_invitro_IFNG_DE <- KRT_invitro_IFNG_DE[,-1]
                KRT_invitro_IFNG_DE$padj <- as.numeric(KRT_invitro_IFNG_DE$padj)
                KRT_invitro_IFNG_DE <- KRT_invitro_IFNG_DE %>% filter(padj <= 0.01)
                KRT_invitro_IFNG_DE <- KRT_invitro_IFNG_DE %>% arrange(desc(log2FoldChange))
                
                KRT_invitro_IFNG_induced_genes <- KRT_invitro_IFNG_DE %>% filter(log2FoldChange >= 1.5) %>% dplyr::select(gene) %>% unlist() %>% as.character()
                KRT_invitro_IFNG_suppressed_genes <- KRT_invitro_IFNG_DE %>% filter(log2FoldChange <= -1.5) %>% dplyr::select(gene) %>% unlist() %>% as.character()
                
                      #Term2gene dataframe for custom Gene lists:
                      KRT_IFNG_induced <- bitr(KRT_invitro_IFNG_induced_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                      KRT_IFNG_induced$ENTREZID <- as.numeric(KRT_IFNG_induced$ENTREZID)
                      KRT_IFNG_induced$term <- "IFNG_induced"
                      
                      KRT_IFNG_suppressed <- bitr(KRT_invitro_IFNG_suppressed_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                      KRT_IFNG_suppressed$ENTREZID <- as.numeric(KRT_IFNG_suppressed$ENTREZID)
                      KRT_IFNG_suppressed$term <- "IFNG_suppressed"
                      
                      #Term2gene.df <- bind_rows(KRT_IFNG_induced, KRT_IFNG_suppressed)[,-1]
                      Term2gene.df <- KRT_IFNG_induced[,-1]
                      Term2gene.df <- data.frame(term = Term2gene.df$term,
                                                 gene = Term2gene.df$ENTREZID)      
              #DE:
                D2Allergy_vs_allcontrols <- fread("Data/DE_mincell0/Semi_Detailed_combo_lesions/KRT_Day2_Allergy_vs_Control.txt")
                D4Allergy_vs_allcontrols <- fread("Data/DE_mincell0/Semi_Detailed_combo_lesions/KRT_Day4_Allergy_vs_Control.txt")
                D2Allergy_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Day2_Allergy_vs_Nonlesional.txt")
                D4Allergy_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Day4_Allergy_vs_Nonlesional.txt")
                Irritant_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Irritant_vs_Nonlesional.txt")
                Acetone_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Acetone_Vehicle_vs_Nonlesional.txt")
                
                    colnames(D2Allergy_vs_allcontrols)[1] <- "Gene"
                    colnames(D4Allergy_vs_allcontrols)[1] <- "Gene"
                    colnames(D2Allergy_vs_NL)[1] <- "Gene"
                    colnames(D4Allergy_vs_NL)[1] <- "Gene"
                    colnames(Irritant_vs_NL)[1] <- "Gene"
                    colnames(Acetone_vs_NL)[1] <- "Gene"
                    
              #GSEA:
                  DE_comparison <- c("D2Allergy_vs_allcontrols", "D4Allergy_vs_allcontrols", "D2Allergy_vs_NL", "D4Allergy_vs_NL", "Irritant_vs_NL", "Acetone_vs_NL")
                  
                  for (c in 1:length(DE_comparison)) {
                    
                    comp <- DE_comparison[c]
                  
                            comp.DE.data <- get(comp)
                            comp.DE.data <- comp.DE.data %>% arrange(desc(logFC))
                            
                            temp.genelist <- comp.DE.data$logFC
                            names(temp.genelist) <- comp.DE.data$Gene
                            
                            universe_id <- bitr(names(temp.genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                            universe_id$ENTREZID <- as.numeric(universe_id$ENTREZID)
                            universe_id$avg_log2FC <- temp.genelist[match(universe_id$SYMBOL, names(temp.genelist))]
                            
                            gsea_set <- as.vector(universe_id$avg_log2FC)
                            names(gsea_set) <- universe_id$ENTREZID
                            
                            gsea_set <- sort(gsea_set, decreasing = TRUE)
                            
                            gsea <- GSEA(geneList = gsea_set,
                                              TERM2GENE = Term2gene.df,
                                              by = "fgsea",
                                              nPerm = 10000,
                                              pvalueCutoff = 1,
                                              minGSSize = 0,
                                              maxGSSize = 25000,
                                              
                                              exponent = 1,
                                              eps = 1e-10,
                                              pAdjustMethod = "BH", 
                                              TERM2NAME = NA,
                                              verbose = TRUE,
                                              seed = FALSE)
                    
                    assign(paste(comp, ".gsea", sep = ""),gsea)
                    rm(gsea)
                  }
                  
                  library(enrichplot)
                  
                  resplots = c()
                  nesplots = c()
                  for (c in 1:length(DE_comparison)) {
                    
                    comp <- DE_comparison[c]
                    
                    gsea <- get(paste(comp, ".gsea", sep = ""))
                    
                    RES <- gseaplot2(gsea, geneSetID = 1, color = c("black", "gray"), subplots = 1) + 
                      geom_hline(yintercept = 0, width = 0.5) +
                      scale_x_continuous(breaks = NULL) +
                      scale_y_continuous(limits = c(-1,1)) +
                      labs(title = paste(comp)) +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                    NES <- ggplot(gsea@result, aes(x = ID, y = NES)) +
                      geom_col(width = 0.7) +
                      theme_classic() +
                      geom_hline(yintercept = 0, width = 0.5) +
                      scale_y_continuous(limits = c(-2.5,2.5)) +
                      labs(y = "Normalized Enrichment Score", title = paste(comp)) +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                    
                    assign(paste(comp, ".res", sep = ""), RES)
                    assign(paste(comp, ".nes", sep = ""), NES)
                    
                    resplots <- c(resplots, paste(comp, ".res", sep = ""))
                    nesplots <- c(nesplots, paste(comp, ".nes", sep = ""))
                    
                    rm(gsea)
                  }  
                  
                  #Single lesions vs NL:
                  Irritant_vs_NL.res + Acetone_vs_NL.res + D2Allergy_vs_NL.res + D4Allergy_vs_NL.res + patchwork::plot_layout(ncol = 4, nrow = 1)


                  #Combo lesion comparisions:
                  D2Allergy_vs_allcontrols.res + D4Allergy_vs_allcontrols.res
                  
                  
  # KRT IL4/IL13 response GSEA:
                  
                  #Load KRT in vitro IL4 IL13 response DE list:
                  KRT_invitro_IL4_IL13_DE <- readxl::read_xlsx("gene_lists/ntert_KRTs_IL4_IL13_stim_DE.xlsx")  
                  
                  KRT_invitro_IL4_IL13_DE$Padj <- str_replace_all(KRT_invitro_IL4_IL13_DE$Padj, pattern = "E", replacement = "e")
                  KRT_invitro_IL4_IL13_DE$Padj <- str_replace_all(KRT_invitro_IL4_IL13_DE$Padj, pattern = "–", replacement = "-")
                  KRT_invitro_IL4_IL13_DE$log2FC <- str_replace_all(KRT_invitro_IL4_IL13_DE$log2FC, pattern = "–", replacement = "-")
                  
                  KRT_invitro_IL4_IL13_DE$log2FC <- as.numeric(KRT_invitro_IL4_IL13_DE$log2FC)
                  KRT_invitro_IL4_IL13_DE$Padj <- as.numeric(KRT_invitro_IL4_IL13_DE$Padj)
                  
                  KRT_invitro_IL4_IL13_DE <- KRT_invitro_IL4_IL13_DE %>% filter(Padj <= 0.01)
                  KRT_invitro_IL4_IL13_DE <- KRT_invitro_IL4_IL13_DE %>% arrange(desc(log2FC))
                  
                  KRT_invitro_IL4_IL13_induced_genes <- KRT_invitro_IL4_IL13_DE %>% filter(log2FC >= 1.5) %>% dplyr::select(Gene) %>% unlist() %>% as.character()
                  KRT_invitro_IL4_IL13_suppressed_genes <- KRT_invitro_IL4_IL13_DE %>% filter(log2FC <= -1.5) %>% dplyr::select(Gene) %>% unlist() %>% as.character()
                  
                  #Term2gene dataframe for custom Gene lists:
                  KRT_IL4_IL13_induced <- bitr(KRT_invitro_IL4_IL13_induced_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                  KRT_IL4_IL13_induced$ENTREZID <- as.numeric(KRT_IL4_IL13_induced$ENTREZID)
                  KRT_IL4_IL13_induced$term <- "IL4_IL13_induced"
                  
                  KRT_IL4_IL13_suppressed <- bitr(KRT_invitro_IL4_IL13_suppressed_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                  KRT_IL4_IL13_suppressed$ENTREZID <- as.numeric(KRT_IL4_IL13_suppressed$ENTREZID)
                  KRT_IL4_IL13_suppressed$term <- "IL4_IL13_suppressed"
                  
                  Term2gene.df <- bind_rows(KRT_IL4_IL13_induced, KRT_IL4_IL13_suppressed)[,-1]
                  Term2gene.df <- data.frame(term = Term2gene.df$term,
                                             gene = Term2gene.df$ENTREZID)      
                  #DE:
                  D2Allergy_vs_allcontrols <- fread("Data/DE_mincell0/Semi_Detailed_combo_lesions/KRT_Day2_Allergy_vs_Control.txt")
                  D4Allergy_vs_allcontrols <- fread("Data/DE_mincell0/Semi_Detailed_combo_lesions/KRT_Day4_Allergy_vs_Control.txt")
                  D2Allergy_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Day2_Allergy_vs_Nonlesional.txt")
                  D4Allergy_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Day4_Allergy_vs_Nonlesional.txt")
                  Irritant_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Irritant_vs_Nonlesional.txt")
                  Acetone_vs_NL <- fread("Data/DE_mincell0/Semi_Detailed_single_lesions/KRT_Acetone_Vehicle_vs_Nonlesional.txt")
                  
                  colnames(D2Allergy_vs_allcontrols)[1] <- "Gene"
                  colnames(D4Allergy_vs_allcontrols)[1] <- "Gene"
                  colnames(D2Allergy_vs_NL)[1] <- "Gene"
                  colnames(D4Allergy_vs_NL)[1] <- "Gene"
                  colnames(Irritant_vs_NL)[1] <- "Gene"
                  colnames(Acetone_vs_NL)[1] <- "Gene"
                  
                  #GSEA:
                  DE_comparison <- c("D2Allergy_vs_allcontrols", "D4Allergy_vs_allcontrols", "D2Allergy_vs_NL", "D4Allergy_vs_NL", "Irritant_vs_NL", "Acetone_vs_NL")
                  
                  for (c in 1:length(DE_comparison)) {
                    
                    comp <- DE_comparison[c]
                    
                    comp.DE.data <- get(comp)
                    comp.DE.data <- comp.DE.data %>% arrange(desc(logFC))
                    
                    temp.genelist <- comp.DE.data$logFC
                    names(temp.genelist) <- comp.DE.data$Gene
                    
                    universe_id <- bitr(names(temp.genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                    universe_id$ENTREZID <- as.numeric(universe_id$ENTREZID)
                    universe_id$avg_log2FC <- temp.genelist[match(universe_id$SYMBOL, names(temp.genelist))]
                    
                    gsea_set <- as.vector(universe_id$avg_log2FC)
                    names(gsea_set) <- universe_id$ENTREZID
                    
                    gsea_set <- sort(gsea_set, decreasing = TRUE)
                    
                    gsea <- GSEA(geneList = gsea_set,
                                 TERM2GENE = Term2gene.df,
                                 by = "fgsea",
                                 nPerm = 10000,
                                 pvalueCutoff = 1,
                                 minGSSize = 0,
                                 maxGSSize = 25000,
                                 
                                 exponent = 1,
                                 eps = 1e-10,
                                 pAdjustMethod = "BH", 
                                 TERM2NAME = NA,
                                 verbose = TRUE,
                                 seed = FALSE)
                    
                    assign(paste(comp, ".gsea", sep = ""),gsea)
                    rm(gsea)
                  }
                  
                  library(enrichplot)
                  
                  resplots = c()
                  nesplots = c()
                  for (c in 1:length(DE_comparison)) {
                    
                    comp <- DE_comparison[c]
                    
                    gsea <- get(paste(comp, ".gsea", sep = ""))
                    
                    RES <- gseaplot2(gsea, geneSetID = 1:2, color = c("black", "gray"), subplots = 1) + 
                      geom_hline(yintercept = 0, width = 0.5) +
                      scale_x_continuous(breaks = NULL) +
                      scale_y_continuous(limits = c(-1,1)) +
                      labs(title = paste(comp)) +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                    NES <- ggplot(gsea@result, aes(x = ID, y = NES)) +
                      geom_col(width = 0.7) +
                      theme_classic() +
                      geom_hline(yintercept = 0, width = 0.5) +
                      scale_y_continuous(limits = c(-2.5,2.5)) +
                      labs(y = "Normalized Enrichment Score", title = paste(comp)) +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                    
                    assign(paste(comp, ".res", sep = ""), RES)
                    assign(paste(comp, ".nes", sep = ""), NES)
                    
                    resplots <- c(resplots, paste(comp, ".res", sep = ""))
                    nesplots <- c(nesplots, paste(comp, ".nes", sep = ""))
                    
                    rm(gsea)
                  }  
                  
                  #Single lesions vs NL:
                  Irritant_vs_NL.nes + Acetone_vs_NL.nes + D2Allergy_vs_NL.nes + D4Allergy_vs_NL.nes + 
                    Irritant_vs_NL.res + Acetone_vs_NL.res + D2Allergy_vs_NL.res + D4Allergy_vs_NL.res + patchwork::plot_layout(ncol = 4, nrow = 2)
                  
                  
                  #Combo lesion comparisions:
                  D2Allergy_vs_allcontrols.nes + D4Allergy_vs_allcontrols.nes + D2Allergy_vs_allcontrols.res + D4Allergy_vs_allcontrols.res
                  
                                  
                                
    
  plot_tsne_metadata_srt(contact_derm, color_by = "Basic_Celltype", size = 0.01, plot_label = FALSE)
  plot_tsne_metadata_srt(contact_derm, color_by = "Detailed_Celltype", size = 0.01, plot_label = TRUE)
  
  plot_violin_srt(contact_derm, gene = "IL4", color_by = "Lesion", facet_by = "Semi_Detailed_Celltype", )

  
  #Build Custom Color Palette:
  flat.palette <- c("#F44336", "#C2185B", "#BA68C8", "#673AB7", "#7986CB", "#1E88E5", "#26C6DA", "#00897B", "#4CAF50",  "#AFB42B", "#FF8F00", "#FF5722", "#8D6E63", "#78909C")
  GetPalette = colorRampPalette(flat.palette)
  
  ColorCount = as.numeric(length(unique(contact_derm@meta.data$Patient)))
  discrete.colors1 <- sample(GetPalette(ColorCount), size = ColorCount, replace = FALSE)
  
  plot_tsne_metadata_srt(contact_derm, color_by = "Patient", facet_by = c("Disease", "Lesion"), colors = discrete.colors1, facet_grid = TRUE, plot_legend = FALSE, plot_labels = FALSE)
  
