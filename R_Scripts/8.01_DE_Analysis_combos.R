#######################################################
# DE Analysis: Counting DE genes by sample permutation
#######################################################
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
library(RColorBrewer)
library (purrr)
library(data.table)
library(readxl)
library(clusterProfiler)

edgeRDE_srt <- function (input, groups, batch, sizefactor, lib_size, minCells, 
                         pVal = 1, contrast = list(c(-1, 1))) {  
  
  rvar <- apply(input@assays$RNA@counts, 1, var)
  idx.keep <- which(rvar > 0)
  input <- input[idx.keep, ]
  
  gfrac <- apply(input@assays$RNA@counts, 1, function(a) length(a[a > 
                                                                    0])/length(a))
  idx.keep <- which(gfrac > minCells)
  input <- input[idx.keep, ]
  
  y <- edgeR::DGEList(counts = input@assays$RNA@counts, group = groups)
  
  y$samples$lib.size <- lib_size
  
  if (is.null(sizefactor)) {
    message(sprintf("Calculating norm factors..."))
    flush.console()
    y <- edgeR::calcNormFactors(y)
  }
  else {
    y$samples$norm.factors = input@meta.data$sizefactor/lib_size
  }
  if (length(unique(batch)) > 1) {
    design <- model.matrix(~0 + batch + group, data = y$samples)
  }
  else {
    design <- model.matrix(~0 + group, data = y$samples)
  }
  message("Estimating dispersion...")
  flush.console()
  y <- edgeR::estimateDisp(y, design)
  message("Doing likelihood ratio fit...")
  flush.console()
  
  if (is.na(y$common.dispersion) == TRUE) {
    tab <- NA
  }
  if (is.na(y$common.dispersion) == FALSE){
    fit <- edgeR::glmFit(y, design)
    tab <- list()
    if (length(unique(batch)) > 1) {
      lrt <- edgeR::glmLRT(fit)
      tab[["contrast_1"]] <- as.data.frame(edgeR::topTags(lrt, 
                                                          p.value = pVal, n = Inf, sort.by = "logFC"))
    }
    else {
      for (j in 1:length(contrast)) {
        cStr <- sprintf("contrast_%d", j)
        lrt <- edgeR::glmLRT(fit, contrast = contrast[[j]])
        tab[[cStr]] <- as.data.frame(edgeR::topTags(lrt, 
                                                    p.value = pVal, n = Inf, sort.by = "logFC"))
      }
    }
    return(tab)
  }
}

findDEgenes_srt <- function (input, DEgroup, contrastID, sizefactor = NULL, 
                             lib_size = NULL, facet_by, minCells = 0, batchID = NULL, 
                             outdir = "DEresults/", outsuffix = "DEresults.tsv", pVal = 1, 
                             contrast = list(c(-1, 1))) {
  
  if (is.null(lib_size)) {
    libsize = colSums(input@assays$RNA@counts)
    input@meta.data$lib_size = libsize
    lib_size = "lib_size"
  }
  for (i in 1:length(unique(input@meta.data[, facet_by]))) {
    name = sort(unique(as.character(input@meta.data[, facet_by])))[i]
    print(paste0("Performing DE for ", name))
    
    idx = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name)]
    cntr = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name & input@meta.data[, 
                                                                                                 DEgroup] == contrastID)]
    z = as.matrix(input@assays$RNA@counts)[, idx]
    if (!is.null(dim(z))) {
      if (ncol(z) > 2) {
        if (is.null(batchID)) {
          batch = as.factor(rep(1, ncol(z)))
        }
        else {
          batch_table = as.matrix(table(input@meta.data[idx, batchID], 
                                        input@meta.data[idx, DEgroup]))
          if (length(batch_table[batch_table == 0]) < 2) {
            batch = as.factor(input@meta.data[idx, batchID])
          }
          else {
            message(sprintf("Warning: not enough cells to include batch in model, setting batches to 1"))
            flush.console()
            batch = as.factor(rep(1, ncol(z)))
          }
        }
        groupList = rep(0, times = ncol(z))
        groupList[which(colnames(z) %in% cntr)] = 1
        if (length(unique(groupList)) > 1) {
          group = factor(groupList)
          
          Idents(input) <- facet_by.var
          z <- subset(input, idents = name)
          z_libsize <- z@meta.data$lib_size
          
          tab = edgeRDE_srt(input = z, groups = group, batch = batch, 
                            sizefactor = sizefactor, lib_size = z_libsize, 
                            minCells = minCells, pVal = pVal, contrast = contrast)
          
          DEtablecntr = tab[["contrast_1"]]
          outfile = paste(name, contrastID, outsuffix, 
                          sep = "_")
          write.table(DEtablecntr, paste(outdir, outfile, 
                                         sep = ""), sep = "\t")
        }
      }
    }
  }
}

#edited for batch:
findDEgenes_srt <- function (input, DEgroup, contrastID, sizefactor = NULL, 
                             lib_size = NULL, facet_by, minCells = 0, batchID = NULL, 
                             outdir = "DEresults/", outsuffix = "DEresults.tsv", pVal = 1, 
                             contrast = list(c(-1, 1))) {
  
  if (is.null(lib_size)) {
    libsize = colSums(input@assays$RNA@counts)
    input@meta.data$lib_size = libsize
    lib_size = "lib_size"
  }
  for (i in 1:length(unique(input@meta.data[, facet_by]))) {
    name = sort(unique(as.character(input@meta.data[, facet_by])))[i]
    print(paste0("Performing DE for ", name))
    
    idx = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name)]
    cntr = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name & input@meta.data[, 
                                                                                                 DEgroup] == contrastID)]
    
    z = as.matrix(input@assays$RNA@counts)[, idx]
    if (!is.null(dim(z))) {
      if (ncol(z) > 2) {
        if (is.null(batchID)) {
          batch = as.factor(rep(1, ncol(z)))
        }
        else {
          batch_table = as.matrix(table(input@meta.data[idx, batchID], 
                                        input@meta.data[idx, DEgroup]))
          #Test full rank:
          groupList = rep(0, times = ncol(z))
          groupList[which(colnames(z) %in% cntr)] = 1
          if (length(unique(groupList)) > 1) {
            group = factor(groupList)
          }
          y <- edgeR::DGEList(counts = as.matrix(z), group = group)
          batch = as.factor(input@meta.data[idx, batchID])
          
          design <- model.matrix(~0 + batch + group, data = y$samples)
          
          if (qr(design)$rank == ncol(design)) {
            batch = as.factor(input@meta.data[idx, batchID])
          }
          else {
            message(sprintf("Warning: not enough cells and/or design matrix not full rank to include batch in model, setting batches to 1"))
            flush.console()
            batch = as.factor(rep(1, ncol(z)))
          }
        }
        groupList = rep(0, times = ncol(z))
        groupList[which(colnames(z) %in% cntr)] = 1
        if (length(unique(groupList)) > 1) {
          group = factor(groupList)
          
          Idents(input) <- facet_by.var
          z <- subset(input, idents = name)
          #Note to self: shouldn't celltype subsetting be done first at the top, before batch design matrix rank calculation?
          
          z_libsize <- z@meta.data$lib_size
          
          tab = edgeRDE_srt(input = z, groups = group, batch = batch, 
                            sizefactor = sizefactor, lib_size = z_libsize, 
                            minCells = minCells, pVal = pVal, contrast = contrast)
          
          DEtablecntr = tab[["contrast_1"]]
          outfile = paste(name, contrastID, outsuffix, 
                          sep = "_")
          write.table(DEtablecntr, paste(outdir, outfile, 
                                         sep = ""), sep = "\t")
        }
      }
    }
  }
}
countDEgenes_srt <- function (input, DEgroup, contrastID, sizefactor = NULL, 
                              lib_size = NULL, facet_by, minCells = 0.1, batchID = NULL, 
                              outdir = "DEresults/", outsuffix = "DEresults.tsv", pVal = 1, 
                              contrast = list(c(-1, 1))) {
  
  if (is.null(lib_size)) {
    libsize = colSums(input@assays$RNA@counts)
    input@meta.data$lib_size = libsize
    lib_size = "lib_size"
  }
  for (i in 1:length(unique(input@meta.data[, facet_by]))) {
    name = sort(unique(as.character(input@meta.data[, facet_by])))[i]
    print(paste0("Performing DE for ", name))
    
    idx = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name)]
    cntr = rownames(input@meta.data)[which(input@meta.data[, facet_by] == name & input@meta.data[, 
                                                                                                 DEgroup] == contrastID)]
    z = as.matrix(input@assays$RNA@counts)[, idx]
    if (ncol(z) > 2) {
      if (is.null(batchID)) {
        batch = as.factor(rep(1, ncol(z)))
      }
      else {
        batch_table = as.matrix(table(input@meta.data[idx, batchID], 
                                      input@meta.data[idx, DEgroup]))
        if (length(batch_table[batch_table == 0]) < 2) {
          batch = as.factor(input@meta.data[idx, batchID])
        }
        else {
          message(sprintf("Warning: not enough cells to include batch in model, setting batches to 1"))
          flush.console()
          batch = as.factor(rep(1, ncol(z)))
        }
      }
      groupList = rep(0, times = ncol(z))
      groupList[which(colnames(z) %in% cntr)] = 1
      if (length(unique(groupList)) > 1) {
        group = factor(groupList)
        
        Idents(input) <- facet_by.var
        z <- subset(input, idents = name)
        z_libsize <- z@meta.data$lib_size
        
        tab = edgeRDE_srt(input = z, groups = group, batch = batch, 
                          sizefactor = sizefactor, lib_size = z_libsize, 
                          minCells = minCells, pVal = pVal, contrast = contrast)
        
        ##Same as findDEgenes_srt function above. Below here counts genes as numbers:
        if (tab == "NA") {
          DE.gene.count <- NA
          output <- data.frame(facet_by = name,
                               DE.gene.count = DE.gene.count)
          output$DE.gene.count <- as.character(DE.gene.count)
          assign(paste(name, ".output", sep = ""), output)
        }
        if (tab != "NA") {
          DEtablecntr = tab[["contrast_1"]]
          DE.gene.count <- as.numeric(nrow(DEtablecntr))
          output <- data.frame(facet_by = name,
                               DE.gene.count = DE.gene.count)
          output$DE.gene.count <- as.character(DE.gene.count)
          assign(paste(name, ".output", sep = ""), output)   
        }
      }
    }
  }
}

setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

#Load data:
load(file = "Data/contact_derm_filtered3.Rdata")

DE.input = contact_derm

#Only Run DE on samples with visual reactoin scores of 1 or greater:
        DE.input@meta.data <- DE.input@meta.data %>% mutate(DE_inclusion_factor = case_when(Lesion == "Nonlesional" ~ "yes",
                                                                                            Lesion == "Acetone_Vehicle" & Lesion_visual_score < 1 ~ "yes",
                                                                                            Lesion == "Irritant" ~ "yes",
                                                                                            Lesion == "Day2_Allergy" ~ "yes",
                                                                                            Lesion == "Day4_Allergy" ~ "yes",
                                                                                            TRUE ~ "no"))

        
        
        Idents(DE.input) <- "DE_inclusion_factor"
        DE.input <- subset(DE.input, idents = "yes")

#Lesion combos for DE:
        DE.input@meta.data$Lesion <- as.character(DE.input@meta.data$Lesion)
        DE.input@meta.data <- DE.input@meta.data %>% mutate(Lesion_combos = case_when(Lesion %in% c("Nonlesional", "Irritant", "Acetone_Vehicle") ~ "Control",
                                                                                      Lesion == "Day2_Allergy" ~ "Day2_Allergy",
                                                                                      Lesion == "Day4_Allergy" ~ "Day4_Allergy"))
        DE.input@meta.data$Lesion_combos <- factor(DE.input@meta.data$Lesion_combos,
                                                   levels = c("Control", "Day2_Allergy", "Day4_Allergy"))
        
        DE.input@meta.data <- DE.input@meta.data %>% mutate(Lesion_combos2 = case_when(Lesion %in% c("Nonlesional") ~ "Nonlesional",
                                                                                      Lesion == "Irritant" ~ "Inflammatory_Reference",
                                                                                      Lesion == "Day2_Allergy" ~ "Inflammatory_Reference",
                                                                                      Lesion == "Day4_Allergy" ~ "Inflammatory_Reference"))
        
        DE.input@meta.data$Lesion_combos2 <- factor(DE.input@meta.data$Lesion_combos2,
                                                   levels = c("Inflammatory_Reference", "Nonlesional"))

        #Batch groups for DE:
              #Grouping Sequencers into 2 groups based on the largest differences: NS500602, NB501205, A00197 vs VH00230 and A00439
              DE.input@meta.data$Sequence_instrument <- as.character(DE.input@meta.data$Sequence_instrument)
              DE.input@meta.data <- DE.input@meta.data %>% mutate(BatchDE = case_when(Sequence_instrument %in% c("NS500602", "NB501205", "A00197") ~ "Batch1",
                                                                                      Sequence_instrument %in% c("VH00230", "A00439") ~ "Batch2"))
              DE.input@meta.data$BatchDE <- factor(DE.input@meta.data$BatchDE,
                                                   levels = c("Batch1", "Batch2"))
      



#Set Comparisons:
            #comparisons <- as.character(sort(unique(DE.input@meta.data$Lesion)))
            #comparisons <- expand.grid(comparisons, comparisons, stringsAsFactors = FALSE)
            #colnames(comparisons) <- c("reference", "comparison")
            #comparisons <- comparisons[c(5,10,15,20,9,14,6,11),]

#For Legion combo comparisons:
      #comparisons <- data.frame(reference = c("Control", "Control"), comparison = c("Day2_Allergy","Day4_Allergy"))
      
      comparisons <- data.frame(reference = "Inflammatory_Reference",
                                comparison = "Nonlesional")    
#Compute DE
DE.var = "Lesion_combos2"
facet_by.var = "Semi_Detailed_Celltype"
batch.var = NULL

Idents(DE.input) <- DE.var

for (i in 1:nrow(comparisons)) {
  int <- comparisons[i,]
  reference <- int[,1]
  compare <- int[,2]
  print(paste0(reference, " vs ", compare))
  
  DE_subset <- subset(DE.input, idents = c(unique(as.character(Idents(DE.input)))[which(unique(as.character(Idents(DE.input))) %in% as.character(int))]))
  
  findDEgenes_srt(input = DE_subset,
                  contrastID = compare,
                  DEgroup = DE.var,
                  facet_by = facet_by.var,
                  batchID = batch.var,
                  sizefactor = NULL,
                  lib_size = NULL,
                  minCells = 0,
                  outsuffix = paste0("vs_",reference,".txt"),
                  outdir = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/DE_less_acetone_rxns/Semi_Detailed_combo_lesions/",
                  contrast = list(c(-1,1)))
  
}






