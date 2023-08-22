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
library(DESeq2)
library(edgeR)


setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/")


#Calc_agg_bulk:
calc_agg_bulk_srt <- function (input, aggregate_by, group_by = FALSE, cutoff_frac = FALSE, 
                               cutoff_num = FALSE, cutoff_cpm = FALSE, expand_with_dummies = FALSE) {
  
  if (class(input) == "ExpressionSet"){
    
    if (group_by != FALSE) {
      ind <- match(group_by, aggregate_by)
      if (ind != 1) {
        stop("Please provide the group_by value first in the aggreggate_by argument")
      }
    }
    
    #Hide previous aggbulk columns from this analysis:
    check <- grep("bulk", colnames(fData(input)))
    if (length(check) > 0) {
      fData(input) <- fData(input)[, -check]
    }
    
    #Temporarily Convert pData Factors to Characters:
    if (any(sapply(pData(input), FUN = class) == "factor")) {
      
      factor.columns <- which(sapply(pData(input), FUN = class) == "factor")
      
      #Store factor columns for later:
      meta.data.factors <- pData(input)[,factor.columns]
      
      pData(input)[,factor.columns] <- sapply(pData(input)[,factor.columns], FUN = as.character)
    }
    
    #Setup to compute bulks:
    to_expand <- vector("list", length(aggregate_by))
    for (i in 1:length(aggregate_by)) {
      var <- aggregate_by[i]
      vars <- unique(pData(input)[, var])
      vars <- sort(vars)
      to_expand[[i]] <- vars
    }
    names(to_expand) <- aggregate_by
    bulks <- expand.grid(to_expand, stringsAsFactors = FALSE)
    colnames(bulks) <- c(aggregate_by)
    groups <- sort(unique(pData(input)[, group_by]))
    upm_vals <- c()
    num_cells_vals <- c()
    num_genes_vals <- c()
    rem_genes <- vector(mode = "list", length = nrow(bulks))
    
    #Compute Bulks:
    for (j in 1:nrow(bulks)) {
      int <- bulks[j, ]
      full_match <- c()
      for (k in 1:length(int)) {
        ind <- which(pData(input)[, colnames(bulks)[k]] == 
                       int[[k]])
        if (k == 1) {
          full_match <- c(full_match, ind)
        }
        else {
          full_match <- intersect(full_match, ind)
        }
      }
      num_cells <- length(full_match)
      num_cells_vals <- c(num_cells_vals, num_cells)
      #mike edit: changed length(full_match) to > "0" instead of "1":
      if (length(full_match) > 0) {
        #mike edit: added as.matrix to tmp line:
        tmp <- as.matrix(exprs(input)[, full_match])
        #mike edit: added if all(tmp == 0) statements for dummy cell treatments:
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == TRUE) {
          upm <- tmp
        }
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == FALSE) {
          upm <- (apply(tmp, 1, sum)/sum(tmp)) * 1e+06
        }
        if (cutoff_frac != FALSE || cutoff_num != FALSE) {
          zero_out_frac <- c()
          zero_out_num <- c()
          tmp2 <- tmp
          tmp2[which(tmp2 > 0)] <- 1
          gSums <- apply(tmp2, 1, sum)
          if (cutoff_frac != FALSE) {
            frac <- gSums/num_cells
            zero_out_frac <- which(frac < cutoff_frac)
          }
          if (cutoff_num != FALSE) {
            zero_out_num <- which(gSums < cutoff_num)
          }
          if (group_by != FALSE) {
            rem_genes[[j]] <- unique(names(c(zero_out_num, 
                                             zero_out_frac)))
          }
          else {
            upm[unique(c(zero_out_num, zero_out_frac))] <- 0
          }
        }
      }
      #mike edit: add dummy values for bulk permutations that don't exist in data:
      if (length(full_match) == 0) {
        tmp <- rep(0, times = nrow(exprs(input)))
        names(tmp) <- rownames(exprs(input))
        upm <- tmp
      }
      expressed <- length(which(upm > 0))
      upm_vals <- c(upm_vals, upm)
      num_genes_vals <- c(num_genes_vals, expressed)
    }
    
    names(rem_genes) <- apply(bulks, 1, FUN = paste, collapse = "__")
    bulks$numcells <- num_cells_vals
    bulks$numgenes <- num_genes_vals
    bulks$proportion <- 0
    if (group_by != FALSE) {
      for (i in 1:length(groups)) {
        ind <- grep(groups[i], bulks[, group_by])
        total <- sum(bulks$numcells[ind])
        bulks$proportion[ind] <- round((bulks$numcells[ind]/total) * 
                                         100, 2)
      }
    } else {
      total <- sum(bulks$numcells)
      bulks$proportion <- round((bulks$numcells/total) * 100, 
                                2)
    }
    bulk <- matrix(upm_vals, nrow = nrow(exprs(input)))
    rownames(bulk) <- rownames(exprs(input))
    colnames(bulk) <- seq(1:ncol(bulk))
    for (l in 1:nrow(bulks)) {
      cname <- bulks[l, -c(match(c("numcells", "proportion", 
                                   "numgenes"), colnames(bulks)))]
      cname2 <- c()
      for (i in 1:length(cname)) {
        cint <- as.character(cname[[i]])
        cname2 <- c(cname2, cint)
      }
      cname <- cname2
      cnum <- bulks[l, "numcells"]
      cpro <- bulks[l, "proportion"]
      cgen <- bulks[l, "numgenes"]
      cname <- paste0(c(cname, "num__genes", cgen, "num__cells", 
                        cnum, "percent", cpro, "bulk"), collapse = "__")
      colnames(bulk)[l] <- cname
    }
    if (group_by != FALSE) {
      if (!is.null(unlist(rem_genes))) {
        for (i in 1:length(vars)) {
          int_cell <- vars[i]
          ind <- grep(int_cell, names(rem_genes))
          vals <- table(unlist(rem_genes[ind]))
          zero_out <- names(which(vals == max(vals)))
          bulk[zero_out, ind] <- 0
        }
      }
    }
    if (cutoff_cpm != FALSE) {
      for (i in 1:length(vars)) {
        int_cell <- vars[i]
        ind <- grep(int_cell, colnames(bulk))
        if (length(ind) > 1) {
          gCount <- apply(bulk[, ind], 1, function(x) length(which(x >= cutoff_cpm)))
          zero_out <- which(gCount == 0)
        }
        else {
          zero_out <- which(bulk[, ind] < cutoff_cpm)
        }
        bulk[zero_out, ind] <- 0
      }
    }
    
    #Return pData Factor columns:
    if (("meta.data.factors" %in% ls()) == TRUE) {
      pData(input)[,factor.columns] <- meta.data.factors
    }
    
    #Remove non-existent dummy combinations:
    if (expand_with_dummies == FALSE) {
      
      no_cell_combos <- which(str_detect(colnames(bulk), pattern = "__num__cells__0__"))
      
      if (length(no_cell_combos) > 0){
        bulk <- bulk[,-no_cell_combos]  
      }
    }
    
    
    #Bind new agg bulk data:
    fData(input) <- cbind(fData(input), bulk)
  }
  
  if (class(input) == "Seurat") {
    
    if (group_by != FALSE) {
      ind <- match(group_by, aggregate_by)
      if (ind != 1) {
        stop("Please provide the group_by value first in the aggreggate_by argument")
      }
    }
    
    #Hide previous aggbulk columns from this analysis:
    check <- grep("bulk", colnames(input@assays$RNA@meta.features))
    if (length(check) > 0) {
      input@assays$RNA@meta.features <- input@assays$RNA@meta.features[, -check]
    }
    
    #Temporarily Convert Metadata Factors to Characters:
    if (any(sapply(input@meta.data, FUN = class) == "factor")) {
      
      factor.columns <- which(sapply(input@meta.data, FUN = class) == "factor")
      
      #Store factor columns for later:
      meta.data.factors <- input@meta.data[,factor.columns]
      
      input@meta.data[,factor.columns] <- sapply(input@meta.data[,factor.columns], FUN = as.character)
    }
    
    #Setup to compute bulks:
    to_expand <- vector("list", length(aggregate_by))
    for (i in 1:length(aggregate_by)) {
      var <- aggregate_by[i]
      vars <- unique(input@meta.data[, var])
      vars <- sort(vars)
      to_expand[[i]] <- vars
    }
    names(to_expand) <- aggregate_by
    bulks <- expand.grid(to_expand, stringsAsFactors = FALSE)
    colnames(bulks) <- c(aggregate_by)
    groups <- sort(unique(input@meta.data[, group_by]))
    upm_vals <- c()
    num_cells_vals <- c()
    num_genes_vals <- c()
    rem_genes <- vector(mode = "list", length = nrow(bulks))
    
    #Compute Bulks:
    for (j in 1:nrow(bulks)) {
      int <- bulks[j, ]
      full_match <- c()
      for (k in 1:length(int)) {
        ind <- which(input@meta.data[, colnames(bulks)[k]] == 
                       int[[k]])
        if (k == 1) {
          full_match <- c(full_match, ind)
        }
        else {
          full_match <- intersect(full_match, ind)
        }
      }
      num_cells <- length(full_match)
      num_cells_vals <- c(num_cells_vals, num_cells)
      
      #mike edit: changed length(full_match) to > "0" instead of "1":
      if (length(full_match) > 0) {
        #mike edit: added as.matrix to tmp line:
        tmp <- as.matrix(input@assays$RNA@data[, full_match])
        #mike edit: added if all(tmp == 0) statements for dummy cell treatments:
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == TRUE) {
          upm <- tmp
        }
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == FALSE) {
          upm <- (apply(tmp, 1, sum)/sum(tmp)) * 1e+06
        }
        if (cutoff_frac != FALSE || cutoff_num != FALSE) {
          zero_out_frac <- c()
          zero_out_num <- c()
          tmp2 <- tmp
          tmp2[which(tmp2 > 0)] <- 1
          gSums <- apply(tmp2, 1, sum)
          if (cutoff_frac != FALSE) {
            frac <- gSums/num_cells
            zero_out_frac <- which(frac < cutoff_frac)
          }
          if (cutoff_num != FALSE) {
            zero_out_num <- which(gSums < cutoff_num)
          }
          if (group_by != FALSE) {
            rem_genes[[j]] <- unique(names(c(zero_out_num, 
                                             zero_out_frac)))
          }
          else {
            upm[unique(c(zero_out_num, zero_out_frac))] <- 0
          }
        }
      }
      #mike edit: add dummy values for bulk permutations that don't exist in data:
      if (length(full_match) == 0) {
        tmp <- rep(0, times = nrow(input@assays$RNA@data))
        names(tmp) <- rownames(input@assays$RNA@data)
        upm <- tmp
      }
      expressed <- length(which(upm > 0))
      upm_vals <- c(upm_vals, upm)
      num_genes_vals <- c(num_genes_vals, expressed)
    }
    
    names(rem_genes) <- apply(bulks, 1, FUN = paste, collapse = "__")
    bulks$numcells <- num_cells_vals
    bulks$numgenes <- num_genes_vals
    bulks$proportion <- 0
    if (group_by != FALSE) {
      for (i in 1:length(groups)) {
        ind <- grep(groups[i], bulks[, group_by])
        total <- sum(bulks$numcells[ind])
        bulks$proportion[ind] <- round((bulks$numcells[ind]/total) * 
                                         100, 2)
      }
    } else {
      total <- sum(bulks$numcells)
      bulks$proportion <- round((bulks$numcells/total) * 100, 
                                2)
    }
    bulk <- matrix(upm_vals, nrow = nrow(input@assays$RNA@data))
    rownames(bulk) <- rownames(input@assays$RNA@data)
    colnames(bulk) <- seq(1:ncol(bulk))
    for (l in 1:nrow(bulks)) {
      cname <- bulks[l, -c(match(c("numcells", "proportion", 
                                   "numgenes"), colnames(bulks)))]
      cname2 <- c()
      for (i in 1:length(cname)) {
        cint <- as.character(cname[[i]])
        cname2 <- c(cname2, cint)
      }
      cname <- cname2
      cnum <- bulks[l, "numcells"]
      cpro <- bulks[l, "proportion"]
      cgen <- bulks[l, "numgenes"]
      cname <- paste0(c(cname, "num__genes", cgen, "num__cells", 
                        cnum, "percent", cpro, "bulk"), collapse = "__")
      colnames(bulk)[l] <- cname
    }
    if (group_by != FALSE) {
      if (!is.null(unlist(rem_genes))) {
        for (i in 1:length(vars)) {
          int_cell <- vars[i]
          ind <- grep(int_cell, names(rem_genes))
          vals <- table(unlist(rem_genes[ind]))
          zero_out <- names(which(vals == max(vals)))
          bulk[zero_out, ind] <- 0
        }
      }
    }
    if (cutoff_cpm != FALSE) {
      for (i in 1:length(vars)) {
        int_cell <- vars[i]
        ind <- grep(int_cell, colnames(bulk))
        if (length(ind) > 1) {
          gCount <- apply(bulk[, ind], 1, function(x) length(which(x >= cutoff_cpm)))
          zero_out <- which(gCount == 0)
        }
        else {
          zero_out <- which(bulk[, ind] < cutoff_cpm)
        }
        bulk[zero_out, ind] <- 0
      }
    }
    
    #Return Metadata Factor Columns:
    if (("meta.data.factors" %in% ls()) == TRUE) {
      input@meta.data[,factor.columns] <- meta.data.factors
    }
    
    #Remove non-existent dummy combinations:
    if (expand_with_dummies == FALSE) {
      
      no_cell_combos <- which(str_detect(colnames(bulk), pattern = "__num__cells__0__"))
      
      if (length(no_cell_combos) > 0){
        bulk <- bulk[,-no_cell_combos]  
      }
      
    }
    
    #Bind new agg bulk data:
    input@assays$RNA@meta.features<- cbind(input@assays$RNA@meta.features, bulk)
    
    
  }
  
  return(input)
  
}



#Load data:
load(file = "Data/contact_derm_filtered.Rdata")

DE.input = contact_derm

#Set Comparisons:
comparisons <- as.character(sort(unique(DE.input@meta.data$Lesion)))
comparisons <- expand.grid(comparisons, comparisons, stringsAsFactors = FALSE)
colnames(comparisons) <- c("reference", "comparison")
comparisons <- comparisons[c(11, 16, 21),]



#Set Variables and outdir:
DE.var = "Lesion"
facet_by.var = "Semi_Detailed_Celltype"
Batch = NULL

outdir = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_Bulk_DE/"

Idents(DE.input) <- DE.var

      # Bulk DE by EdgeR
      DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c(facet_by.var, DE.var, "Patient"))
      
      Bulk_data <- DE.input@assays$RNA@meta.features[,which(str_detect(colnames(DE.input@assays$RNA@meta.features), pattern = "bulk"))]
      
      facets <- as.character(unique(DE.input@meta.data[,facet_by.var]))
      
      reference = unique(comparisons$reference)
      
      #Filter Bulk Data by cell numbers:
      Bulk_data <- as.data.frame(t(Bulk_data))
      
      Bulk_data$cell.num <- as.numeric(sapply(rownames(Bulk_data), FUN = function(x) str_split(x, pattern = "__")[[1]][which(str_split(x, pattern = "__")[[1]] == "cells")+1]))
      
      bulk_gene_count_df <- Bulk_data[,1:ncol(Bulk_data)-1]
      bulk_gene_count_df[bulk_gene_count_df>0] <- 1
      Bulk_data$gene.sum <- apply(bulk_gene_count_df, MARGIN = 1, FUN = sum)
      
      Bulk_data <- Bulk_data %>% filter(cell.num >= 5)
      
      Bulk_data <- as.data.frame(t(Bulk_data[,-which(colnames(Bulk_data) %in% c("cell.num", "gene.sum"))]))
      
  for (c in 1:length(facets)) {
        
        facet <- facets[c]
        
        # Subset by facet:
        facet_string <- as.character(sapply(colnames(Bulk_data), FUN = function(x) str_split(x, pattern = "__")[[1]][which(str_split(x, pattern = "__")[[1]] %in% facets)]))
        facet_subset <- Bulk_data[,which(facet_string == facet)] 
        
        #Create Metadata Dataframe:
        DE.vars <- as.character(unique(DE.input@meta.data[,DE.var]))
        DE.var_string <- as.character(sapply(colnames(Bulk_data), FUN = function(x) str_split(x, pattern = "__")[[1]][which(str_split(x, pattern = "__")[[1]] %in% DE.vars)]))
        
        Patients <- as.character(unique(DE.input@meta.data$Patient))
        Patient_string <- as.character(sapply(colnames(Bulk_data), FUN = function(x) str_split(x, pattern = "__")[[1]][which(str_split(x, pattern = "__")[[1]] %in% Patients)]))
        
        metadata <- data.frame(facet = facet_string[which(facet_string == facet)],
                               DE.var = DE.var_string[which(facet_string == facet)],
                               Patient = Patient_string[which(facet_string == facet)])
        colnames(metadata) <- c(facet_by.var, DE.var, "Patient")
        rownames(metadata) <- colnames(facet_subset)
        
        if (levels(DE.input@meta.data[,DE.var])[1] != reference) {
          stop("set levels of DE.var so that reference condition is position 1")
        }
        
        metadata[,DE.var] <- factor(metadata[,DE.var], 
                                    levels = levels(DE.input@meta.data[,DE.var]))
        
        #Run EdgeR:
          #Filter Genes First:
              #filter genes with zero variance:
              gene.vars <- apply(facet_subset, 1, var)
              facet_subset <- facet_subset[which(gene.vars > 0), ]
              
              gene.sums <- apply(facet_subset, 1, sum) # total gene counts per sample
              facet_subset <- facet_subset[which(gene.sums >= 10), ]

          #Create edgeR list object:
          d <- DGEList(counts=facet_subset, genes=rownames(facet_subset), group = metadata[,DE.var])
        
          #Normalize:
          d$samples$lib.size <- colSums(d$counts)
          
          d <- calcNormFactors(d)
          
          if (!is.null(Batch)) {
            design.data <- data.frame(group = metadata[,DE.var],
                                      batch = metadata[,Batch])
            
            design <- model.matrix(~ 0 + batch + group, data=design.data)
            
            colnames(design) <- str_replace(colnames(design), pattern = "group", replacement = "")
          }

          if (is.null(Batch)) {
            design <- model.matrix(~0 + group, data=d$samples)
            colnames(design) <- levels(d$samples$group)
          }
          
          d <- edgeR::estimateDisp(d, design)
          
          fit <- glmFit(d, design)
          
          DE.conditions <- DE.vars[which(DE.vars %in% comparisons$comparison)]
          DE.conditions <- DE.conditions[which(DE.conditions %in% metadata$Lesion)]
          
            for (condition in DE.conditions) {
              
              pair_vector = sprintf("%s-%s", condition, reference) # Samples to be compared
              pair_contrast = makeContrasts(contrasts=pair_vector, levels=design) # Make contrast
              
              lrt <- edgeR::glmLRT(fit, contrast=pair_contrast)

              results <- as.data.frame(edgeR::topTags(lrt, p.value = 1, n = Inf, sort.by = "logFC"))
              
              colnames(results)[which(colnames(results) == "genes")] <- "Gene"
              colnames(results)[which(colnames(results) == "PValue")] <- "pvalue"
              
              results <- results[order(results$pvalue),]
              
              results <- results %>% mutate(Reference = rep(reference, n = nrow(results)),
                                            Comparison = rep(condition, n = nrow(results)))
              
              outfile <- paste(facet, "_", condition, "_vs_", reference,".txt", sep = "")
              
              write.table(results, paste(outdir, outfile, 
                                         sep = ""), sep = "\t")
              rm(results)
              
            } 
          
}





###########################################################################################################################
###########################################################################################################################
###########################################################################################################################



#Compute bulk DE by DESeq2:
DE.var = "Lesion"
facet_by.var = "Semi_Detailed_Celltype"

outdir = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/12_23_21_Analysis/Data/DE/Semi_Detailed_Celltype_Bulk_DE/"

Idents(DE.input) <- DE.var


  #Calc_agg_bulk:
  calc_agg_bulk_srt <- function (input, aggregate_by, group_by = FALSE, cutoff_frac = FALSE, 
                               cutoff_num = FALSE, cutoff_cpm = FALSE, expand_with_dummies = FALSE) {
  
  if (class(input) == "ExpressionSet"){
    
    if (group_by != FALSE) {
      ind <- match(group_by, aggregate_by)
      if (ind != 1) {
        stop("Please provide the group_by value first in the aggreggate_by argument")
      }
    }
    
    #Hide previous aggbulk columns from this analysis:
    check <- grep("bulk", colnames(fData(input)))
    if (length(check) > 0) {
      fData(input) <- fData(input)[, -check]
    }
    
    #Temporarily Convert pData Factors to Characters:
    if (any(sapply(pData(input), FUN = class) == "factor")) {
      
      factor.columns <- which(sapply(pData(input), FUN = class) == "factor")
      
      #Store factor columns for later:
      meta.data.factors <- pData(input)[,factor.columns]
      
      pData(input)[,factor.columns] <- sapply(pData(input)[,factor.columns], FUN = as.character)
    }
    
    #Setup to compute bulks:
    to_expand <- vector("list", length(aggregate_by))
    for (i in 1:length(aggregate_by)) {
      var <- aggregate_by[i]
      vars <- unique(pData(input)[, var])
      vars <- sort(vars)
      to_expand[[i]] <- vars
    }
    names(to_expand) <- aggregate_by
    bulks <- expand.grid(to_expand, stringsAsFactors = FALSE)
    colnames(bulks) <- c(aggregate_by)
    groups <- sort(unique(pData(input)[, group_by]))
    upm_vals <- c()
    num_cells_vals <- c()
    num_genes_vals <- c()
    rem_genes <- vector(mode = "list", length = nrow(bulks))
    
    #Compute Bulks:
    for (j in 1:nrow(bulks)) {
      int <- bulks[j, ]
      full_match <- c()
      for (k in 1:length(int)) {
        ind <- which(pData(input)[, colnames(bulks)[k]] == 
                       int[[k]])
        if (k == 1) {
          full_match <- c(full_match, ind)
        }
        else {
          full_match <- intersect(full_match, ind)
        }
      }
      num_cells <- length(full_match)
      num_cells_vals <- c(num_cells_vals, num_cells)
      #mike edit: changed length(full_match) to > "0" instead of "1":
      if (length(full_match) > 0) {
        #mike edit: added as.matrix to tmp line:
        tmp <- as.matrix(exprs(input)[, full_match])
        #mike edit: added if all(tmp == 0) statements for dummy cell treatments:
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == TRUE) {
          upm <- tmp
        }
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == FALSE) {
          upm <- (apply(tmp, 1, sum)/sum(tmp)) * 1e+06
        }
        if (cutoff_frac != FALSE || cutoff_num != FALSE) {
          zero_out_frac <- c()
          zero_out_num <- c()
          tmp2 <- tmp
          tmp2[which(tmp2 > 0)] <- 1
          gSums <- apply(tmp2, 1, sum)
          if (cutoff_frac != FALSE) {
            frac <- gSums/num_cells
            zero_out_frac <- which(frac < cutoff_frac)
          }
          if (cutoff_num != FALSE) {
            zero_out_num <- which(gSums < cutoff_num)
          }
          if (group_by != FALSE) {
            rem_genes[[j]] <- unique(names(c(zero_out_num, 
                                             zero_out_frac)))
          }
          else {
            upm[unique(c(zero_out_num, zero_out_frac))] <- 0
          }
        }
      }
      #mike edit: add dummy values for bulk permutations that don't exist in data:
      if (length(full_match) == 0) {
        tmp <- rep(0, times = nrow(exprs(input)))
        names(tmp) <- rownames(exprs(input))
        upm <- tmp
      }
      expressed <- length(which(upm > 0))
      upm_vals <- c(upm_vals, upm)
      num_genes_vals <- c(num_genes_vals, expressed)
    }
    
    names(rem_genes) <- apply(bulks, 1, FUN = paste, collapse = "__")
    bulks$numcells <- num_cells_vals
    bulks$numgenes <- num_genes_vals
    bulks$proportion <- 0
    if (group_by != FALSE) {
      for (i in 1:length(groups)) {
        ind <- grep(groups[i], bulks[, group_by])
        total <- sum(bulks$numcells[ind])
        bulks$proportion[ind] <- round((bulks$numcells[ind]/total) * 
                                         100, 2)
      }
    } else {
      total <- sum(bulks$numcells)
      bulks$proportion <- round((bulks$numcells/total) * 100, 
                                2)
    }
    bulk <- matrix(upm_vals, nrow = nrow(exprs(input)))
    rownames(bulk) <- rownames(exprs(input))
    colnames(bulk) <- seq(1:ncol(bulk))
    for (l in 1:nrow(bulks)) {
      cname <- bulks[l, -c(match(c("numcells", "proportion", 
                                   "numgenes"), colnames(bulks)))]
      cname2 <- c()
      for (i in 1:length(cname)) {
        cint <- as.character(cname[[i]])
        cname2 <- c(cname2, cint)
      }
      cname <- cname2
      cnum <- bulks[l, "numcells"]
      cpro <- bulks[l, "proportion"]
      cgen <- bulks[l, "numgenes"]
      cname <- paste0(c(cname, "num__genes", cgen, "num__cells", 
                        cnum, "percent", cpro, "bulk"), collapse = "__")
      colnames(bulk)[l] <- cname
    }
    if (group_by != FALSE) {
      if (!is.null(unlist(rem_genes))) {
        for (i in 1:length(vars)) {
          int_cell <- vars[i]
          ind <- grep(int_cell, names(rem_genes))
          vals <- table(unlist(rem_genes[ind]))
          zero_out <- names(which(vals == max(vals)))
          bulk[zero_out, ind] <- 0
        }
      }
    }
    if (cutoff_cpm != FALSE) {
      for (i in 1:length(vars)) {
        int_cell <- vars[i]
        ind <- grep(int_cell, colnames(bulk))
        if (length(ind) > 1) {
          gCount <- apply(bulk[, ind], 1, function(x) length(which(x >= cutoff_cpm)))
          zero_out <- which(gCount == 0)
        }
        else {
          zero_out <- which(bulk[, ind] < cutoff_cpm)
        }
        bulk[zero_out, ind] <- 0
      }
    }
    
    #Return pData Factor columns:
    if (("meta.data.factors" %in% ls()) == TRUE) {
      pData(input)[,factor.columns] <- meta.data.factors
    }
    
    #Remove non-existent dummy combinations:
    if (expand_with_dummies == FALSE) {
      
      no_cell_combos <- which(str_detect(colnames(bulk), pattern = "__num__cells__0__"))
      
      if (length(no_cell_combos) > 0){
        bulk <- bulk[,-no_cell_combos]  
      }
    }
    
    
    #Bind new agg bulk data:
    fData(input) <- cbind(fData(input), bulk)
  }
  
  if (class(input) == "Seurat") {
    
    if (group_by != FALSE) {
      ind <- match(group_by, aggregate_by)
      if (ind != 1) {
        stop("Please provide the group_by value first in the aggreggate_by argument")
      }
    }
    
    #Hide previous aggbulk columns from this analysis:
    check <- grep("bulk", colnames(input@assays$RNA@meta.features))
    if (length(check) > 0) {
      input@assays$RNA@meta.features <- input@assays$RNA@meta.features[, -check]
    }
    
    #Temporarily Convert Metadata Factors to Characters:
    if (any(sapply(input@meta.data, FUN = class) == "factor")) {
      
      factor.columns <- which(sapply(input@meta.data, FUN = class) == "factor")
      
      #Store factor columns for later:
      meta.data.factors <- input@meta.data[,factor.columns]
      
      input@meta.data[,factor.columns] <- sapply(input@meta.data[,factor.columns], FUN = as.character)
    }
    
    #Setup to compute bulks:
    to_expand <- vector("list", length(aggregate_by))
    for (i in 1:length(aggregate_by)) {
      var <- aggregate_by[i]
      vars <- unique(input@meta.data[, var])
      vars <- sort(vars)
      to_expand[[i]] <- vars
    }
    names(to_expand) <- aggregate_by
    bulks <- expand.grid(to_expand, stringsAsFactors = FALSE)
    colnames(bulks) <- c(aggregate_by)
    groups <- sort(unique(input@meta.data[, group_by]))
    upm_vals <- c()
    num_cells_vals <- c()
    num_genes_vals <- c()
    rem_genes <- vector(mode = "list", length = nrow(bulks))
    
    #Compute Bulks:
    for (j in 1:nrow(bulks)) {
      int <- bulks[j, ]
      full_match <- c()
      for (k in 1:length(int)) {
        ind <- which(input@meta.data[, colnames(bulks)[k]] == 
                       int[[k]])
        if (k == 1) {
          full_match <- c(full_match, ind)
        }
        else {
          full_match <- intersect(full_match, ind)
        }
      }
      num_cells <- length(full_match)
      num_cells_vals <- c(num_cells_vals, num_cells)
      
      #mike edit: changed length(full_match) to > "0" instead of "1":
      if (length(full_match) > 0) {
        #mike edit: added as.matrix to tmp line:
        tmp <- as.matrix(input@assays$RNA@data[, full_match])
        #mike edit: added if all(tmp == 0) statements for dummy cell treatments:
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == TRUE) {
          upm <- tmp
        }
        if (any(is.nan(apply(tmp, 1, sum)/sum(tmp))) == FALSE) {
          upm <- (apply(tmp, 1, sum)/sum(tmp)) * 1e+06
        }
        if (cutoff_frac != FALSE || cutoff_num != FALSE) {
          zero_out_frac <- c()
          zero_out_num <- c()
          tmp2 <- tmp
          tmp2[which(tmp2 > 0)] <- 1
          gSums <- apply(tmp2, 1, sum)
          if (cutoff_frac != FALSE) {
            frac <- gSums/num_cells
            zero_out_frac <- which(frac < cutoff_frac)
          }
          if (cutoff_num != FALSE) {
            zero_out_num <- which(gSums < cutoff_num)
          }
          if (group_by != FALSE) {
            rem_genes[[j]] <- unique(names(c(zero_out_num, 
                                             zero_out_frac)))
          }
          else {
            upm[unique(c(zero_out_num, zero_out_frac))] <- 0
          }
        }
      }
      #mike edit: add dummy values for bulk permutations that don't exist in data:
      if (length(full_match) == 0) {
        tmp <- rep(0, times = nrow(input@assays$RNA@data))
        names(tmp) <- rownames(input@assays$RNA@data)
        upm <- tmp
      }
      expressed <- length(which(upm > 0))
      upm_vals <- c(upm_vals, upm)
      num_genes_vals <- c(num_genes_vals, expressed)
    }
    
    names(rem_genes) <- apply(bulks, 1, FUN = paste, collapse = "__")
    bulks$numcells <- num_cells_vals
    bulks$numgenes <- num_genes_vals
    bulks$proportion <- 0
    if (group_by != FALSE) {
      for (i in 1:length(groups)) {
        ind <- grep(groups[i], bulks[, group_by])
        total <- sum(bulks$numcells[ind])
        bulks$proportion[ind] <- round((bulks$numcells[ind]/total) * 
                                         100, 2)
      }
    } else {
      total <- sum(bulks$numcells)
      bulks$proportion <- round((bulks$numcells/total) * 100, 
                                2)
    }
    bulk <- matrix(upm_vals, nrow = nrow(input@assays$RNA@data))
    rownames(bulk) <- rownames(input@assays$RNA@data)
    colnames(bulk) <- seq(1:ncol(bulk))
    for (l in 1:nrow(bulks)) {
      cname <- bulks[l, -c(match(c("numcells", "proportion", 
                                   "numgenes"), colnames(bulks)))]
      cname2 <- c()
      for (i in 1:length(cname)) {
        cint <- as.character(cname[[i]])
        cname2 <- c(cname2, cint)
      }
      cname <- cname2
      cnum <- bulks[l, "numcells"]
      cpro <- bulks[l, "proportion"]
      cgen <- bulks[l, "numgenes"]
      cname <- paste0(c(cname, "num__genes", cgen, "num__cells", 
                        cnum, "percent", cpro, "bulk"), collapse = "__")
      colnames(bulk)[l] <- cname
    }
    if (group_by != FALSE) {
      if (!is.null(unlist(rem_genes))) {
        for (i in 1:length(vars)) {
          int_cell <- vars[i]
          ind <- grep(int_cell, names(rem_genes))
          vals <- table(unlist(rem_genes[ind]))
          zero_out <- names(which(vals == max(vals)))
          bulk[zero_out, ind] <- 0
        }
      }
    }
    if (cutoff_cpm != FALSE) {
      for (i in 1:length(vars)) {
        int_cell <- vars[i]
        ind <- grep(int_cell, colnames(bulk))
        if (length(ind) > 1) {
          gCount <- apply(bulk[, ind], 1, function(x) length(which(x >= cutoff_cpm)))
          zero_out <- which(gCount == 0)
        }
        else {
          zero_out <- which(bulk[, ind] < cutoff_cpm)
        }
        bulk[zero_out, ind] <- 0
      }
    }
    
    #Return Metadata Factor Columns:
    if (("meta.data.factors" %in% ls()) == TRUE) {
      input@meta.data[,factor.columns] <- meta.data.factors
    }
    
    #Remove non-existent dummy combinations:
    if (expand_with_dummies == FALSE) {
      
      no_cell_combos <- which(str_detect(colnames(bulk), pattern = "__num__cells__0__"))
      
      if (length(no_cell_combos) > 0){
        bulk <- bulk[,-no_cell_combos]  
      }
      
    }
    
    #Bind new agg bulk data:
    input@assays$RNA@meta.features<- cbind(input@assays$RNA@meta.features, bulk)
    
    
  }
  
  return(input)
  
}

  
  DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c(facet_by.var, DE.var, "Patient"))
  
  Bulk_data <- DE.input@assays$RNA@meta.features[,which(str_detect(colnames(DE.input@assays$RNA@meta.features), pattern = "bulk"))]
  
  facets <- as.character(unique(DE.input@meta.data[,facet_by.var]))
  
  reference = unique(comparisons$reference)
  
  
  for (c in 1:length(facets)) {
    
    facet <- facets[c]
    
          # Subset by facet:
          facet_string <- as.character(sapply(colnames(Bulk_data), FUN = function(x) str_split(x, pattern = "__")[[1]][which(str_split(x, pattern = "__")[[1]] %in% facets)]))
          facet_subset <- Bulk_data[,which(facet_string == facet)] 
          
          #Create Metadata Dataframe:
          DE.vars <- as.character(unique(DE.input@meta.data[,DE.var]))
          DE.var_string <- as.character(sapply(colnames(Bulk_data), FUN = function(x) str_split(x, pattern = "__")[[1]][which(str_split(x, pattern = "__")[[1]] %in% DE.vars)]))
          
          Patients <- as.character(unique(DE.input@meta.data$Patient))
          Patient_string <- as.character(sapply(colnames(Bulk_data), FUN = function(x) str_split(x, pattern = "__")[[1]][which(str_split(x, pattern = "__")[[1]] %in% Patients)]))
          
          metadata <- data.frame(facet = facet_string[which(facet_string == facet)],
                                 DE.var = DE.var_string[which(facet_string == facet)],
                                 Patient = Patient_string[which(facet_string == facet)])
          colnames(metadata) <- c(facet_by.var, DE.var, "Patient")
          rownames(metadata) <- colnames(facet_subset)
          
          if (levels(DE.input@meta.data[,DE.var])[1] != reference) {
            stop("set levels of DE.var so that reference condition is position 1")
          }
          
          metadata[,DE.var] <- factor(metadata[,DE.var], 
                                      levels = levels(DE.input@meta.data[,DE.var]))
          
          #Run DESeq2:
          dds <-DESeqDataSetFromMatrix(countData = as.matrix(round(facet_subset)),
                                       colData = metadata,
                                       design = formula(paste("~", DE.var))) # Note: add Patient to design formula to include patient as confounding variable
          
          dds <- DESeq(dds)
          
          DE.conditions <- DE.vars[which(DE.vars %in% comparisons$comparison)]
          
              for (condition in DE.conditions) {
                  results <- results(dds,
                                     contrast = c(DE.var, condition, reference), 
                                     alpha = 0.05)
                  results <- as.data.frame(results[order(results$pvalue),])
                  results <- rownames_to_column(results, var = "Gene")
                  results <- results %>% filter(!is.na(padj))
                  
                  results <- results %>% mutate(Reference = rep(reference, n = nrow(results)),
                                                Comparison = rep(condition, n = nrow(results)))
                  
                  outfile <- paste(facet, "_", condition, "_vs_", reference,".txt", sep = "")
                  
                  write.table(results, paste(outdir, outfile, 
                                                 sep = ""), sep = "\t")
                rm(results)
              }
          
        }
  




