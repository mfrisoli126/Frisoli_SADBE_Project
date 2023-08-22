#######################################################
#Prep Data for Network Calc by CellphoneDB:
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
library(clusterProfiler)
library(magrittr)
library(aliases2entrez)
library(data.table)
library(Matrix)
library(ggrepel)
library(purrr)
library(readxl)
library(RColorBrewer)

setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

# Load Data:
load(file = "Data/contact_derm_filtered3.Rdata")

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
        DE.input <- calc_agg_bulk_srt(DE.input, aggregate_by = c("Semi_Detailed_Celltype", "Lesion"))
      
        aggregate.counts <- DE.input@assays$RNA@meta.features[,which(str_detect(string = colnames(DE.input@assays$RNA@meta.features), pattern = "bulk"))]

            #Gene Fix:
            gene.fix.df <- as.data.frame(rownames(aggregate.counts))
            colnames(gene.fix.df) = "Original_symbol"
            
            file <- system.file("extdata", "HGNC.txt", package = "aliases2entrez")
            HGNC <- fread("gene_lists/HGNC.tsv")
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
            
            gene.fix.df <- full_join(gene.fix.df, entrez.to.ensembl)
    
    
             #Convert to Counts by ENSEMBL ID:  
                
                #Filter to genes with positive ENSEMBL ID:
                positive_ids <- gene.fix.df[-which(is.na(gene.fix.df$ENSEMBL)),]
                
                #Many genes have several ENTREZ IDs and ENSEMBL IDs. Filter to only first line for each gene symbol:
                positive_ids <- positive_ids[-which(duplicated(positive_ids$Original_symbol)),]
                
                aggregate.counts <- aggregate.counts[which(rownames(aggregate.counts) %in% positive_ids$Original_symbol),]
                
                #For some reason, different gene symbols have same ENSEMBL ID:
                dup.ensembleIDs <- positive_ids[which(duplicated(positive_ids$ENSEMBL)),]
                dup.ensembleIDs <- positive_ids[which(positive_ids$ENSEMBL %in% dup.ensembleIDs$ENSEMBL),]
                
                counts.by.ensembleID <- aggregate.counts[-which(rownames(aggregate.counts) %in% dup.ensembleIDs$Original_symbol),]
                new.rownames <- positive_ids[which(positive_ids$Original_symbol %in% rownames(counts.by.ensembleID)),]
                rownames(counts.by.ensembleID) <- new.rownames$ENSEMBL[match(new.rownames$Original_symbol, rownames(counts.by.ensembleID))]
                
                      for (ID in 1:length(unique(dup.ensembleIDs$ENSEMBL))) {
                        
                        id <- unique(dup.ensembleIDs$ENSEMBL)[ID]
                        
                        symbols <- positive_ids %>% filter(ENSEMBL == id) %>% dplyr::select(Original_symbol) %>% unlist()
                        ensembleID <- positive_ids %>% filter(ENSEMBL == id) %>% dplyr::select(ENSEMBL) %>% unlist() %>% unique()
                        
                        data <- aggregate.counts[which(rownames(aggregate.counts) %in% symbols),]
                        merged.data <- bind_rows(data, apply(data, MARGIN = 2, FUN = sum))
                        rownames(merged.data)[nrow(merged.data)] <- ensembleID
                        
                        counts.by.ensembleID<- bind_rows(counts.by.ensembleID, merged.data[nrow(merged.data),])
                        
          }
                      
                      counts.by.ensembleID <- rownames_to_column(counts.by.ensembleID, var = "Gene")
                      
                      colnames <- str_split(colnames(counts.by.ensembleID)[-1], pattern = "__", n =3, simplify = TRUE)
                      colnames <- as.data.frame(colnames, ncol = 3)[,-3]
                      colnames <- apply(colnames, MARGIN = 1, FUN = function(x) str_c(string = x, collapse = "__"))
                      
                      colnames(counts.by.ensembleID)[-1] <- colnames
    
          write.table(counts.by.ensembleID, file = "Data/RL_network2/aggregate.grouped.counts.by.ensembleID.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    
          #Create Metadata txt file:
                data_IDs <- colnames(counts.by.ensembleID)[-1] 
                
                les.string <- c()
                cell.string <- c()
                for (i in 1:length(data_IDs)) {
                  
                  id <- data_IDs[i]
                  
                  string <- str_split(id, pattern = "__",) %>% unlist()
                  les <- string[1]
                  cell <- string[2]
                  
                  les.string <- c(les.string, les)
                  cell.string <- c(cell.string, cell)
                }
                
                meta.data <- data.frame(Cell = data_IDs,
                                        cell_type = cell.string,
                                        Lesion = les.string)
                
                meta.data <- meta.data %>% mutate(cell_type = paste(cell_type, "__", Lesion, sep = ""))
                meta.data <- meta.data[,-3]
          
          write.table(meta.data, file = "Data/RL_network2/aggregate.metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)
          
  # Individual Cell Data:
              counts <- as.data.frame(contact_derm@assays$RNA@counts)
              
              #Gene Fix:
              gene.fix.df <- as.data.frame(rownames(counts))
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
              
              gene.fix.df <- full_join(gene.fix.df, entrez.to.ensembl)
              
              
              #Convert to Counts  by ENSEMBL ID:  
                    
                    #Filter to genes with positive ENSEMBL ID:
                    positive_ids <- gene.fix.df[-which(is.na(gene.fix.df$ENSEMBL)),]
                    
                    #Many genes have several ENTREZ IDs and ENSEMBL IDs. Filter to only first line for each gene symbol:
                    positive_ids <- positive_ids[-which(duplicated(positive_ids$Original_symbol)),]
                    
                    counts <- counts[which(rownames(counts) %in% positive_ids$Original_symbol),]
                    
                    #For some reason, different gene symbols have same ENSEMBL ID:
                    dup.ensembleIDs <- positive_ids[which(duplicated(positive_ids$ENSEMBL)),]
                    dup.ensembleIDs <- positive_ids[which(positive_ids$ENSEMBL %in% dup.ensembleIDs$ENSEMBL),]
                    
                    counts.by.ensembleID <- counts[-which(rownames(counts) %in% dup.ensembleIDs$Original_symbol),]
                    new.rownames <- positive_ids[which(positive_ids$Original_symbol %in% rownames(counts.by.ensembleID)),]
                    rownames(counts.by.ensembleID) <- new.rownames$ENSEMBL[match(new.rownames$Original_symbol, rownames(counts.by.ensembleID))]
                    
                    for (ID in 1:length(unique(dup.ensembleIDs$ENSEMBL))) {
                      
                      id <- unique(dup.ensembleIDs$ENSEMBL)[ID]
                      
                      symbols <- positive_ids %>% filter(ENSEMBL == id) %>% dplyr::select(Original_symbol) %>% unlist()
                      ensembleID <- positive_ids %>% filter(ENSEMBL == id) %>% dplyr::select(ENSEMBL) %>% unlist() %>% unique()
                      
                      data <- counts[which(rownames(counts) %in% symbols),]
                      merged.data <- bind_rows(data, apply(data, MARGIN = 2, FUN = sum))
                      rownames(merged.data)[nrow(merged.data)] <- ensembleID
                      
                      counts.by.ensembleID<- bind_rows(counts.by.ensembleID, merged.data[nrow(merged.data),])
                      
                    }
                    
                    counts.by.ensembleID <- rownames_to_column(counts.by.ensembleID, var = "Gene")
                    
              write.table(counts.by.ensembleID, file = "Data/RL_network2/counts_by_ensembleID.txt", sep = "\t", row.names = FALSE, quote = FALSE)
                          
                    #Create Metadata txt file:
                      metadata <- contact_derm@meta.data
                      
                      metadata <- metadata[,which(colnames(metadata) %in% c("Semi_Detailed_Celltype", "Lesion"))]
                      metadata <- rownames_to_column(metadata, var = "Cell")
                      colnames(metadata)[3] <- "cell_type"
                      
                      metadata <- metadata %>% mutate(cell_type = paste(cell_type, "__", Lesion, sep = ""))
                      metadata <- metadata[,-2]
              
              write.table(metadata, file = "Data/RL_network2/metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)
              
    














