###########################################################################################################
### Import without pre-filtering

###########################################################################################################

setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/")

library(dplyr)
library(data.table)
library(purrr)
library(tibble)
library(stringr)
library(ggplot2)
library(readxl)


#Sequencing Sample record:
sequencing_record <- read_xlsx(path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/Filtered_raw_data/sequencing_sample_record.xlsx")
sequencing_record$Sample <- str_replace_all(sequencing_record$Sample, pattern = "[:space:]", replacement = "_")

sample_list <- c()

#Contact Derm Patient Data:
    path = "/nl/umw_manuel_garber/yuqing/inDrop/Skin/dolphin2/Mikedata/"
    
    reports = list.files(path = path, pattern = "report")
    
        #Create tmp files for each sample:
        
        for (r in 1:length(reports)) {
          
          #Get File list for each report:
          report <- reports[r]

          UMIfolder <- "/UMI_count_final_after_star/"
          
          filelist <- list.files(path = paste(path, report, UMIfolder, sep = ""), pattern = "umiClean")
          
          if (length(filelist) == 0) {
            UMIfolder <- "/UMI_count_after_star/"
            filelist <- list.files(path = paste(path, report, UMIfolder, sep = ""), pattern = "umiClean")
          }
          
                #Import file data:
                for (f in 1:length(filelist)) {
                    sample_input <- filelist[f]
                    
                    sample <- strsplit(sample_input, split = "_")[[1]]
                    sample <- str_c(sample[-which(str_detect(sample, pattern = "umiClean"))], collapse = "_")
                    
                    sequence_date <- sequencing_record %>% filter(Report_num == report, Sample == sample) %>% select(Sequence_Date) %>% unlist() %>% unique()
                    
                    sample <- paste(sample, "_", sequence_date, sep = "") 
                    
                    print(paste("Starting", "Sample ", f, ":", sample))
                    
                    sample.tmp <- fread(file = paste(path, report, UMIfolder, sample_input, sep = ""), sep = "\t", header = T)
                    sample.tmp <- as.data.frame(sample.tmp)
                    rownames(sample.tmp) <- sample.tmp[,1]
                    sample.tmp <- as.data.frame(t(sample.tmp[,2:ncol(sample.tmp)]))
                    
                    #Add sample label 
                    sample.tmp$Sample <- replicate(n = nrow(sample.tmp), expr = paste(sample))
                    
                    #Add cell summary stats:
                    sample.tmp$UMI.Sum <- apply(sample.tmp[, -which(names(sample.tmp) %in% c("Sample"))], MARGIN = 1, FUN = sum)
                    sample.tmp$gene.Sum <- apply(sample.tmp[, -which(names(sample.tmp) %in% c("Sample", "UMI.Sum"))], MARGIN = 1, FUN = function(x) length(which(x != 0)))
                    
                    mito.genes <- grep(pattern = "^MT-", x = colnames(x = sample.tmp), value = TRUE)
                    sample.tmp$mito.UMI.Sum <- apply(sample.tmp[, which(names(sample.tmp) %in% mito.genes)],1,sum)
                    sample.tmp$percent.mt <- (sample.tmp$mito.UMI.Sum/sample.tmp$UMI.Sum)*100
                    
                    #Create new Sample_Cell_ID Label:
                    sample.tmp <- rownames_to_column(sample.tmp)
                    sample.tmp <- sample.tmp %>% rename("Cell_ID" = rowname)
                    sample.tmp$Cell_ID <- paste(sample.tmp$Sample, sample.tmp$Cell_ID, sep = "_")
                    
                    #Keep only summary stats:
                    sample.tmp <- sample.tmp[,which(colnames(sample.tmp) %in% c("Cell_ID", "UMI.Sum", "gene.Sum", "mito.UMI.Sum", "percent.mt"))]
                    
                    #Assign Data Label:
                    assign(sample, sample.tmp)
                    
                    sample_list <- c(sample_list, sample)
                    
                    rm(sample.tmp)
                }
          
        }


#Healthy Control Data:

      path = "/nl/umw_manuel_garber/human/skin/reprocessing/process/"
      
      yuqing_sample_folders = list.files(path = path)
      
      #exclude folders that don't seem to be related to samples:
      yuqing_sample_folders <- yuqing_sample_folders[which(str_detect(string = yuqing_sample_folders, pattern = "_"))]
      
      #Only include healthy controls:
      yuqing_sample_folders <- yuqing_sample_folders[which(str_detect(string = yuqing_sample_folders, pattern = "^CB"))]
            
            for (y in 1:length(yuqing_sample_folders)) {
              
              report_folder <- list.files(path = paste(path, yuqing_sample_folders[y], sep = ""), pattern = "report")
              report_folder <- report_folder[1]
              
              report_sub_folders <- list.files(path = paste(path, yuqing_sample_folders[y], "/", report_folder, sep = ""))
              
              if (any(str_detect(string = report_sub_folders, pattern = "UMI_count_final_after_star")) == TRUE) {
                
                filelist <- list.files(path = paste(path, yuqing_sample_folders[y], "/", report_folder, "/UMI_count_final_after_star", sep = ""))
                
                      #Import file data:
                      for (f in 1:length(filelist)) {
                        sample_input <- filelist[f]
                        sample <- strsplit(sample_input, split = "_")[[1]]
                        sample <- paste0(sample[1], "_", sample[2], "_", sample[3],"_", sample[7])
                        
                        print(paste("Starting", "Sample ", y, ":", sample))
                        
                        sample.tmp <- fread(file = paste(path, yuqing_sample_folders[y], "/", report_folder, "/UMI_count_final_after_star/", sample_input, sep = ""), sep = "\t", header = T)
                        
                        sample.tmp <- as.data.frame(sample.tmp)
                        rownames(sample.tmp) <- sample.tmp[,1]
                        sample.tmp <- as.data.frame(t(sample.tmp[,2:ncol(sample.tmp)]))
                        
                        #Add sample label 
                        sample.tmp$Sample <- replicate(n = nrow(sample.tmp), expr = paste(sample))
                        
                        #Add cell summary stats:
                        sample.tmp$UMI.Sum <- apply(sample.tmp[, -which(names(sample.tmp) %in% c("Sample"))], MARGIN = 1, FUN = sum)
                        sample.tmp$gene.Sum <- apply(sample.tmp[, -which(names(sample.tmp) %in% c("Sample", "UMI.Sum"))], MARGIN = 1, FUN = function(x) length(which(x != 0)))
                        
                        mito.genes <- grep(pattern = "^MT-", x = colnames(x = sample.tmp), value = TRUE)
                        sample.tmp$mito.UMI.Sum <- apply(sample.tmp[, which(names(sample.tmp) %in% mito.genes)],1,sum)
                        sample.tmp$percent.mt <- (sample.tmp$mito.UMI.Sum/sample.tmp$UMI.Sum)*100
                        
                        #Create new Sample_Cell_ID Label:
                        sample.tmp <- rownames_to_column(sample.tmp)
                        sample.tmp <- sample.tmp %>% rename("Cell_ID" = rowname)
                        sample.tmp$Cell_ID <- paste(sample.tmp$Sample, sample.tmp$Cell_ID, sep = "_")
                        
                        
                        #Keep only summary stats:
                        sample.tmp <- sample.tmp[,which(colnames(sample.tmp) %in% c("Cell_ID", "UMI.Sum", "gene.Sum", "mito.UMI.Sum", "percent.mt"))]
                        
                        
                        #Assign Data Label:
                        assign(sample, sample.tmp)
                        
                        sample_list <- c(sample_list, sample)
                        
                        rm(sample.tmp)
                      }
            }
      }

## Merge All data:
      
  data.summary.stats <- reduce(mget(sample_list), full_join)
  print(paste("Merge Complete, Summary Stats Data Size is...", dim(data.summary.stats)[1], "rows   x ", dim(data.summary.stats)[2], "columns"))
  
  data.summary.stats$Cell_ID <- str_replace_all(data.summary.stats$Cell_ID, pattern = "__", replacement = "_")
  
  #Add Metadata:
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
  
    meta.data <- as.data.frame(matrix(unlist(t(sapply(data.summary.stats$Cell_ID, FUN = split_cellID_to_elements))), ncol = 3))
    colnames(meta.data) <- c("Patient", "Lesion", "Sequence_date")
    
    data.summary.stats <- cbind(data.summary.stats, meta.data)
    
  
  #Save Folder Location:
  save.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/"
  
  write.table(data.summary.stats, file = paste(save.path, "data.summary.stats.tsv", sep = ""), row.names=FALSE, sep="\t")
  
#######################
# Plot Results
#######################
library(dplyr)
library(data.table)
library(purrr)
library(tibble)
library(stringr)
library(ggplot2)
library(readxl)

setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/")

#Sequencing Sample record:
sequencing_record <- read_xlsx(path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/Filtered_raw_data/sequencing_sample_record.xlsx")
sequencing_record$Sample <- str_replace_all(sequencing_record$Sample, pattern = "[:space:]", replacement = "_")

file.dims <- read.delim("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/file_dims_5_6_22.tsv")

data.summary.stats <- read.delim("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/data.summary.stats.tsv")

#Plot overall file stats:

merged.file.dims <- file.dims[-which(str_detect(file.dims$Sample, pattern = "_A_|_B_")),]
split.samples <- file.dims[which(str_detect(file.dims$Sample, pattern = "_A_|_B_")),]
  
  for (r in 1:nrow(split.samples)) {
    
    row <- split.samples[r,]
    
    if (str_detect(row$Sample, pattern = "_B_") == TRUE) {
      
      sample.name <- str_replace(row$Sample, pattern = "_B_", replacement = "_")
      
      name.list.length <- as.numeric(lengths(str_split(sample.name, pattern = "_")))
      sequence_date <- str_c(unlist(str_split(sample.name, pattern = "_"))[(name.list.length-2):(name.list.length)], collapse = "_")
      sample.patient.lesion <- str_c(unlist(str_split(sample.name, pattern = "_"))[1:(name.list.length-3)], collapse = "_")

      patient.lesion.rows <- split.samples[which(str_detect(split.samples$Sample, pattern = sample.patient.lesion)),]      
      
            if (nrow(patient.lesion.rows) == 2) {
              
              merged.data <- data.frame(Sample = str_c(c(sample.patient.lesion, sequence_date), collapse = "_"),
                                        rows = sum(patient.lesion.rows$rows),
                                        columns = max(patient.lesion.rows$columns))
              
              merged.file.dims <- bind_rows(merged.file.dims, merged.data)
              
            }
            
            if (nrow(patient.lesion.rows) > 2) {
              
              patient.lesion.date.rows <- patient.lesion.rows[which(str_detect(patient.lesion.rows$Sample, pattern = sequence_date)),]  
              
              merged.data <- data.frame(Sample = str_c(c(sample.patient.lesion, sequence_date), collapse = "_"),
                                        rows = sum(patient.lesion.date.rows$rows),
                                        columns = max(patient.lesion.date.rows$columns))
              
              merged.file.dims <- bind_rows(merged.file.dims, merged.data)
              
            }
      
      
    }
  }



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
  
  meta.data <- as.data.frame(matrix(unlist(t(sapply(merged.file.dims$Sample, FUN = split_sample_to_elements))), ncol = 3))
  colnames(meta.data) <- c("Patient", "Lesion", "Sequence_date")
  
  merged.file.dims <- cbind(merged.file.dims, meta.data)
  
  colnames(merged.file.dims)[which(colnames(merged.file.dims) %in% c("rows", "columns"))] <- c("passable_cells", "genes")
  
  merged.file.dims$Sequence_date <- factor(merged.file.dims$Sequence_date, 
                                    levels = c("4_16_20", "1_20_21", "12_3_21", "4_28_22"))

  #Reformat:
  merged.file.dims$Lesion[which(merged.file.dims$Lesion == "D2_SADBE")] <- "Day2_SADBE"
  merged.file.dims$Lesion[which(merged.file.dims$Lesion == "D4_SADBE")] <- "Day4_SADBE"
  merged.file.dims$Lesion[which(merged.file.dims$Lesion == "SLS")] <- "Day2_SLS"
  merged.file.dims$Lesion[which(merged.file.dims$Lesion == "Acetone")] <- "Acetone_Vehicle"
  

  merged.file.dims$Lesion <- factor(merged.file.dims$Lesion, 
                                           levels = c("NL", "Day2_SLS", "Acetone_Vehicle", "Day2_SADBE", "Day4_SADBE"))
  
  merged.file.dims <- merged.file.dims %>% mutate(homogenous_sample_name = paste(Patient, Lesion, Sequence_date, sep = "_"))
  
  #Add Cell Count Data:
  blister.cell.counts <- read_xlsx(path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/Filtered_raw_data/master_SADBE_sample_info.xlsx")

      #Reformat:
      blister.cell.counts$Blister_Date <- as.Date(blister.cell.counts$Blister_Date)
      blister.cell.counts$Blister_Date <- str_replace_all(format(as.Date(blister.cell.counts$Blister_Date), "%m_%d_%y"), pattern = "^0", replacement = "")
      blister.cell.counts$Blister_Date <- str_replace_all(blister.cell.counts$Blister_Date, pattern = "_0", replacement = "_")
      blister.cell.counts$Lesion <- str_replace(blister.cell.counts$Lesion, pattern = "Day ", replacement = "Day")
      blister.cell.counts$Lesion <- str_replace(blister.cell.counts$Lesion, pattern = "D2", replacement = "Day2")
      blister.cell.counts$Lesion <- str_replace(blister.cell.counts$Lesion, pattern = "D4", replacement = "Day4")
      blister.cell.counts$Lesion <- str_replace(blister.cell.counts$Lesion, pattern = "\\s", replacement = "_")
      blister.cell.counts$Lesion <- str_replace(blister.cell.counts$Lesion, pattern = "Day2_Acetone", replacement = "Acetone_Vehicle")
      
      
  blister.cell.counts <- blister.cell.counts %>% filter(Sequenced == "yes")
  blister.cell.counts <- blister.cell.counts %>% mutate(homogenous_sample_name = paste(Patient, Lesion, `Sequence Date`, sep = "_"))
  
  merged.file.dims <- full_join(merged.file.dims, 
                                blister.cell.counts[,which(colnames(blister.cell.counts) %in% c("Blister_Date", "sequence_file_name", "Sample Cell Concentration per mL", "Total cells", "homogenous_sample_name"))])
  
  colnames(merged.file.dims)[which(colnames(merged.file.dims) %in% c("Sample Cell Concentration per mL", "Total cells"))] <- c("Blister_cell_concentration","Blister_cell_count")
  
  for (r in 1:nrow(merged.file.dims)) {
    if (is.na(merged.file.dims$Sample[r])) {
      merged.file.dims$Sample[r] <- merged.file.dims$homogenous_sample_name[r]
      merged.file.dims$passable_cells[r] <- 0
      merged.file.dims$genes[r] <- 0
      
      name.list.length <- as.numeric(lengths(str_split(merged.file.dims$homogenous_sample_name[r], pattern = "_")))
      merged.file.dims$Sequence_date[r] <- str_c(unlist(str_split(merged.file.dims$homogenous_sample_name[r], pattern = "_"))[(name.list.length-2):(name.list.length)], collapse = "_")
      merged.file.dims$Patient[r] <- unlist(str_split(merged.file.dims$homogenous_sample_name[r], pattern = "_"))[1]
      merged.file.dims$Lesion[r] <- str_c(unlist(str_split(merged.file.dims$homogenous_sample_name[r], pattern = "_"))[2:(name.list.length-3)], collapse = "_")
    }
  }
  
  
  #Add Duplication Rates:
  
    overall_summary <- data.frame(matrix(nrow = 0, ncol = 12))
    colnames(overall_summary) <- c("sample", "Total.Reads", "Valid.Reads", "Multimapped.Reads.Aligned..STAR.", "Unique.Reads.Aligned..STAR.",
                                   "Total.aligned.UMIs..ESAT.", "Total.deduped.UMIs..ESAT.", "Duplication.Rate", "Number.of.Barcodes", "Mean.UMIs.per.Barcode",
                                   "Number.of.Genes", "Mean.Genes.per.Barcode")
    overall_summary[,1] <- as.character(overall_summary[,1])
    overall_summary[,c(2:7,9,11)] <- lapply(overall_summary[,c(2:7,9,11)], FUN = as.integer)
    overall_summary[,c(8,10,12)] <- lapply(overall_summary[,c(8,10,12)], FUN = as.numeric)
    
          #Contact Derm Patient Data:
          path = "/nl/umw_manuel_garber/yuqing/inDrop/Skin/dolphin2/Mikedata/"
          
          reports = list.files(path = path, pattern = "report")
          
              for (r in 1:length(reports)) {
                report <- reports[r]
                
                report_summary <- read.delim(file = paste(path, reports[r],"/summary/overall_summary.tsv", sep = ""))
              
                samplelist <- report_summary$sample
                
                reformated.samplelist <- c()
                      for (s in 1:length(samplelist)) {
                        sample <- samplelist[s]
                        
                        sequence_date <- sequencing_record %>% filter(Report_num == report, Sample == sample) %>% select(Sequence_Date) %>% unlist() %>% unique()
                        
                        sample <- paste(sample, "_", sequence_date, sep = "") 
                        
                        reformated.samplelist <- c(reformated.samplelist, sample)
                        rm(sample)
                      }
                
                report_summary$sample <- reformated.samplelist
                
                overall_summary <- bind_rows(overall_summary, report_summary)
                
                rm(report_summary)
              }
          
          #Healthy Control Data:
          
          path = "/nl/umw_manuel_garber/human/skin/reprocessing/process/"
          
          yuqing_sample_folders = list.files(path = path)
          
          #exclude folders that don't seem to be related to samples:
          yuqing_sample_folders <- yuqing_sample_folders[which(str_detect(string = yuqing_sample_folders, pattern = "_"))]
          
          #Only include healthy controls:
          yuqing_sample_folders <- yuqing_sample_folders[which(str_detect(string = yuqing_sample_folders, pattern = "^CB"))]
          
              for (y in 1:length(yuqing_sample_folders)) {
                
                report_folder <- list.files(path = paste(path, yuqing_sample_folders[y], sep = ""), pattern = "report")
                report_folder <- report_folder[1]
                
                report_sub_folders <- list.files(path = paste(path, yuqing_sample_folders[y], "/", report_folder, sep = ""))
                
                    if (any(str_detect(string = report_sub_folders, pattern = "UMI_count_final_after_star")) == TRUE) {
                      
                      report_summary <- read.delim(file = paste(path, yuqing_sample_folders[y], "/", report_folder, "/summary/overall_summary.tsv", sep = ""))
    
                      report_summary$sample <- str_replace(report_summary$sample, pattern = "_Bst_sc_iD_", replacement = "_")
                        
                      overall_summary <- bind_rows(overall_summary, report_summary)
                      
                      rm(report_summary)
                      
                    }
              }
  
  
        # Merge split emulsions for comparison analysis:
          merged.sequencing.summary <- overall_summary[-which(str_detect(overall_summary$sample, pattern = "_A_|_B_")),]
          split.samples <- overall_summary[which(str_detect(overall_summary$sample, pattern = "_A_|_B_")),]
          
          for (r in 1:nrow(split.samples)) {
            
            row <- split.samples[r,]
            
            if (str_detect(row$sample, pattern = "_B_") == TRUE) {
              
              sample.name <- str_replace(row$sample, pattern = "_B_", replacement = "_")
              
              name.list.length <- as.numeric(lengths(str_split(sample.name, pattern = "_")))
              sequence_date <- str_c(unlist(str_split(sample.name, pattern = "_"))[(name.list.length-2):(name.list.length)], collapse = "_")
              sample.patient.lesion <- str_c(unlist(str_split(sample.name, pattern = "_"))[1:(name.list.length-3)], collapse = "_")
              
              patient.lesion.rows <- split.samples[which(str_detect(split.samples$sample, pattern = sample.patient.lesion)),]      
              
              if (nrow(patient.lesion.rows) == 2) {
                
                merged.data <- data.frame(sample = str_c(c(sample.patient.lesion, sequence_date), collapse = "_"),
                                          Total.Reads = sum(patient.lesion.rows$Total.Reads),
                                          Valid.Reads = sum(patient.lesion.rows$Valid.Reads),
                                          Multimapped.Reads.Aligned..STAR. = sum(patient.lesion.rows$Multimapped.Reads.Aligned..STAR.),
                                          Unique.Reads.Aligned..STAR. = mean(patient.lesion.rows$Unique.Reads.Aligned..STAR.),
                                          Total.aligned.UMIs..ESAT. = mean(patient.lesion.rows$Total.aligned.UMIs..ESAT.),
                                          Total.deduped.UMIs..ESAT. = mean(patient.lesion.rows$Total.deduped.UMIs..ESAT.),
                                          Duplication.Rate = mean(patient.lesion.rows$Duplication.Rate),
                                          Number.of.Barcodes = mean(patient.lesion.rows$Number.of.Barcodes),
                                          Mean.UMIs.per.Barcode = mean(patient.lesion.rows$Mean.UMIs.per.Barcode),
                                          Number.of.Genes = mean(patient.lesion.rows$Number.of.Genes),
                                          Mean.Genes.per.Barcode = mean(patient.lesion.rows$Mean.Genes.per.Barcode))
                

                merged.sequencing.summary <- bind_rows(merged.sequencing.summary, merged.data)
                
              }
              
              if (nrow(patient.lesion.rows) > 2) {
                
                patient.lesion.date.rows <- patient.lesion.rows[which(str_detect(patient.lesion.rows$sample, pattern = sequence_date)),]  
                
                merged.data <- data.frame(sample = str_c(c(sample.patient.lesion, sequence_date), collapse = "_"),
                                          Total.Reads = sum(patient.lesion.date.rows$Total.Reads),
                                          Valid.Reads = sum(patient.lesion.date.rows$Valid.Reads),
                                          Multimapped.Reads.Aligned..STAR. = sum(patient.lesion.date.rows$Multimapped.Reads.Aligned..STAR.),
                                          Unique.Reads.Aligned..STAR. = mean(patient.lesion.date.rows$Unique.Reads.Aligned..STAR.),
                                          Total.aligned.UMIs..ESAT. = mean(patient.lesion.date.rows$Total.aligned.UMIs..ESAT.),
                                          Total.deduped.UMIs..ESAT. = mean(patient.lesion.date.rows$Total.deduped.UMIs..ESAT.),
                                          Duplication.Rate = mean(patient.lesion.date.rows$Duplication.Rate),
                                          Number.of.Barcodes = mean(patient.lesion.date.rows$Number.of.Barcodes),
                                          Mean.UMIs.per.Barcode = mean(patient.lesion.date.rows$Mean.UMIs.per.Barcode),
                                          Number.of.Genes = mean(patient.lesion.date.rows$Number.of.Genes),
                                          Mean.Genes.per.Barcode = mean(patient.lesion.date.rows$Mean.Genes.per.Barcode))
                
                
                merged.sequencing.summary <- bind_rows(merged.sequencing.summary, merged.data)
                
              }
              
              
            }
          }
  
  
          merged.sequencing.summary <- full_join(merged.file.dims, merged.sequencing.summary, by = c("Sample" = "sample"))
  
          merged.sequencing.summary <- merged.sequencing.summary %>% mutate(barcode.percent.pass.filter = round((passable_cells/Number.of.Barcodes)*100, digits = 3))
          
  
  ggplot(merged.sequencing.summary, aes(x = passable_cells, y = genes, color = Sequence_date)) + 
              geom_point(size = 2, alpha = 0.7) + 
              theme_classic() +
              labs(title = "Sample Sequencing Stats") +
              theme(plot.title = element_text(hjust = 0.5))
  
  ggplot(merged.sequencing.summary, aes(x = passable_cells, y = genes, color = Lesion)) + 
              geom_point(size = 2, alpha = 0.7) + 
              theme_classic() +
              labs(title = "Sample Sequencing Stats") +
              theme(plot.title = element_text(hjust = 0.5))
  
  ggplot(merged.sequencing.summary, aes(x = passable_cells, y = genes, color = Sequence_date)) + 
              geom_point(size = 2, alpha = 0.7) + 
              theme_classic() +
              labs(title = "Sample Sequencing Stats") +
              theme(plot.title = element_text(hjust = 0.5)) + 
              facet_wrap(.~Lesion)
  
  lm <- lm(passable_cells ~ Blister_cell_count + Blister, data = merged.sequencing.summary)
  
  ggplot(merged.sequencing.summary, aes(x = Blister_cell_count, y = passable_cells, color = Sequence_date)) + 
    geom_point(size = 2, alpha = 0.7) + 
    theme_classic() +
    labs(title = "Sample Sequencing Stats", y = "Passable Cells after Sequencing", x = "Initial Blister Cell Count") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0,3000), breaks = seq(0, 3000, by = 500)) +
    scale_x_continuous(limits = c(0, 400000), breaks = seq(0,400000, by = 100000), labels = scales::comma) +
    geom_abline(intercept = as.numeric(lm$coefficients[1]), slope = as.numeric(lm$coefficients[2]), linetype = "dashed")  


  ggplot(merged.sequencing.summary, aes(x = passable_cells, y = Duplication.Rate, color = Sequence_date)) + 
    geom_point(size = 2, alpha = 0.7) + 
    theme_classic() +
    labs(title = "Sample Sequencing Stats", y = "Duplication Rate", x = "Passable Cells after Sequencing") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggplot(merged.sequencing.summary, aes(x = Number.of.Barcodes, y = barcode.percent.pass.filter, color = Sequence_date)) + 
    geom_point(size = 2, alpha = 0.7) + 
    theme_classic() +
    labs(title = "Sample Sequencing Stats", y = "Droplet Pass Filter Perecent", x = "Droplet Number") +
    theme(plot.title = element_text(hjust = 0.5))


##########################
# Data Summary Stat Plots:  
##########################
    
  data.summary.stats$Cell_ID_mutated <- str_replace_all(data.summary.stats$Cell_ID, pattern = "_[:upper:]*$", replacement = "")
  
  sample.summary.stats <- data.summary.stats[0,]
    
  for (s in 1:length(unique(data.summary.stats$Cell_ID_mutated))) {
    
        sample.rows <- data.summary.stats %>% filter(Cell_ID_mutated == unique(data.summary.stats$Cell_ID_mutated)[s])
        
        if (nrow(sample.rows) <= 1000) {
          
          sample.summary.stats <- bind_rows(sample.summary.stats, sample.rows)
          rm(sample.rows)
        }
        
        if (nrow(sample.rows) > 1000) {
      
          sample.rows <- sample.rows[sample(c(1:nrow(sample.rows)), size = 1000, replace = FALSE),]
          
          sample.summary.stats <- bind_rows(sample.summary.stats, sample.rows)
          rm(sample.rows)
        }
    
  }
  
  
  sample.summary.stats$Sequence_date <- factor(sample.summary.stats$Sequence_date, 
                                           levels = c("4_16_20", "1_20_21", "12_3_21", "4_28_22"))
  
  sample.summary.stats$Lesion[which(sample.summary.stats$Lesion == "SLS")] <- "Day2_SLS"
  sample.summary.stats$Lesion <- factor(sample.summary.stats$Lesion, 
                                    levels = c("NL", "Day2_SLS","Acetone", "Day2_SADBE", "Day4_SADBE"))
  

  ggplot(sample.summary.stats, aes(x = UMI.Sum, y = percent.mt, color = Sequence_date)) + 
    geom_point(alpha = 0.7) + 
    theme_classic() +
    labs(title = "Sample Sequencing Stats", y = "Percent Mitochondrial", x = "Droplet UMI sum") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Cell_ID_mutated, scales = "free") + 
    geom_vline(xintercept = 300, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 7.5, linetype = "dashed", color = "red") 
    
  ggplot(sample.summary.stats, aes(x = UMI.Sum, y = percent.mt, color = Sequence_date)) + 
    geom_point(alpha = 0.7) + 
    theme_classic() +
    labs(title = "Sample Sequencing Stats", y = "Percent Mitochondrial", x = "Droplet UMI sum") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Sequence_date, scales = "free") + 
    geom_vline(xintercept = 300, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 7.5, linetype = "dashed", color = "red") 
  
  ggplot(sample.summary.stats, aes(x = UMI.Sum, y = percent.mt)) + 
    geom_point(alpha = 0.7) + 
    theme_classic() +
    labs(title = "Cell Sequencing Stats", y = "Percent Mitochondrial", "Droplet UMI sum") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = 300, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 7.5, linetype = "dashed", color = "red") 
  
  #UMI density plots:
  ggplot(filter(sample.summary.stats, UMI.Sum > 10), aes(x = log10(UMI.Sum), color = Sequence_date)) + 
    geom_density(alpha = 0.7) + 
    theme_classic() +
    labs(title = "Sequencing Stats", x = "Log10 (Droplet UMI sum)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = log10(300), linetype = "dashed", color = "red")
  
  ggplot(sample.summary.stats, aes(x = log10(UMI.Sum), color = Sequence_date)) + 
      geom_density(alpha = 0.7) + 
      theme_classic() +
      labs(title = "Sequencing Stats", x = "Log10 (Droplet UMI sum)") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(xintercept = log10(300), linetype = "dashed", color = "red")
    
            
  
