###########################################################################################################
### Data Pre-Filtering

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


#Filter Parameters:
percent.mt_cutoff <- 7.5
UMI_cutoff <- 300
gene_cutoff <- 150

#Save Folder Location:
save.path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/Filtered_raw_data/"


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
                    
                    #Filter cells:
                    sample.tmp <-  sample.tmp %>% filter(UMI.Sum > as.numeric(UMI_cutoff) & percent.mt < as.numeric(percent.mt_cutoff) & gene.Sum > as.numeric(gene_cutoff))
                    
                    #Save as tsv:
                    write.table(sample.tmp, file = paste(save.path, sample, ".tsv", sep = ""), row.names=FALSE, sep="\t")
                    
                    rm(sample.tmp)
              }
            
          }

      
#All other Skin Samples:

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
                            
                            #Filter cells:
                            sample.tmp <-  sample.tmp %>% filter(UMI.Sum > as.numeric(UMI_cutoff) & percent.mt < as.numeric(percent.mt_cutoff) & gene.Sum > as.numeric(gene_cutoff))
                            
                            #Save as tsv:
                            write.table(sample.tmp, file = paste(save.path, sample, ".tsv", sep = ""), row.names=FALSE, sep="\t")
                            
                            rm(sample.tmp)
                          }
                    }
          }