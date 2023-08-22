
###########################################################################################################
### Data Import

###########################################################################################################
setwd("/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

library(dplyr)
library(data.table)
library(purrr)
library(tibble)
library(stringr)

#Save File dimensions as summary dataframe:
file_dims_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(file_dims_df) <- c("Sample", "rows", "columns")
class(file_dims_df$Sample) <- "character"
class(file_dims_df$rows) <- "numeric"
class(file_dims_df$columns) <- "numeric"


#Get File List of pre-filtered raw data:
path = "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/Filtered_raw_data/"

filelist <- list.files(path = path, pattern = ".tsv")

sample_list <- c()

    #Import file data:
      for (f in 1:length(filelist)) {
        
        file.data <- read.delim(file = paste(path, filelist[f], sep = ""))
        
        if ((nrow(file.data) > 0) == TRUE) {
          
            file.data$Cell_ID <- as.character(file.data$Cell_ID)
            file.data$Sample <- as.character(file.data$Sample)
            
            sample_name <- str_replace(filelist[f], pattern = ".tsv", replacement = "")
            
            file.dims <- data.frame(Sample = sample_name,
                                    rows = nrow(file.data),
                                    columns = ncol(file.data))
            file_dims_df <- bind_rows(file_dims_df, file.dims)
            
            assign(sample_name, file.data)
            
            sample_list <- c(sample_list, sample_name)
            
            rm(sample_name)
        }
        
        rm(file.data)
      }

#Save dimension stats:
save_path <- "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/"
write.table(file_dims_df, file = paste(save_path, "file_dims_5_6_22.tsv", sep = ""), row.names=FALSE, sep="\t")



#Merge sample data to a single master_data file:
  master_data <- reduce(mget(sample_list), full_join)
  print(paste("Merge Complete, Master Data Size is...", dim(master_data)[1], "rows   x ", dim(master_data)[2], "columns"))
  
  
  #Convert NAs to 0:
  master_data <- as.data.frame(master_data)
  master_data[is.na(master_data)] <- 0  

  #Save:
  save_path <- "/project/umw_john_harris/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/Data/"
  
  write.table(master_data, file = paste(save_path, "master_data_5_6_22.tsv", sep = ""), row.names=FALSE, sep="\t")
  
  
#Summary Table:
summary_table <- master_data %>% 
  group_by(Sample) %>%
  summarize(Number_cells = n(),
            Mean_UMIs = mean(UMI.Sum),
            Median_UMIs = median(UMI.Sum),
            Max_Gene_Num = max(gene.Sum),
            Mean_Genes = mean(gene.Sum),
            Median_Genes = median(gene.Sum),
            Mean_Percent.MT = mean(percent.mt))

write.table(summary_table, file = paste(save_path, "Sample_Summary_table_5_6_22.tsv", sep = ""), row.names=FALSE, sep="\t")
