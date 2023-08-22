###############################################################
# Create a small data subset for testing code:
###############################################################
library(dplyr)
library(data.table)
library(purrr)
library(tibble)
library(stringr)
library(SignallingSingleCell)

setwd("/project/umw_john_harris/frisoli/RStudio/Multi_disease_analysis/8_23_21_Analysis/")

#Load Data
master_data <- fread(file = "Data/master_data_8_25_21.tsv", header = TRUE)

#Filter contact derm, GVHD, and SJS out of this analysis:
master_data <- master_data[-which(str_detect(string = master_data$Cell_ID, pattern = "SLS|SADBE|GH|SJS|HC")),]

#Create sample subset with 25,000 cells:
sample <- sample(seq(1:nrow(master_data)), 25000)

sample_data <- master_data[which(rownames(master_data) %in% sample),]

#Save Sample Data:
save_path <- "/project/umw_john_harris/frisoli/RStudio/Multi_disease_analysis/8_23_21_Analysis/Data/"

write.table(sample_data, file = paste(save_path, "sample_data_8_30_21.tsv", sep = ""), row.names=FALSE, sep="\t")