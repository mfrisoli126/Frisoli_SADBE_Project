#################################    
# Load and Merge all Olink Data:
#################################    
#load packages:
library(readxl)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(rstatix)
library(ggrepel)
library(lubridate)
library(umap)
library(ggnewscale) 
library(Seurat)
library(pheatmap)
library(glmnet)

#load data:
setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")
        
        blister.protein <- read_xlsx(path = "Data/Olink/Olink_6_15_22_Plate1_blister.xlsx", sheet = "Sheet1")
        MN.protein <- read_xlsx(path = "Data/Olink/Olink_6_15_22_Plate1_MN.xlsx", sheet = "Sheet1")
        template <- read_xlsx(path = "Data/Olink/samplemanifest_template_MIKE.xlsx")
        template_tube_labels <- read_xlsx(path = "Data/Olink/6_15_22_Olink_Layout.xlsx")
        colnames(template_tube_labels)[which(colnames(template_tube_labels) == "Well_Name")] <- "WellID"
        
        template <- left_join(template, template_tube_labels)
        
        MN.protein2 <- read_xlsx(path = "Data/Olink/9_13_22_olink_NPX.xlsx", sheet = "Sheet1")
        template2 <- read_xlsx(path = "Data/Olink/9_13_22_Olink_template.xlsx")
        
        MN.protein3 <- read_xlsx(path = "Data/Olink/02032023_Inflammation_microneedles_NPX.xlsx", sheet = "Sheet1")
        
        lesion.scores <- read_xlsx(path = "Data/Olink/master_SADBE_sample_info.xlsx")
    
      #Reformat and merge:
        blister.protein$Sample <- str_replace_all(blister.protein$Sample, pattern = "Sample", replacement = "")
        blister.protein$Sample <- as.numeric(blister.protein$Sample)
        
        MN.protein$Sample <- str_replace_all(MN.protein$Sample, pattern = "Sample", replacement = "")
        MN.protein$Sample <- as.numeric(MN.protein$Sample)
        
        template$SampleID <- str_replace_all(template$SampleID, pattern = "Sample", replacement = "")
        template$SampleID <- as.numeric(template$SampleID)
        template$Disease <- "Contact_derm"
        template$sample_num <- paste(template$sample_type, template$tube_num, sep = "_")
        template$Olink_batch <- "Batch1"
        
        colnames(lesion.scores) <- tolower(colnames(lesion.scores))
        
        template2 <- template2 %>% filter(sample_type == "Microneedle")
        template2$Disease <- template2$patient
        template2$Disease[which(str_detect(template2$Disease, pattern = "^HC[:digit:]"))] <- "Contact_derm"
        template2$Disease[which(str_detect(template2$Disease, pattern = "^V"))] <- "Vitiligo"
        template2$Disease[which(is.na(template2$Disease))] <- "Technical_Control"
        
        colnames(template2)[which(colnames(template2) == "sample_num")] <- "SampleID"
        template2$sample_num <- paste(template2$sample_type, template2$tube_num, sep = "_")
        template2$Olink_batch <- "Batch2"
        
        MN.protein3$Olink_batch <- "Batch3"
        
        
        #merge templates with lesion score info:
        template$blister_date <- as.Date(as.numeric(template$blister_date), origin = "1899-12-30")
        template2$blister_date <- as.Date(as.numeric(template2$blister_date), origin = "1899-12-30")
        
        lesion.scores$blister_date <- as_date(lesion.scores$blister_date)
        
        template <- left_join(template, lesion.scores)
        template2 <- left_join(template2, lesion.scores)
        
        blister.protein <- right_join(template, blister.protein, by = c("SampleID" = "Sample"))
        MN.protein <- right_join(template, MN.protein, by = c("SampleID" = "Sample"))
        MN.protein2 <- right_join(template2, MN.protein2, by = c("WellID" = "Well"))
        
        MN.protein <- MN.protein[,-which(colnames(MN.protein) %in% c("SampleID"))]
        MN.protein2 <- MN.protein2[,-which(colnames(MN.protein2) %in% c("SampleID"))]
        
        blister.protein <- blister.protein[,which(colnames(blister.protein) %in% colnames(MN.protein2))]
        MN.protein <- MN.protein[,which(colnames(MN.protein) %in% colnames(MN.protein2))]
        MN.protein2$tube_num <- as.character(MN.protein2$tube_num)
        
        olink.data <- bind_rows(blister.protein, MN.protein, MN.protein2)
        
        #Add Wei's Microneedle data:
        template3 <- as.data.frame(matrix(ncol=length(colnames(olink.data)[-which(colnames(olink.data) %in% colnames(MN.protein3))]), nrow = nrow(MN.protein3)))
        colnames(template3) <- colnames(olink.data)[-which(colnames(olink.data) %in% colnames(MN.protein3))]
        
        MN.protein3 <- cbind(template3, MN.protein3)
        MN.protein3$sample_type <- "Microneedle"
        MN.protein3$Disease <- "Contact_derm" 
        
        #Filter out "IPC" samples:
        MN.protein3 <- MN.protein3[-which(str_detect(MN.protein3$sample_num, pattern = "IPC")),]
        MN.protein3 <- MN.protein3[-which(MN.protein3$sample_num %in% c("S12 Sample control", "S24 Sample control")),]
        
        #Annotate Wei's samples:
        wei.samples <- MN.protein3$sample_num[-which(str_detect(MN.protein3$sample_num, pattern = "mf"))]
        
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "blistfluid")),"sample_type"] <- "Natural Blister"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "ctrl")),"lesion"] <- "NL"
        
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "Neg control")),"lesion"] <- "Never in skin"
        
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "bb00")),"patient"] <- "bb00"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "bb06")),"patient"] <- "bb06"
        
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "48h-[:digit:]")),"lesion"] <- "Day 2 SADBE"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "48h$")),"lesion"] <- "Day 2 SADBE"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "48h-mem")),"lesion"] <- "Day 2 SADBE"
        
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "96h-[:digit:]")),"lesion"] <- "Day 4 SADBE"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "96h$")),"lesion"] <- "Day 4 SADBE"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "96h-mem")),"lesion"] <- "Day 4 SADBE"
        
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "carbamix")),"lesion"] <- "carbamix"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "nickel")),"lesion"] <- "nickel"
        MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "cobalt")),"lesion"] <- "cobalt"
        
        wei.sadbe.samples <- MN.protein3[which(str_detect(MN.protein3$sample_num, pattern = "bb[:digit:]*-[:digit:]*h")),]
        wei.sadbe.samples <- wei.sadbe.samples[,which(colnames(wei.sadbe.samples) %in% c("sample_num", "lesion visual score"))]
        wei.sadbe.samples$dose_label <- sapply(wei.sadbe.samples$sample_num, FUN = function (x) as.numeric(unlist(str_split(x, pattern = "-"))[3]))
        wei.sadbe.samples$timepoint <- sapply(wei.sadbe.samples$sample_num, FUN = function (x) unlist(str_split(x, pattern = "-"))[2])
        wei.sadbe.samples$Patient <- sapply(wei.sadbe.samples$sample_num, FUN = function (x) unlist(str_split(x, pattern = "-"))[1])
        wei.sadbe.samples$Patient <- str_replace_all(wei.sadbe.samples$Patient, pattern = "^S[:digit:]* ", replacement = "")
        wei.sadbe.samples <- wei.sadbe.samples[-which(is.na(wei.sadbe.samples$dose_label)),]
        
        wei.sadbe.samples <- wei.sadbe.samples %>% arrange(Patient, timepoint, dose_label)
        
        wei.sadbe.samples[which(wei.sadbe.samples$Patient == "bb00"),"lesion visual score"] <- c(1,1,rep(2, times = 7),
                                                                                                 2, rep(3, times = 9))
        wei.sadbe.samples[which(wei.sadbe.samples$Patient == "bb06"),"lesion visual score"] <- c(1,1,1,1,1,1,1,1,2,2,
                                                                                                 1,2,2,2,1,2,2,2,2,2)
        
        for (r in 1:nrow(MN.protein3)) {
          
          patient <- MN.protein3$patient[r]
          sam_num <- MN.protein3$sample_num[r]
          
          if (is.na(patient)) {
            next()
          }
          
          if (patient %in% c("bb00", "bb06")) {
            dl <- sapply(sam_num, FUN = function (x) unlist(str_split(x, pattern = "-"))[3])
            tp <- sapply(sam_num, FUN = function (x) unlist(str_split(x, pattern = "-"))[2])
            
            if (dl %in% c("mem", "ctrl") | !tp %in% c("48h", "96h")) {
              next()
            }
            
            score <- wei.sadbe.samples %>% filter(timepoint == tp & dose_label == dl & Patient == patient) %>% dplyr::pull(`lesion visual score`)
            
            MN.protein3$`lesion visual score`[r] <- score
            
            rm(dl)
            rm(tp)
          }
          
        }
        
        
        olink.data <- bind_rows(olink.data, MN.protein3)
        olink.data <- olink.data[,-which(colnames(olink.data) %in% c("Inc Ctrl 1","Ext Ctrl","Inc Ctrl 2","Det Ctrl"))]
        
        olink.data <- olink.data %>% mutate(Lesion_severity_score = case_when(`lesion visual score` == 3 ~ "+++",
                                                                              `lesion visual score` == 2 ~ "++",
                                                                              `lesion visual score` %in% c(0.5, 1) ~ "+",
                                                                              `lesion visual score` %in% c(0.25) ~ "?/+",
                                                                              `lesion visual score` == 0 ~ "-"))
        
        olink.data$Lesion_severity_score <- factor(olink.data$Lesion_severity_score,
                                                   levels = c("-", "?/+", "+", "++", "+++"))
        
        olink.data <- olink.data[,c(1:which(colnames(olink.data) == "lesion visual score"), which(colnames(olink.data) == "Lesion_severity_score"), c(c(which(colnames(olink.data) == "lesion visual score")+1):c(ncol(olink.data)-1)))]
        
        #Add sample info for Mike's samples that were re-analyzed to test concentration columns:
        mike.concentration.test.samples <- olink.data$sample_num[which(str_detect(olink.data$sample_num, pattern = "mf-[:digit:]"))]
        
        mike.concentration.test.sample.layout <- data.frame(sample = seq(1:6),
                                                            batch = c(2,1,1,2,2,2),
                                                            well = c("F1", "D9", "C9", "H3", "D5", "C3"))
        
        
        for (s in 1:length(mike.concentration.test.samples)) {
          
          sample_num_label <- mike.concentration.test.samples[s]
          sample_num <- str_extract(sample_num_label, pattern = "-[:digit:]$")
          sample_num <- as.numeric(str_replace(sample_num, pattern = "-", replacement = ""))
          
          sample.layout.info <- mike.concentration.test.sample.layout %>% filter(sample == sample_num)
          
          sample.olink.row.data <- olink.data %>% filter(Olink_batch == paste("Batch", sample.layout.info$batch, sep = "") & WellID == sample.layout.info$well)
          sample.olink.meta.info <- sample.olink.row.data[,c(1:which(colnames(sample.olink.row.data) == "IL8")-1)]
          
          
          sample.olink.meta.info$WellID <- NA
          sample.olink.meta.info$sample_num <- sample_num_label
          sample.olink.meta.info$Olink_batch <- "Batch3"
          
          #Fill into Olink data:
          olink.data[which(olink.data$sample_num == sample_num_label),which(colnames(olink.data) %in% colnames(sample.olink.meta.info))] <- sample.olink.meta.info
          
        }
        
        
        olink.data <- olink.data %>% filter(sample_type %in% c("Blister", "Microneedle", "Natural Blister"))
        
        
        #Add column for microneedle dilution method:
        olink.data <- olink.data %>% mutate(dilution_method = case_when(sample_type %in% c("Blister", "Natural Blister") ~ "Neat",
                                                                        sample_type == "Microneedle" & Olink_batch == "Batch3" & str_detect(sample_num, pattern = "dilute") == TRUE ~ "Old method",
                                                                        sample_type == "Microneedle" & Olink_batch != "Batch3"  ~ "Old method",
                                                                        TRUE ~ "New method"))
        
        olink.data <- olink.data[,c(c(1:which(colnames(olink.data) == "Olink_batch")), which(colnames(olink.data) == "dilution_method"), c(c(which(colnames(olink.data) == "Olink_batch")+1):c(length(colnames(olink.data))-1)))]
        
        olink.data <- olink.data[,-which(colnames(olink.data) %in% c("WellID", "blister_bate","SampleID", "PlateID","freezer_box", "Plate ID", "QC Warning", "bioanalyzer", "lesion photo", "patient note"))]
        
        olink.data <- pivot_longer(olink.data, cols = seq(which(str_detect(colnames(olink.data), pattern = "IL8")), ncol(olink.data)), names_to = "Protein", values_to = "Concentration")
        
        olink.data$lesion_proximity <- olink.data$lesion
        
        olink.data$lesion_proximity[which(olink.data$lesion_proximity %in% c("Day 4 upper arm rxn", "Day 2 torso rxn"))] <- "distant"
        olink.data$lesion_proximity[-which(olink.data$lesion_proximity == "distant")] <- "local"
        
        olink.data$lesion[which(olink.data$lesion == "Day 4 upper arm rxn")] <- "Day4 Allergy"
        olink.data$lesion[which(olink.data$lesion == "Day 2 torso rxn")] <- "Day2 Allergy"
        
        olink.data$lesion[which(olink.data$lesion == "Day 4 SADBE")] <- "Day4 Allergy"
        olink.data$lesion[which(olink.data$lesion == "Day 2 SADBE")] <- "Day2 Allergy"
        olink.data$lesion[which(olink.data$lesion == "Day 2 SLS")] <- "Irritant"
        olink.data$lesion[which(olink.data$lesion == "Day 2 Acetone")] <- "Acetone Vehicle"
        olink.data$lesion[which(olink.data$lesion == "NL")] <- "Nonlesional"
        olink.data$lesion[which(str_detect(olink.data$lesion, pattern = "Never in skin"))] <- "Negative Control"
        olink.data$lesion[which(str_detect(olink.data$lesion, pattern = "^L[:digit:]"))] <- "Vitiligo_Lesion"
        olink.data$lesion[which(olink.data$lesion == "L")] <- "Vitiligo_Lesion"
        
        viitligo.olink <- olink.data %>% filter(Disease == "Vitiligo")
        
        olink.data <- olink.data %>% filter(!Disease == "Vitiligo")
        
        olink.data$lesion <- factor(olink.data$lesion,
                                    levels = c("Negative Control", "Nonlesional", "Irritant", "Acetone Vehicle", "Day2 Allergy" , "Day4 Allergy", "nickel", "cobalt", "carbamix"))
        
        viitligo.olink$lesion <- factor(viitligo.olink$lesion,
                                        levels = c("Nonlesional", "Vitiligo_Lesion"))
        
        olink.data$`lesion visual score`[which(is.na(olink.data$`lesion visual score`))] <- 0
        olink.data$`lesion visual score` <- factor(olink.data$`lesion visual score`,
                                                   levels = c(0, 0.25, 0.5, 1, 2))
        
        #olink.data <- olink.data %>% filter(Protein != "FGF-23")
        
        #convert back to raw data from log2 values:
        olink.data$Concentration <- 2^(olink.data$Concentration)
        
        #Rename IFNG:
        olink.data$Protein[which(olink.data$Protein == "IFN-gamma")] <- "IFNG"
        
        
      #Filter out 2 positive acetone reactions:
        olink.data <- olink.data %>% mutate(DE.inclusion.factor = case_when(lesion != "Acetone Vehicle" ~ "yes",
                                                                            lesion == "Acetone Vehicle" & as.numeric(as.character(`lesion visual score`)) < 1 ~ "yes",
                                                                            lesion == "Acetone Vehicle" & as.numeric(as.character(`lesion visual score`)) >= 1 ~ "no"))
        
        olink.data$DE.inclusion.factor[which(olink.data$Olink_batch == "Batch3")] <- "yes"
        
        filtered.olink.data <- olink.data %>% filter(DE.inclusion.factor == "yes")
        
        filtered.olink.data <- filtered.olink.data[-which(filtered.olink.data$lesion %in% c("nickel", "cobalt", "carbamix")),]
        
        #Filter out Wei's samples testing memory placement in same location:
        samples.to.filter <- unique(filtered.olink.data$sample_num[which(str_detect(filtered.olink.data$sample_num, pattern = "-ctrl|-mem"))])
        filtered.olink.data <- filtered.olink.data %>% filter(!sample_num %in%samples.to.filter)
        
        
        #save(filtered.olink.data, file = "Data/Olink/filtered.olink.data.Rdata")
        
#################################   
# Analyze Olink Data:
#################################
#load packages:
library(readxl)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(rstatix)
library(ggrepel)
library(lubridate)
library(umap)
library(ggnewscale) 
library(Seurat)
library(pheatmap)
library(glmnet)
library(purrr)
library(writexl)
        
    #load data:
    setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")
    
    load("Data/Olink/filtered.olink.data.Rdata")    
        
          #filter out attempted concentrated samples since it didn't even do anything:
              filtered.olink.data <- filtered.olink.data[-which(str_detect(filtered.olink.data$sample_num, pattern = "mf-")),]
          

    #DE analysis:
              #Aggregate for Fold change calculation:
              olink.blister.aggregate <- filtered.olink.data %>% 
                filter(sample_type == "Blister") %>%
                group_by(lesion, Protein) %>%
                summarise(mean(Concentration)) %>%
                rename_at("mean(Concentration)", ~ "protein_value")
              
              olink.blister.aggregate <- olink.blister.aggregate %>% spread(key = lesion, value = protein_value)
              
              #Set up DE comparisons:           
              comparisons <- data.frame(reference = c("Nonlesional", "Nonlesional",   "Nonlesional",   "Nonlesional",   "Irritant",     "Irritant"),
                                        comparison = c("Irritant", "Acetone Vehicle", "Day2 Allergy", "Day4 Allergy", "Day2 Allergy", "Day4 Allergy"))
              
              for (c in 1:nrow(comparisons)) {
                
                ref = comparisons$reference[c]
                comp = comparisons$comparison[c]
                
                #Paired T tests for significance:
                comparison.blister.DE <- filtered.olink.data %>% filter(sample_type == "Blister" & lesion %in% c(ref, comp)) %>% group_by(Protein) %>% pairwise_t_test(
                  Concentration ~ lesion, paired = FALSE, 
                  p.adjust.method = "BH") %>%
                  select(-p, -p.signif) # Remove details
                
                #Fold Change:
                comparison.blister.aggregate <- olink.blister.aggregate %>% mutate(foldchange = eval(as.symbol(comp))/eval(as.symbol(ref)))
                
                comparison.blister.DE <- cbind(comparison.blister.DE, comparison.blister.aggregate$foldchange[match(comparison.blister.DE$Protein, comparison.blister.aggregate$Protein)])
                
                colnames(comparison.blister.DE)[9] <- "Protein_fold_change"
                
                comparison.blister.DE <- comparison.blister.DE %>% filter(!(is.nan(p.adj)))
                comparison.blister.DE$p.adj.signif <- factor(comparison.blister.DE$p.adj.signif, levels = c("ns", "*", "**", "***", "****"))
                
                comparison.blister.DE <- comparison.blister.DE %>% mutate(change_direction = ifelse(Protein_fold_change > 1, "Up", "Down"))
                comparison.blister.DE <- comparison.blister.DE %>% mutate(Color_factor = case_when(change_direction == "Up" & p.adj <= 0.05 ~ "Significantly Up", 
                                                                                                   change_direction == "Down" & p.adj <= 0.05 ~ "Significantly Down",
                                                                                                   p.adj > 0.05 ~ "Not Significant"))
                comparison.blister.DE$Color_factor <- factor(comparison.blister.DE$Color_factor, levels = c("Significantly Up", "Significantly Down", "Not Significant"))
                
                
                #Label comparison conditions:
                comparison.blister.DE$Reference <- ref
                comparison.blister.DE$Comparison <- comp
                
                if (c == 1) {
                  olink.blister.DE <- comparison.blister.DE
                }
                
                if (c > 1) {
                  olink.blister.DE <- bind_rows(olink.blister.DE, comparison.blister.DE)
                }
                
              }
              
              olink.blister.DE$ref_comp <- paste(olink.blister.DE$Reference, "vs",olink.blister.DE$Comparison)
              olink.blister.DE$ref_comp <- factor(olink.blister.DE$ref_comp,
                                                  levels = c("Nonlesional vs Irritant", "Nonlesional vs Acetone Vehicle", "Nonlesional vs Day2 Allergy", "Nonlesional vs Day4 Allergy", "Irritant vs Day2 Allergy", "Irritant vs Day4 Allergy"))
              
              olink.blister.DE %>% filter(!ref_comp %in% c("Irritant vs Day2 Allergy", "Irritant vs Day4 Allergy")) %>%
                ggplot(., aes(x = log10(Protein_fold_change), y = -log10(p.adj), color = Color_factor)) + 
                facet_wrap(~ref_comp) +
                geom_point(size = 2) + 
                theme_classic() + 
                labs(title = "Suction Blister Proteins", x = "log10(Fold Change)", color = "Protein Regulation") + 
                theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = c("red", "blue", "gray")) +
                scale_x_continuous(limits = c(-1.5,1.5), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5)) +
                geom_text_repel(data = . %>% filter(p.adj.signif != "ns"), 
                                mapping = aes(x = log10(Protein_fold_change), 
                                              y = -log10(p.adj), label = Protein), size = 3, color = "black") 
              
              olink.blister.DE %>% filter(ref_comp %in% c("Irritant vs Day2 Allergy", "Irritant vs Day4 Allergy")) %>%
                ggplot(., aes(x = log10(Protein_fold_change), y = -log10(p.adj), color = Color_factor)) + 
                facet_wrap(~ref_comp) +
                geom_point(size = 2) + 
                theme_classic() + 
                labs(title = "Suction Blister Proteins", x = "log10(Fold Change)", color = "Protein Regulation") + 
                theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = c("red", "blue", "gray")) +
                scale_x_continuous(limits = c(-1.5,2), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5,2.0)) + 
                geom_text_repel(data = . %>% filter(p.adj.signif != "ns"), 
                                mapping = aes(x = log10(Protein_fold_change), 
                                              y = -log10(p.adj), label = Protein), size = 3, color = "black") 
              
              colnames(olink.blister.DE)
              olink.blister.DE %>% filter(ref_comp %in% c("Irritant vs Day2 Allergy")) %>%
                ggplot(., aes(x = log10(Protein_fold_change), y = -log10(p.adj), color = Color_factor)) + 
                facet_wrap(~ref_comp) +
                geom_point(size = 2) + 
                theme_classic() + 
                labs(title = "Suction Blister Proteins", x = "log10(Fold Change)", color = "Protein Regulation") + 
                theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = c("red", "blue", "gray")) +
                scale_x_continuous(limits = c(-1,2), breaks = c(-1,-0.5,0,0.5,1,1.5,2.0)) + 
                geom_text_repel(data = . %>% filter(ref_comp %in% c("Irritant vs Day2 Allergy") & Color_factor != "Not Significant"), 
                                mapping = aes(x = log10(Protein_fold_change), 
                                              y = -log10(p.adj), label = Protein), size = 3, color = "black") 
              
              blister.Olink.volcano.plot <- last_plot()

              #ggsave(blister.Olink.volcano.plot, width = 20, height = 15,units = "cm",filename = "Plots/olink/DE_volcano.pdf", limitsize = FALSE)  
              
    #Heatmaps:

              #Identify proteins reliably captured above convincing threshold:
              protein.stats <- filtered.olink.data %>% group_by(Protein)%>% summarize(min = min(Concentration),
                                                                                      max = max(Concentration),
                                                                                      log10max = log10(max(Concentration)),
                                                                                      average = mean(Concentration),
                                                                                      log10range = log10(max(Concentration))-log10(min(Concentration)))
              
              
              
              lesion.protein.stats <- filtered.olink.data %>% group_by(Protein, lesion)%>% summarize(min = min(Concentration),
                                                                                                     max = max(Concentration),
                                                                                                     log10max = log10(max(Concentration)),
                                                                                                     average = mean(Concentration),
                                                                                                     log10range = log10(max(Concentration))-log10(min(Concentration)))
              
              
              protein.stats$nl_vs_control_log10max_diff <- NA
              protein.stats$allergy_vs_control_log10max_diff <- NA
              
              for (p in 1:length(protein.stats$Protein)) {
                
                prot = protein.stats$Protein[p]
                
                prot.les.stats <- lesion.protein.stats %>% filter(Protein ==prot)
                
                neg.control.max <- prot.les.stats %>% filter(lesion == "Negative Control") %>% dplyr::pull(log10max)
                nonlesional.max <- prot.les.stats %>% filter(lesion == "Nonlesional") %>% dplyr::pull(log10max)
                allergy.max <- prot.les.stats %>% filter(lesion %in% c("Day2 Allergy", "Day4 Allergy")) %>% dplyr::pull(log10max) %>% mean()
                
                nl.diff = nonlesional.max-neg.control.max
                allergy.diff = allergy.max-neg.control.max
                
                protein.stats[p, "nl_vs_control_log10max_diff"] <-nl.diff
                protein.stats[p, "allergy_vs_control_log10max_diff"] <-allergy.diff
                
              }
              
              Proteins.reliably.captured <- protein.stats %>% filter(log10max >=2 & allergy_vs_control_log10max_diff>=1) %>% dplyr::pull(Protein)
              Proteins.reliably.captured <- unique(c(Proteins.reliably.captured, "IL5", "IL4", "IL13", "TNFS14", "TNFB", "TRAIL"))
              
              positive.olink.data <- filtered.olink.data %>% filter(Protein %in% Proteins.reliably.captured)
              
              positive.olink.data.new_dilution_method <- positive.olink.data %>% filter(dilution_method != "Old method")
              
              merged_data_spread <- positive.olink.data %>% pivot_wider(names_from = Protein, values_from = Concentration)
              
              merged.counts = log2(t(merged_data_spread[,c(which(colnames(merged_data_spread) == "IL8"):length(colnames(merged_data_spread)))]))
              colnames(merged.counts) <- merged_data_spread$sample_num
              
              Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
              
              merged.pheat.annotations <- column_to_rownames(merged_data_spread[,which(colnames(merged_data_spread) %in% c("sample_num", "lesion","sample_type", "dilution_method"))], var = "sample_num")
              merged_colour = list(
                sample_type = c(Blister = "#F85656", `Natural Blister` = "#FF0000" , Microneedle = "#5B4BFF"),
                dilution_method = c(Neat = "#FD5DE2", `New method` = "#A459FF" , `Old method` = "#BBB46D"),
                lesion = c(`Negative Control` = "#000000", Nonlesional = "#9C9C9C", Irritant = "#FFB659", `Acetone Vehicle` = "#9ABFBA", `Day2 Allergy` = "#52B6FF", `Day4 Allergy` = "#4840FF")
              )
              
              # All data:
              olink.heatmap <- pheatmap(merged.counts,
                                        annotation_col = merged.pheat.annotations,
                                        annotation_colors = merged_colour,
                                        cutree_cols = 1,
                                        cutree_rows = 1,
                                        color = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                        border_color=NA,
                                        treeheight_row = 0,
                                        treeheight_col = 0)
              
                
              #Just Suction Blister data:
                    blister.data <- positive.olink.data %>% filter(sample_type == "Blister")
                    blister_data_spread <- blister.data %>% pivot_wider(names_from = Protein, values_from = Concentration)
                    
                    #Convert back to log2 NPX values:
                    blister.counts = log2(t(blister_data_spread[,c(which(colnames(blister_data_spread) == "IL8"):length(colnames(blister_data_spread)))]))
                    colnames(blister.counts) <- blister_data_spread$sample_num
                    
                    meta.data <- blister_data_spread[,c(1:which(colnames(blister_data_spread) == "IL8")-1)]
                    
                    
                    blister.pheat.annotations <- column_to_rownames(blister_data_spread[,which(colnames(blister_data_spread) %in% c("sample_num","lesion", "Lesion_severity_score"))], var = "sample_num")
                    
                    colnames(blister.pheat.annotations)[which(colnames(blister.pheat.annotations) == "Lesion_severity_score")] <- "Lesion_severity_score"
                    

                    color.palette <- c("#080808","#F50000")
                    GetPalette = colorRampPalette(color.palette)
                    
                    ColorCount = as.numeric(length(unique(blister.pheat.annotations$Lesion_severity_score)))
                    custom.colors <- GetPalette(ColorCount)
                    
                    my_colour = list(
                      lesion = c(Nonlesional = "#9C9C9C", Irritant = "#FFB659", `Acetone Vehicle` = "#9ABFBA", `Day2 Allergy` = "#52B6FF", `Day4 Allergy` = "#4840FF"),
                      Lesion_severity_score = c(`-` = "#080808", `?/+` = "#430606", `+` = "#7E0404", `++` = "#B90202", `+++` = "#F50000")
                    )
                    
                    set.seed(100)
                    olink.heatmap <- pheatmap(blister.counts,
                                              annotation_col = blister.pheat.annotations,
                                              annotation_colors = my_colour,
                                              cutree_cols = 1,
                                              cutree_rows = 1,
                                              color = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                              border_color=NA,
                                              treeheight_row = 0,
                                              treeheight_col = 0)
                    
                    #reorder and filter to reliable proteins:
                    reordered.proteins <- c("MCP-1", "MCP-4", "MMP-1",
                                            "uPA", "CD40", "HGF", "4E-BP1", "VEGFA",
                                            "Flt3L", "IL18",
                                            "CCL23","CSF-1",
                                            "IL8", "CXCL6", "CCL3", "CCL4",
                                            "CCL19", "CCL20", "OPG","MCP-2","CXCL1","IL6", "OSM",
                                            "MMP-10","CXCL10", "CXCL9", "CXCL11", "IFNG",
                                            "TRAIL", "TNFB", "TNFSF14", "IL4", "IL5", "IL13")
                    
                    
                    blister.olink.heatmap <- pheatmap(blister.counts[reordered.proteins,],
                                                      annotation_col = blister.pheat.annotations,
                                                      annotation_colors = my_colour,
                                                      cutree_cols = 1,
                                                      cutree_rows = 1,
                                                      color = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                      border_color=NA,
                                                      treeheight_row = 0,
                                                      treeheight_col = 0,
                                                      show_colnames = F,
                                                      cluster_rows = F)
                    
                    
              #Just Microneedle data that looks good:
                    positive.MN.data <- positive.olink.data %>% filter(dilution_method != "Old method" & sample_type == "Microneedle")
                    sample_stats <- positive.MN.data %>% group_by(sample_num) %>% summarize(average_concentration = mean(Concentration),
                                                                                            max_concentration = max(Concentration))
                    
                    bad.samples <- sample_stats %>% filter(max_concentration <1000) %>% dplyr::pull(sample_num)
                    
                    merged_MN_data_spread <- positive.MN.data %>% filter(!sample_num %in% bad.samples) %>% pivot_wider(names_from = Protein, values_from = Concentration)
                                  
                    #Convert back to log2 NPX values:
                    MN.counts = log2(t(merged_MN_data_spread[,c(which(colnames(merged_MN_data_spread) == "IL8"):length(colnames(merged_MN_data_spread)))]))
                    colnames(MN.counts) <- merged_MN_data_spread$sample_num
                    
                    meta.data <- merged_MN_data_spread[,c(1:which(colnames(merged_MN_data_spread) == "IL8")-1)]
                    
                    
                    MN.pheat.annotations <- column_to_rownames(merged_MN_data_spread[,which(colnames(merged_MN_data_spread) %in% c("sample_num","lesion", "Lesion_severity_score"))], var = "sample_num")
                    MN_colour = list(
                      lesion = c(`Negative Control` = "#000000", Nonlesional = "#9C9C9C", Irritant = "#FFB659", `Acetone Vehicle` = "#9ABFBA", `Day2 Allergy` = "#52B6FF", `Day4 Allergy` = "#4840FF"),
                      Lesion_severity_score = c(`-` = "#080808", `?/+` = "#430606", `+` = "#7E0404", `++` = "#B90202", `+++` = "#F50000")
                      
                    )
                    
                    MN.olink.heatmap <- pheatmap(MN.counts[reordered.proteins,],
                                                 annotation_col = MN.pheat.annotations,
                                                 annotation_colors = MN_colour,
                                                 cutree_cols = 1,
                                                 cutree_rows = 1,
                                                 color = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                 border_color=NA,
                                                 treeheight_row = 0,
                                                 treeheight_col = 0,
                                                 show_colnames = F,
                                                 cluster_rows = F)
        
              #Row-scaled Heatmap:
                    positive.MN.data <- positive.olink.data %>% filter(dilution_method != "Old method" & sample_type == "Microneedle")
                    sample_stats <- positive.MN.data %>% group_by(sample_num) %>% summarize(average_concentration = mean(Concentration),
                                                                                            max_concentration = max(Concentration))
                    
                    bad.samples <- sample_stats %>% filter(max_concentration <1000) %>% dplyr::pull(sample_num)
                    
                    Merged_Olink_data <- positive.olink.data %>% filter(sample_type == "Blister" | dilution_method != "Old method" & sample_type == "Microneedle" & !sample_num %in% bad.samples)

                    Merged_Olink_data_spread <- Merged_Olink_data %>% pivot_wider(names_from = Protein, values_from = Concentration)
                                  
                    #Convert back to log2 NPX values:
                    Merged.counts = log2(t(Merged_Olink_data_spread[,c(which(colnames(Merged_Olink_data_spread) == "IL8"):length(colnames(Merged_Olink_data_spread)))]))
                    colnames(Merged.counts) <- Merged_Olink_data_spread$sample_num
                    
                    meta.data <- Merged_Olink_data_spread[,c(1:which(colnames(Merged_Olink_data_spread) == "IL8")-1)]
                    
                    
                    Merged.pheat.annotations <- column_to_rownames(Merged_Olink_data_spread[,which(colnames(Merged_Olink_data_spread) %in% c("sample_num","lesion", "Lesion_severity_score", "sample_type"))], var = "sample_num")
                    Merged_colour = list(
                      lesion = c(`Negative Control` = "#000000", Nonlesional = "#9C9C9C", Irritant = "#FFB659", `Acetone Vehicle` = "#9ABFBA", `Day2 Allergy` = "#52B6FF", `Day4 Allergy` = "#4840FF"),
                      Lesion_severity_score = c(`-` = "#080808", `?/+` = "#430606", `+` = "#7E0404", `++` = "#B90202", `+++` = "#F50000"),
                      sample_type = c(`Blister` = "#07B1AF", `Microneedle` = "#E5F000")
                    )
                    
                    Merged.olink.heatmap <- pheatmap(Merged.counts[reordered.proteins,],
                                                 annotation_col = Merged.pheat.annotations,
                                                 annotation_colors = Merged_colour,
                                                 cutree_cols = 1,
                                                 cutree_rows = 1,
                                                 color = colorRampPalette(c("#005775","#003A4E","#000000","#792B45","#EC1F63"))(256),
                                                 border_color=NA,
                                                 treeheight_row = 0,
                                                 treeheight_col = 0,
                                                 show_colnames = F,
                                                 cluster_rows = T,
                                                 scale = "row")
                    
                    Merged.olink.heatmap
                    
                    #ggsave(Merged.olink.heatmap, width = 40, height = 20,units = "cm",filename = "Plots/olink/rowscaled_olink_heatmap.pdf", limitsize = FALSE)  
                    
                    
#PCA Analysis:
          #Blister:
          blister_data_for_pca <- filtered.olink.data %>% filter(sample_type == "Blister") %>% pivot_wider(names_from = Protein, values_from = Concentration)
          
          #PCA using prcomp:
          meta.data <- blister_data_for_pca[,c(1:which(colnames(blister_data_for_pca) == "IL8")-1)]
          set.seed(100)
          pca <- prcomp(log2(blister_data_for_pca[,c(which(colnames(blister_data_for_pca) == "IL8"):length(colnames(blister_data_for_pca)))]))
          
          blister_pca <- cbind(meta.data, pca$x)
         
          
          #Build Custom Color Palette:
          color.palette <- c("#080808","#F50000")
          GetPalette = colorRampPalette(color.palette)
          
          ColorCount = as.numeric(length(unique(blister_pca$Lesion_severity_score)))
          custom.colors <- GetPalette(ColorCount)
          
          ggplot(blister_pca, aes(x = PC1, y = PC2, color = Lesion_severity_score, shape = lesion)) +
            geom_point(size = 3) +
            theme_classic()+
            scale_color_manual(values = custom.colors) +
            scale_shape_manual(values = c(16, 8, 0, 7, 12))
          
          blister_pc_plot <- last_plot()
          
          
          Lesion.colors <- c("#9C9C9C", "#FFB659", "#9ABFBA", "#52B6FF", "#4840FF")
          
          ggplot(blister_pca, aes(x = PC1, y = PC2, color = lesion)) +
            geom_point(size = 3) +
            theme_classic() +
            scale_color_manual(values = Lesion.colors) 
          
          blister_pc_plot <- last_plot()
          
          #ggsave(blister_pc_plot, width = 15, height = 10, units = "cm",filename = "Plots/olink/olink.pc.plot.pdf", limitsize = FALSE)  
          
          
          #PC loadings:
          blister.loadings <- as.data.frame(pca$rotation[,1:2])
          blister.loadings <- rownames_to_column(blister.loadings, var = "Protein")
          
          DE.proteins <- olink.blister.DE %>% filter(ref_comp %in% c("Irritant vs Day2 Allergy"))
          
          blister.loadings <- left_join(blister.loadings, DE.proteins[,which(colnames(DE.proteins) %in% c("Protein", "Color_factor"))])
          
          blister.loadings$Color_factor <- as.character(blister.loadings$Color_factor)
          blister.loadings$Color_factor[which(blister.loadings$Color_factor == "Significantly Up")] <- "Day2 Allergy enriched" 
          blister.loadings$Color_factor[which(blister.loadings$Color_factor == "Significantly Down")] <- "Irritant enriched"
          blister.loadings$Color_factor[which(blister.loadings$Color_factor == "Not Significant")] <- "Not Significantly DE"
          
          blister.loadings$Color_factor <- factor(blister.loadings$Color_factor,
                                                  levels = c("Day2 Allergy enriched","Irritant enriched", "Not Significantly DE"))
          
          ggplot(blister.loadings, aes(x = PC1, y = PC2, label = Protein, color = Color_factor)) +
            geom_text_repel() +
            scale_color_manual(values = c("red", "blue", "black")) +
            theme_classic()
          
          blister_pc_loadings <- last_plot()
          
          
          ggplot(blister_data_for_pca, aes(y = log2(IL4), x = log2(IFNG), color = Lesion_severity_score)) +
            geom_point()+
            theme_classic() +
            scale_color_manual(values = custom.colors)
           
          blister_IFNG_IL4_dotplot <- last_plot()
          
          #ggsave(blister_IFNG_IL4_dotplot, width = 15, height = 10, units = "cm",filename = "Plots/olink/olink.IFNG.IL4.dotplot.pdf", limitsize = FALSE)  
          
          
          
          
          
          
#Protein Classifier Comparison:
          all.olink.proteins <- colnames(blister_data_for_pca)[c(which(colnames(blister_data_for_pca) == "IL8"):length(colnames(blister_data_for_pca)))]      
          
      ## Regularized Regression Classifier:
          blister_data_for_pca <- blister_data_for_pca %>% mutate(lesion_classified = case_when(lesion %in% c("Day4 Allergy", "Day2 Allergy") ~ 1,
                                                                                               TRUE ~ 0))


          y <- blister_data_for_pca %>% select(lesion_classified) %>% as.matrix()
          X <- blister_data_for_pca[,which(colnames(blister_data_for_pca) %in% all.olink.proteins)] %>% as.matrix()
          
          lambdas_to_try <- 10^seq(-3, 3, length.out = 100)

            #AOC Measure:
              set.seed(100)
              lasso_cv <- cv.glmnet(X, y, alpha = 1, lambda = lambdas_to_try, family = "binomial", type.measure = "auc",
                                    standardize = TRUE, nfolds = 10)
    
              # Plot cross-validation results
              plot(lasso_cv)
            
              lambda.1se <- lasso_cv$lambda.1se
              regularized_cv_lasso_coefficients <- as.data.frame(coef(lasso_cv)) %>% filter(s1 != 0)
              colnames(regularized_cv_lasso_coefficients) <- "lambda.1se"
              regularized_cv_lasso_coefficients <- rownames_to_column(regularized_cv_lasso_coefficients, var = "Parameter")
              
              lambda.min <- lasso_cv$lambda.min
              min_cv_error_lasso_coefficients <- as.data.frame(coef(lasso_cv, s = "lambda.min")) %>% filter(s1 != 0)
              colnames(min_cv_error_lasso_coefficients) <- "lambda.min"
              min_cv_error_lasso_coefficients <- rownames_to_column(min_cv_error_lasso_coefficients, var = "Parameter")
              
                    #Find Lambdas that give lower degrees of freedom for comparison:
                    lasso_cv$nzero
                    df4.lambda <- lasso_cv$lambda[65]
                    df3.lambda <- lasso_cv$lambda[62]
                    df1.lambda <- lasso_cv$lambda[60]
              
                    lasso_df4 <- glmnet(X,y, alpha = 1, lambda = df4.lambda, family = "binomial", type.measure = "auc")
                    lasso_df3 <- glmnet(X,y, alpha = 1, lambda = df3.lambda, family = "binomial", type.measure = "auc")
                    lasso_df1 <- glmnet(X,y, alpha = 1, lambda = df1.lambda, family = "binomial", type.measure = "auc")
                    
              
                    df4.coefs <- as.data.frame(coef(lasso_df4)) %>% filter(s0 != 0)
                    df3.coefs <- as.data.frame(coef(lasso_df3)) %>% filter(s0 != 0)
                    df1.coefs <- as.data.frame(coef(lasso_df1)) %>% filter(s0 != 0)
                    
                    colnames(df4.coefs) <- "df4.coefs"
                    df4.coefs <- rownames_to_column(df4.coefs, var = "Parameter")
                    
                    colnames(df3.coefs) <- "df3.coefs"
                    df3.coefs <- rownames_to_column(df3.coefs, var = "Parameter")
                    
                    colnames(df1.coefs) <- "df1.coefs"
                    df1.coefs <- rownames_to_column(df1.coefs, var = "Parameter")
                    
                    
                    coef.df <- purrr::reduce(list(min_cv_error_lasso_coefficients,regularized_cv_lasso_coefficients, df4.coefs, df3.coefs, df1.coefs), full_join)
                    coef.df[is.na(coef.df)]<-0
          
                              
          ## Dot Plots for each model:
                   
                    for (c in 1:c(ncol(coef.df)-1)) {
                      
                      model <- colnames(coef.df)[1+c]
                      
                      model.coefs <- coef.df[,c(1,1+c)]
                      
                              if (c == 1) {
                                
                              predicted_data <- blister_data_for_pca %>% mutate(model = model.coefs[1,2] + eval(as.symbol(model.coefs$Parameter[2]))*model.coefs[2,2] +
                                                                                                           eval(as.symbol(model.coefs$Parameter[3]))*model.coefs[3,2] +
                                                                                                           eval(as.symbol(model.coefs$Parameter[4]))*model.coefs[4,2] +
                                                                                                           eval(as.symbol(model.coefs$Parameter[5]))*model.coefs[5,2] +
                                                                                                           eval(as.symbol(model.coefs$Parameter[6]))*model.coefs[6,2] +
                                                                                                           eval(as.symbol(model.coefs$Parameter[7]))*model.coefs[7,2] +
                                                                                                           eval(as.symbol(model.coefs$Parameter[8]))*model.coefs[8,2])
                            
                              colnames(predicted_data)[which(colnames(predicted_data) == "model")] <- paste(model, "_model", sep = "")
                              
                              }
                      
                              if (c > 1) {
                                
                                predicted_data <- predicted_data %>% mutate(model = model.coefs[1,2] + eval(as.symbol(model.coefs$Parameter[2]))*model.coefs[2,2] +
                                                                                    eval(as.symbol(model.coefs$Parameter[3]))*model.coefs[3,2] +
                                                                                    eval(as.symbol(model.coefs$Parameter[4]))*model.coefs[4,2] +
                                                                                    eval(as.symbol(model.coefs$Parameter[5]))*model.coefs[5,2] +
                                                                                    eval(as.symbol(model.coefs$Parameter[6]))*model.coefs[6,2] +
                                                                                    eval(as.symbol(model.coefs$Parameter[7]))*model.coefs[7,2] +
                                                                                    eval(as.symbol(model.coefs$Parameter[8]))*model.coefs[8,2])
                                
                                colnames(predicted_data)[which(colnames(predicted_data) == "model")] <- paste(model, "_model", sep = "")
                                
                              }
                      
                      
                    } 
                    
                    
                    predicted_data$lambda.min_model
                    
                    p1 <- ggplot(predicted_data, aes(x = lesion, y = lambda.min_model, color = Lesion_severity_score)) +
                      geom_jitter(width = 0.2, height = 0)+
                      theme_classic() +
                      scale_color_manual(values = custom.colors) +
                      labs(y = "Classifier Value", title = "lambda.min_model") +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                    p2 <- ggplot(predicted_data, aes(x = lesion, y = lambda.1se_model, color = Lesion_severity_score)) +
                      geom_jitter(width = 0.2, height = 0)+
                      theme_classic() +
                      scale_color_manual(values = custom.colors) +
                      labs(y = "Classifier Value", title = "lambda.1se_model") +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                    p3 <- ggplot(predicted_data, aes(x = lesion, y = df4.coefs_model, color = Lesion_severity_score)) +
                      geom_jitter(width = 0.2, height = 0)+
                      theme_classic() +
                      scale_color_manual(values = custom.colors) +
                      labs(y = "Classifier Value", title = "4-protein model") +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                    p4 <- ggplot(predicted_data, aes(x = lesion, y = df3.coefs_model, color = Lesion_severity_score)) +
                      geom_jitter(width = 0.2, height = 0)+
                      theme_classic() +
                      scale_color_manual(values = custom.colors) +
                      labs(y = "Classifier Value", title = "3-protein model") +
                      theme(plot.title = element_text(hjust = 0.5))
                    
                   p5 <- ggplot(predicted_data, aes(x = lesion, y = df1.coefs_model, color = Lesion_severity_score)) +
                      geom_jitter(width = 0.2, height = 0)+
                      theme_classic() +
                      scale_color_manual(values = custom.colors) +
                      labs(y = "Classifier Value", title = "1-protein model") +
                      theme(plot.title = element_text(hjust = 0.5))
                   
                   
                plot.list <- list(p1, p2, p3, p4, p5)
                
                LASSO.plots <- ggpubr::ggarrange(plotlist = plot.list, common.legend = TRUE, legend = "bottom", nrow = 1)
                LASSO.plots
                 
          ## AOC Curves for each model:
                   
                model.colnames <- colnames(predicted_data)[which(str_detect(colnames(predicted_data), pattern = "model"))]
                
                roc.input.data <- predicted_data[,which(colnames(predicted_data) %in% c("lesion",model.colnames))]
                
                n_steps = 500
                
                for (c in 1:c(ncol(roc.input.data)-1)) {
                  
                  model <- colnames(roc.input.data)[1+c]
                  
                  input.data <- roc.input.data[,c(1,1+c)]
                  
                  range = range(input.data[,2])[2] - range(input.data[,2])[1]
                  
                  step.size = range/n_steps
                  
                  roc.data <- data.frame(threshold = c(range(input.data[,2])[1]-step.size,seq(range(input.data[,2])[1], range(input.data[,2])[2], by = step.size)),
                                         sensitivity = NA,
                                         specificity = NA,
                                         TPR = NA,
                                         FPR = NA)
                  
                          for (r in 1:nrow(roc.data)) {
                          
                            thres <- roc.data$threshold[r]
                            
                            prediction.df <- input.data %>% mutate(predicted_lesion = case_when(eval(as.symbol(colnames(input.data)[2])) > thres ~ "allergy",
                                                                                                eval(as.symbol(colnames(input.data)[2])) <= thres ~ "non_allergy"))
                            
                            prediction.df <- prediction.df %>% mutate(lesion_category = case_when(lesion %in% c("Day2 Allergy", "Day4 Allergy") ~ "allergy",
                                                                                                  TRUE ~ "non_allergy"))
                            
                            table <- table(prediction.df$predicted_lesion, prediction.df$lesion_category)
                            
                            if (nrow(table) ==2) {
                              sens = table[1,1]/sum(table[1,1], table[2,1])
                              spec = table[2,2]/sum(table[1,2], table[2,2])  
                              TPR = table[1,1]/sum(table[1,1], table[2,1])
                              FPR = table[1,2]/sum(table[1,2], table[2,2]) 
                            }
                            
                            if (nrow(table) == 1 & rownames(table)[1] == "non_allergy") {
                              sens = 0
                              spec = 1
                              TPR = 0
                              FPR = 0
                            }
                            if (nrow(table) == 1 & rownames(table)[1] == "allergy") {
                              sens = 1
                              spec = 0
                              TPR = 1
                              FPR = 1
                            }
                            
                            roc.data$sensitivity[r] <- sens
                            roc.data$specificity[r] <- spec
                            roc.data$TPR[r] <- TPR
                            roc.data$FPR[r] <- FPR
                            
                            rm(sens)
                            rm(spec)
                            rm(TPR)
                            rm(FPR)
                          }
                  
                  roc.data$model <- model
                  
                  if (c == 1) {
                    roc.model.data <- roc.data
                  }
                
                  if (c > 1) {
                    roc.model.data <- bind_rows(roc.model.data, roc.data)
                  }
                    
                } 
                
                roc.model.data$model <- factor(roc.model.data$model,
                                               levels = c("lambda.min_model", "lambda.1se_model", "df4.coefs_model", "df3.coefs_model", "df1.coefs_model"))
                
                ggplot(roc.model.data, aes(x = FPR, y = TPR, color = model)) +
                  geom_path(size = 1.5) +
                  theme_classic() +
                  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
                  labs(x = "False Positive Rate", y = "True Positive Rate") +
                  scale_color_manual(values = c("#FF0000", "#FFA200", "#CFD200", "#359F03", "#03709F")) +
                  labs(title = "LASSO Model ROC Curves")+
                  theme(plot.title = element_text(hjust=0.5))+
                  geom_point(x = 0.06451613, y =0.9310345, size = 4, color = "black")
    
              ROC.curve.plot <- last_plot()
              
              best.thresholds <- roc.model.data %>% filter(model == "lambda.min_model" & TPR > 0.75 & FPR <0.1)
              best.threshold <- roc.model.data %>% filter(model == "lambda.min_model" & sensitivity > 0.93 & specificity > 0.93)
                    
              ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`IL-1 alpha`), color = Lesion_severity_score)) +
                geom_jitter(width = 0.2, height = 0)+
                theme_classic() +
                scale_color_manual(values = custom.colors) +
                labs(title = "IL1a") +
                theme(plot.title = element_text(hjust =0.5))

              IL1a.protein.plot <- last_plot()
              
              ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`TSLP`), color = Lesion_severity_score)) +
                geom_jitter(width = 0.2, height = 0)+
                theme_classic() +
                scale_color_manual(values = custom.colors) +
                labs(title = "TSLP") +
                theme(plot.title = element_text(hjust =0.5))
              
              TSLP.protein.plot <- last_plot()
              
              
              ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`CCL19`), color = Lesion_severity_score)) +
                geom_jitter(width = 0.2, height = 0)+
                theme_classic() +
                scale_color_manual(values = custom.colors) +
                labs(title = "CCL19") +
                theme(plot.title = element_text(hjust =0.5))
              
              CCL19.protein.plot <- last_plot()
          
              
              ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`CXCL10`), color = Lesion_severity_score)) +
                geom_jitter(width = 0.2, height = 0)+
                theme_classic() +
                scale_color_manual(values = custom.colors) +
                labs(title = "CXCL10") +
                theme(plot.title = element_text(hjust =0.5))
              
              CXCL10.protein.plot <- last_plot()
              
              ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`IL4`), color = Lesion_severity_score)) +
                geom_jitter(width = 0.2, height = 0)+
                theme_classic() +
                scale_color_manual(values = custom.colors) +
                labs(title = "IL4") +
                theme(plot.title = element_text(hjust =0.5))
              
              IL4.protein.plot <- last_plot()
              
              ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`CCL20`), color = Lesion_severity_score)) +
                geom_jitter(width = 0.2, height = 0)+
                theme_classic() +
                scale_color_manual(values = custom.colors) +
                labs(title = "CCL20") +
                theme(plot.title = element_text(hjust =0.5))
              
              CCL20.protein.plot <- last_plot()
              
              ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`OPG`), color = Lesion_severity_score)) +
                geom_jitter(width = 0.2, height = 0)+
                theme_classic() +
                scale_color_manual(values = custom.colors) +
                labs(title = "OPG") +
                theme(plot.title = element_text(hjust =0.5))
              
              OPG.protein.plot <- last_plot()
              
      plot.list = list(CXCL10.protein.plot,IL4.protein.plot, OPG.protein.plot, CCL19.protein.plot,CCL20.protein.plot,TSLP.protein.plot,IL1a.protein.plot)
      
      Protein.plots <- ggpubr::ggarrange(plotlist = plot.list, common.legend = TRUE, legend = "bottom", nrow = 1)
      Protein.plots
       
      
      #ggsave(Protein.plots, width = 50, height = 15, units = "cm",filename = "Plots/olink/protein.biomarkers.pdf", limitsize = FALSE)  
      #ggsave(LASSO.plots, width = 80, height = 15, units = "cm",filename = "Plots/olink/lasso.model.predictions.pdf", limitsize = FALSE)  
      #ggsave(ROC.curve.plot, width = 25, height = 20, units = "cm",filename = "Plots/olink/lasso.model.ROC.curves.pdf", limitsize = FALSE)  
      
      #writexl::write_xlsx(coef.df, path = "Plots/olink/lasso_model_coefficients.xlsx")
      
      
      
      ### Multiplicative Classifiers:
          all.olink.proteins <- colnames(blister_data_for_pca)[c(which(colnames(blister_data_for_pca) == "IL8"):length(colnames(blister_data_for_pca)))]      
          
          
          #Single proteins as reliable classifiers:
              for (p in 1:length(all.olink.proteins)) {
                
                prot <- all.olink.proteins[p]
              
                prot.data <- cbind(meta.data[,which(colnames(meta.data) %in% c("lesion", "Lesion_severity_score"))], blister_data_for_pca[,which(colnames(blister_data_for_pca) == prot)])
              
                #100 percent specificity threshold:
                    max.lesion.vals <- prot.data %>% group_by(lesion) %>% summarize(max(eval(as.symbol(prot)))) %>% as.data.frame() %>% 
                                                    filter(!lesion %in% c("Day2 Allergy", "Day4 Allergy"))
                    
                    colnames(max.lesion.vals)[2] <- "max_protein_val"
                    
                    max.control.val <- max(max.lesion.vals$max_protein_val)
                
                #Day 2 and Day 4 Allergy Sensitivity:
                
                   Day2sample.count <- prot.data %>% filter(lesion == "Day2 Allergy") %>% summarize(sample.count = n()) %>% 
                                            as.data.frame() %>% dplyr::pull(sample.count)
                   
                   Day4sample.count <- prot.data %>% filter(lesion == "Day4 Allergy") %>% summarize(sample.count = n()) %>% 
                     as.data.frame() %>% dplyr::pull(sample.count)
                   
                   
                   Day2sensitivity <- prot.data %>% filter(lesion == "Day2 Allergy" & eval(as.symbol(prot)) > max.control.val) %>% summarize(sample.count = n()) %>% 
                     as.data.frame() %>% dplyr::pull(sample.count)
                   
                   Day2sensitivity <- round(Day2sensitivity/Day2sample.count, digits = 3)
                   
                   Day4sensitivity <- prot.data %>% filter(lesion == "Day4 Allergy" & eval(as.symbol(prot)) > max.control.val) %>% summarize(sample.count = n()) %>% 
                                                as.data.frame() %>% dplyr::pull(sample.count)
                   
                   Day4sensitivity <- round(Day4sensitivity/Day4sample.count, digits = 3)
                  
               #Output:
                   if (p ==1) {
                     single.prot.classifiers <- data.frame(Protein = prot,
                                                           Day2sensitivity = Day2sensitivity,
                                                           Day4sensitivity = Day4sensitivity)
                   }
                   if (p > 1) {
                     prot.classifier <- data.frame(Protein = prot,
                                                           Day2sensitivity = Day2sensitivity,
                                                           Day4sensitivity = Day4sensitivity)
                     
                     single.prot.classifiers <- bind_rows(single.prot.classifiers, prot.classifier)
                   }
                   rm(prot)
                   rm(Day2sensitivity)
                   rm(Day4sensitivity)
                   
              }
          
              ggplot(single.prot.classifiers, aes(x = Day2sensitivity, y = Day4sensitivity, label = Protein))+
                theme_classic() +
                geom_text_repel() +
                scale_x_continuous(limits = c(0,1)) +
                scale_y_continuous(limits = c(0,1)) +
                labs(title = "Single Protein Allergy Sensitivity at 100% Specificity", x = "Day 2 Allergy Sensitivity", y = "Day 4 Allergy Sensitivity") +
                theme(plot.title = element_text(hjust = 0.5))
              
          #Dual proteins as reliable classifiers:
            Two.prot.combos <- as.data.frame(t(as.data.frame(combn(all.olink.proteins, 2), cols = 2)))

            ranked.single.prots <- single.prot.classifiers %>% arrange(desc(Day2sensitivity)) %>% dplyr::pull(Protein)
            best.single.prots <- ranked.single.prots[1:10]
            
            Two.prot.combos <- Two.prot.combos %>% filter(V1 %in% best.single.prots | V2 %in% best.single.prots) %>% t() %>% as.data.frame()
         
              for (r in 1:ncol(Two.prot.combos)) {
                
                prots <- Two.prot.combos[,r]
                
                prot.data <- cbind(meta.data[,which(colnames(meta.data) %in% c("lesion", "Lesion_severity_score"))], blister_data_for_pca[,which(colnames(blister_data_for_pca) %in% prots)])
                
                prot.data <- prot.data %>% mutate(multiplicative = eval(as.symbol(prots[1]))*eval(as.symbol(prots[2])),
                                                   divisional1 = eval(as.symbol(prots[1]))/eval(as.symbol(prots[2])),
                                                   divisional2 = eval(as.symbol(prots[2]))/eval(as.symbol(prots[1])))

                colnames(prot.data)[-c(1:4)] <- c(paste(prots[1],"*",prots[2], sep = ""),
                                                  paste(prots[1],"/",prots[2], sep = ""),                
                                                  paste(prots[2],"/",prots[1], sep = ""))
                
                
                
                #100 percent specificity threshold:
                max.lesion.vals <- prot.data %>% group_by(lesion) %>% summarise_at(vars(c(paste(prots[1],"*",prots[2], sep = ""),
                                                                                        paste(prots[1],"/",prots[2], sep = ""),                
                                                                                        paste(prots[2],"/",prots[1], sep = ""))), funs(max)) %>% 
                                    as.data.frame() %>%
                                    filter(!lesion %in% c("Day2 Allergy", "Day4 Allergy"))
                
                max.lesion.vals <- column_to_rownames(max.lesion.vals, var = "lesion")
                
                for (c in 1:ncol(max.lesion.vals)) {
              
                  combo = colnames(max.lesion.vals)[c]
                  
                  max.control.val <- max(max.lesion.vals[,c])
                  
                          #Day 2 and Day 4 Allergy Sensitivity:
                          Day2sample.count <- prot.data %>% filter(lesion == "Day2 Allergy") %>% summarize(sample.count = n()) %>% 
                            as.data.frame() %>% dplyr::pull(sample.count)
                          
                          Day4sample.count <- prot.data %>% filter(lesion == "Day4 Allergy") %>% summarize(sample.count = n()) %>% 
                            as.data.frame() %>% dplyr::pull(sample.count)
                          
                          
                          Day2sensitivity <- prot.data %>% filter(lesion == "Day2 Allergy" & eval(as.symbol(combo)) > max.control.val) %>% summarize(sample.count = n()) %>% 
                            as.data.frame() %>% dplyr::pull(sample.count)
                          
                          Day2sensitivity <- round(Day2sensitivity/Day2sample.count, digits = 3)
                          
                          Day4sensitivity <- prot.data %>% filter(lesion == "Day4 Allergy" & eval(as.symbol(combo)) > max.control.val) %>% summarize(sample.count = n()) %>% 
                            as.data.frame() %>% dplyr::pull(sample.count)
                          
                          Day4sensitivity <- round(Day4sensitivity/Day4sample.count, digits = 3)
                          
                          
                          #Output:
                          if (r ==1 & c == 1) {
                            dual.prot.classifiers <- data.frame(Protein = combo,
                                                                  Day2sensitivity = Day2sensitivity,
                                                                  Day4sensitivity = Day4sensitivity)
                          }
                          if (r > 1 | c >1) {
                            prot.classifier <- data.frame(Protein = combo,
                                                          Day2sensitivity = Day2sensitivity,
                                                          Day4sensitivity = Day4sensitivity)
                            
                            dual.prot.classifiers <- bind_rows(dual.prot.classifiers, prot.classifier)
                          }
                          rm(combo)
                          rm(Day2sensitivity)
                          rm(Day4sensitivity)
                          
                }
                
              }
          
         #Three-protein combinations as reliable classifiers:
             Three.prot.combos <- as.data.frame(t(as.data.frame(combn(all.olink.proteins, 3), cols = 3)))
             
             Three.prot.combos <- Three.prot.combos %>% filter(V1 %in% best.single.prots | V2 %in% best.single.prots | V3 %in% best.single.prots) %>% t() %>% as.data.frame()

             for (r in 1:ncol(Three.prot.combos)) {
               
               prots <- Three.prot.combos[,r]
               
               prot.data <- cbind(meta.data[,which(colnames(meta.data) %in% c("lesion", "Lesion_severity_score"))], blister_data_for_pca[,which(colnames(blister_data_for_pca) %in% prots)])
               
               prot.data <- prot.data %>% mutate(multiplicative = eval(as.symbol(prots[1]))*eval(as.symbol(prots[2]))*eval(as.symbol(prots[3])),
                                                 divisional1 = eval(as.symbol(prots[1]))*eval(as.symbol(prots[2]))/eval(as.symbol(prots[3])),
                                                 divisional2 = eval(as.symbol(prots[1]))*eval(as.symbol(prots[3]))/eval(as.symbol(prots[2])),
                                                 divisional3 = eval(as.symbol(prots[2]))*eval(as.symbol(prots[3]))/eval(as.symbol(prots[1])),
                                                 
                                                 divisional4 = eval(as.symbol(prots[1]))/eval(as.symbol(prots[2]))*eval(as.symbol(prots[3])),
                                                 divisional5 = eval(as.symbol(prots[2]))/eval(as.symbol(prots[1]))*eval(as.symbol(prots[3])),
                                                 divisional6 = eval(as.symbol(prots[3]))/eval(as.symbol(prots[1]))*eval(as.symbol(prots[2])))
                                                 
               
               colnames(prot.data)[-c(1:5)] <- c(paste(prots[1],"*",prots[2],"*", prots[3], sep = ""),
                                                 paste(prots[1],"*",prots[2],"/", prots[3], sep = ""),
                                                 paste(prots[1],"*",prots[3],"/", prots[2], sep = ""),
                                                 paste(prots[2],"*",prots[3],"/", prots[1], sep = ""),
                                                 
                                                 paste(prots[1],"/",prots[2],"*", prots[3], sep = ""),
                                                 paste(prots[2],"/",prots[1],"*", prots[3], sep = ""),
                                                 paste(prots[3],"/",prots[1],"*", prots[2], sep = ""))

               
               #100 percent specificity threshold:
               max.lesion.vals <- prot.data %>% group_by(lesion) %>% summarise_at(vars(c(paste(prots[1],"*",prots[2],"*", prots[3], sep = ""),
                                                                                         paste(prots[1],"*",prots[2],"/", prots[3], sep = ""),
                                                                                         paste(prots[1],"*",prots[3],"/", prots[2], sep = ""),
                                                                                         paste(prots[2],"*",prots[3],"/", prots[1], sep = ""),
                                                                                         
                                                                                         paste(prots[1],"/",prots[2],"*", prots[3], sep = ""),
                                                                                         paste(prots[2],"/",prots[1],"*", prots[3], sep = ""),
                                                                                         paste(prots[3],"/",prots[1],"*", prots[2], sep = ""))), funs(max)) %>% 
                 as.data.frame() %>%
                 filter(!lesion %in% c("Day2 Allergy", "Day4 Allergy"))
               
               max.lesion.vals <- column_to_rownames(max.lesion.vals, var = "lesion")
               
               for (c in 1:ncol(max.lesion.vals)) {
                 
                 combo = colnames(max.lesion.vals)[c]
                 
                 max.control.val <- max(max.lesion.vals[,c])
                 
                 #Day 2 and Day 4 Allergy Sensitivity:
                 Day2sample.count <- prot.data %>% filter(lesion == "Day2 Allergy") %>% summarize(sample.count = n()) %>% 
                   as.data.frame() %>% dplyr::pull(sample.count)
                 
                 Day4sample.count <- prot.data %>% filter(lesion == "Day4 Allergy") %>% summarize(sample.count = n()) %>% 
                   as.data.frame() %>% dplyr::pull(sample.count)
                 
                 
                 Day2sensitivity <- prot.data %>% filter(lesion == "Day2 Allergy" & eval(as.symbol(combo)) > max.control.val) %>% summarize(sample.count = n()) %>% 
                   as.data.frame() %>% dplyr::pull(sample.count)
                 
                 Day2sensitivity <- round(Day2sensitivity/Day2sample.count, digits = 3)
                 
                 Day4sensitivity <- prot.data %>% filter(lesion == "Day4 Allergy" & eval(as.symbol(combo)) > max.control.val) %>% summarize(sample.count = n()) %>% 
                   as.data.frame() %>% dplyr::pull(sample.count)
                 
                 Day4sensitivity <- round(Day4sensitivity/Day4sample.count, digits = 3)
                 
                 
                 #Output:
                 if (r ==1 & c == 1) {
                   triple.prot.classifiers <- data.frame(Protein = combo,
                                                       Day2sensitivity = Day2sensitivity,
                                                       Day4sensitivity = Day4sensitivity)
                 }
                 if (r > 1 | c >1) {
                   prot.classifier <- data.frame(Protein = combo,
                                                 Day2sensitivity = Day2sensitivity,
                                                 Day4sensitivity = Day4sensitivity)
                   
                   triple.prot.classifiers <- bind_rows(triple.prot.classifiers, prot.classifier)
                 }
                 rm(combo)
                 rm(Day2sensitivity)
                 rm(Day4sensitivity)
                 
               }
               
             }
             
             
             single.prot.classifiers$Protein_number <- "single"
             dual.prot.classifiers$Protein_number <- "dual"
             triple.prot.classifiers$Protein_number <- "triple"
             multi.prot.classifiers <- bind_rows(single.prot.classifiers, dual.prot.classifiers, triple.prot.classifiers)
             multi.prot.classifiers$Protein_number <- factor(multi.prot.classifiers$Protein_number,
                                                             levels = c("single", "dual", "triple"))
             
             multi.prot.classifiers <- multi.prot.classifiers %>% arrange(desc(Protein_number))
             
             #save(multi.prot.classifiers, file = "Data/Olink/olink.protein.classifier.sensitivity.Rdata")
             
             load("Data/Olink/olink.protein.classifier.sensitivity.Rdata")
             
             single.prot.classifiers <-multi.prot.classifiers %>% filter(Protein_number == "single")
             
             ggplot(single.prot.classifiers, aes(x = Day2sensitivity, y = Day4sensitivity, label = Protein))+
               theme_classic() +
               geom_text_repel(size = 3) +
               scale_x_continuous(limits = c(0,1)) +
               scale_y_continuous(limits = c(0,1)) +
               labs(title = "Single Protein Allergy Sensitivity at 100% Specificity", x = "Day 2 Allergy Sensitivity", y = "Day 4 Allergy Sensitivity") +
               theme(plot.title = element_text(hjust = 0.5))

             single.protein.biomarker.sensitivity <- last_plot()
             
             ggplot(multi.prot.classifiers, aes(x = Day2sensitivity, y = Day4sensitivity, color = Protein_number))+
               theme_minimal() +
               geom_jitter(width = 0.02, height = 0.02, alpha = 0.8, size = 2) +
               scale_x_continuous(limits = c(0,1)) +
               scale_y_continuous(limits = c(0,1)) +
               labs(title = "Single Protein Allergy Sensitivity at 100% Specificity", x = "Day 2 Allergy Sensitivity", y = "Day 4 Allergy Sensitivity") +
               theme(plot.title = element_text(hjust = 0.5)) +
               scale_color_manual(values = c("black", "#FFAE00", "#FF4040")) +
               theme(axis.line = element_line(color="black", size = 0.5))
             
             classifier.protein.num.comparison <- last_plot()
             
             ggplot(blister_data_for_pca, aes(x = lesion, y = log2(IFNG), color = Lesion_severity_score)) +
               geom_jitter(width = 0.2, height = 0)+
               theme_classic() +
               scale_color_manual(values = custom.colors) +
               labs(y = "Classifier Value")  +
               geom_hline(yintercept = blister_data_for_pca %>% filter(!lesion %in% c("Day2 Allergy", "Day4 Allergy")) %>% dplyr::pull(IFNG) %>% log2(.) %>% max(),
                          linetype = "dashed", color = "gray")
             
             IFNG_classifier <- last_plot()
             
             ggplot(blister_data_for_pca, aes(x = lesion, y = log2(IL4), color = Lesion_severity_score)) +
               geom_jitter(width = 0.2, height = 0)+
               theme_classic() +
               scale_color_manual(values = custom.colors) +
               labs(y = "Classifier Value")  +
               geom_hline(yintercept = blister_data_for_pca %>% filter(!lesion %in% c("Day2 Allergy", "Day4 Allergy")) %>% dplyr::pull(IL4) %>% log2(.) %>% max(),
                          linetype = "dashed", color = "gray")
             
             IL4_classifier <- last_plot()
             
             multi.prot.classifiers <- multi.prot.classifiers %>% mutate(Average_sensitivity = (Day2sensitivity+Day4sensitivity)/2)
             
             ggplot(blister_data_for_pca, aes(x = lesion, y = log2(CCL11*IFNG/VEGFA), color = Lesion_severity_score)) +
               geom_jitter(width = 0.2, height = 0)+
               theme_classic() +
               scale_color_manual(values = custom.colors) +
               geom_hline(yintercept = blister_data_for_pca %>% filter(!lesion %in% c("Day2 Allergy", "Day4 Allergy")) %>% mutate(classifier = log2(CCL11*IFNG/VEGFA)) %>% dplyr::pull(classifier) %>% max(),
                          linetype = "dashed", color = "gray") +
               ggtitle("log2(CCL11*IFNG/VEGFA)") +
               theme(plot.title = element_text(hjust = 0.5))
             
             best_3prot_classifier <- last_plot()
             
             ggplot(blister_data_for_pca, aes(x = lesion, y = log2(`IL-1 alpha`*IFNG/CXCL1), color = Lesion_severity_score)) +
               geom_jitter(width = 0.2, height = 0)+
               theme_classic() +
               scale_color_manual(values = custom.colors) +
               geom_hline(yintercept = blister_data_for_pca %>% filter(!lesion %in% c("Day2 Allergy", "Day4 Allergy")) %>% mutate(classifier = log2(`IL-1 alpha`*IFNG/CXCL1)) %>% dplyr::pull(classifier) %>% max(),
                          linetype = "dashed", color = "gray") +
               ggtitle("log2(`IL-1 alpha`*IFNG/CXCL1)") +
               theme(plot.title = element_text(hjust = 0.5))
             
             second_3prot_classifier <- last_plot()
             
             #ggsave(single.protein.biomarker.sensitivity, width = 15, height = 15,units = "cm",filename = "Plots/olink/single.protein.biomarkersensitivity.pdf", limitsize = FALSE)  
             #ggsave(IFNG_classifier, width = 15, height = 15,units = "cm",filename = "Plots/olink/IFNG_classifier.pdf", limitsize = FALSE)  
             #ggsave(IL4_classifier, width = 15, height = 15,units = "cm",filename = "Plots/olink/IL4_classifier.pdf", limitsize = FALSE)  
             #ggsave(classifier.protein.num.comparison, width = 20, height = 15,units = "cm",filename = "Plots/olink/biomarker_protein_num_comparison.png", limitsize = FALSE)  
             #ggsave(best_3prot_classifier, width = 15, height = 15,units = "cm",filename = "Plots/olink/best_3prot_classifier.pdf", limitsize = FALSE)  
             #ggsave(second_3prot_classifier, width = 15, height = 15,units = "cm",filename = "Plots/olink/second_3prot_classifier.pdf", limitsize = FALSE)  
             
             
        ggplot(blister_data_for_pca, aes(x = lesion, y = log2(IFNG/IL8*`MCP-1`), color = Lesion_severity_score)) +
          geom_jitter(width = 0.2, height = 0)+
          theme_classic() +
          scale_color_manual(values = custom.colors) +
          labs(y = "Classifier Value") 
       
        
          
          ggplot(blister_data_for_pca, aes(x = lesion, y = log2(((IFNG+CXCL11+CXCL10+CXCL9)*IL4*IL8/(CCL20))), color = Lesion_severity_score)) +
            geom_jitter(width = 0.2, height = 0)+
            theme_classic() +
            scale_color_manual(values = custom.colors) +
            labs(y = "Classifier Value") +
            geom_hline(yintercept = 14.25, linetype = "dashed", color = "gray")
          
          blister_classifier_score <- last_plot()

#################################   
# Comparison to Bulk RNA Analysis:
#################################
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
load(file = "Data/contact_derm_filtered3.Rdata")
    
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
    
# Bulk Data:
    DE.input <- calc_agg_bulk_srt_new(DE.input, aggregate_by = c("Lesion", "Patient", "Blister_date", "Sequence_instrument"))
    bulk_data <- DE.input@assays$RNA@misc$aggregate.bulk

    
    sample.lesion.info <- DE.input@meta.data %>% group_by_at(c("Lesion", "Patient", "Blister_date", "Lesion_severity_score", "Lesion", "Sequence_instrument")) %>%
                                summarise(n()) %>%
                                rename_at("n()", ~ "cellcount") %>% as.data.frame()
    
    bulk_data <- bulk_data[,-which(colnames(bulk_data) %in% c("numcells", "numgenes", "proportion"))] 
     
    bulk_data <- full_join(bulk_data, sample.lesion.info)
    bulk_data$Lesion_severity_score[which(is.na(bulk_data$Lesion_severity_score))] <- "-"
    
    bulk_data <- cbind(bulk_data[,-which(colnames(bulk_data) %in% rownames(DE.input))], bulk_data[,which(colnames(bulk_data) %in% rownames(DE.input))])
    bulk_data$Lesion <- factor(bulk_data$Lesion,
                               levels = levels(DE.input@meta.data$Lesion))
    
    metadata <- bulk_data[,-which(colnames(bulk_data) %in% rownames(DE.input))]
    bulk_counts <- bulk_data[,which(colnames(bulk_data) %in% rownames(DE.input))]
    
    #PCA using prcomp:
          set.seed(100)
          pca <- prcomp(bulk_counts)
          
          bulk_RNA_pca <- cbind(metadata, pca$x)
          
              p1 <- ggplot(bulk_RNA_pca, aes(x = PC1, y = PC2, color = Sequence_instrument, shape = Lesion)) +
                      geom_jitter(size = 3, width = (max(bulk_RNA_pca$PC1)-min(bulk_RNA_pca$PC1))*0.05, height = (max(bulk_RNA_pca$PC2)-min(bulk_RNA_pca$PC2))*0.05) +
                      theme_classic()+
                      scale_shape_manual(values = c(16, 8, 0, 7, 12))
              
              p2 <- ggplot(bulk_RNA_pca, aes(x = PC1, y = PC3, color = Sequence_instrument, shape = Lesion)) +
                geom_jitter(size = 3, width = (max(bulk_RNA_pca$PC1)-min(bulk_RNA_pca$PC1))*0.05, height = (max(bulk_RNA_pca$PC3)-min(bulk_RNA_pca$PC3))*0.05) +
                theme_classic()+
                scale_shape_manual(values = c(16, 8, 0, 7, 12))
              
              p1+p2
              rawbulk_RNA_pc_plot <- last_plot()    
      
          #UMAP:
          set.seed(100)
          bulk.umap <- umap(bulk_RNA_pca[,which(colnames(bulk_RNA_pca) %in% c(paste("PC", seq(1:20), sep = "")))])
          colnames(bulk.umap$layout) <- c("UMAP_x", "UMAP_y")
          
          bulk_pca_umap <- cbind(bulk_RNA_pca, bulk.umap$layout)
          
          ggplot(bulk_pca_umap, aes(x = UMAP_x, y = UMAP_y, color = Sequence_instrument, shape = Lesion)) +
            geom_point(size = 3) +
            theme_classic()+
            scale_shape_manual(values = c(16, 8, 0, 7, 12))
          
          rawbulk_RNA_umap <- last_plot()    
          
    #Batch corrected: 
    library(Harman)
            
              harman.object <- harman(datamatrix = t(bulk_counts), 
                                      expt = as.vector(metadata$Lesion), 
                                      batch = as.vector(metadata$Sequence_instrument), 
                                      limit = 0.95)
              
              harman.corrected.df <- bulk_data
              harman.corrected.df[,which(colnames(harman.corrected.df) %in% rownames(DE.input))] <- t(reconstructData(harman.object))
              
              #PCA using prcomp:
              set.seed(100)
              pca.harman.corrected <- prcomp(harman.corrected.df[,which(colnames(harman.corrected.df) %in% rownames(DE.input))])
              
              bulk_RNA_pca.harman.corrected <- cbind(metadata, pca.harman.corrected$x)
              
              p1 <- ggplot(bulk_RNA_pca.harman.corrected, aes(x = PC1, y = PC2, color = Sequence_instrument, shape = Lesion)) +
                        geom_jitter(size = 3, width = (max(bulk_RNA_pca.harman.corrected$PC1)-min(bulk_RNA_pca.harman.corrected$PC1))*0.05, height = (max(bulk_RNA_pca.harman.corrected$PC2)-min(bulk_RNA_pca.harman.corrected$PC2))*0.05) +
                        theme_classic()+
                        scale_shape_manual(values = c(16, 8, 0, 7, 12))
              
              p2 <- ggplot(bulk_RNA_pca.harman.corrected, aes(x = PC1, y = PC3, color = Sequence_instrument, shape = Lesion)) +
                        geom_jitter(size = 3, width = (max(bulk_RNA_pca.harman.corrected$PC1)-min(bulk_RNA_pca.harman.corrected$PC1))*0.05, height = (max(bulk_RNA_pca.harman.corrected$PC3)-min(bulk_RNA_pca.harman.corrected$PC3))*0.05) +
                        theme_classic()+
                        scale_shape_manual(values = c(16, 8, 0, 7, 12))
              p1+p2
              
              bulk_RNA_pc_plot_harman_corrected <- last_plot()    
              
              #UMAP:
              set.seed(100)
              bulk.umap.harman.corrected <- umap(bulk_RNA_pca.harman.corrected[,which(colnames(bulk_RNA_pca.harman.corrected) %in% c(paste("PC", seq(1:20), sep = "")))])
              colnames(bulk.umap.harman.corrected$layout) <- c("UMAP_x", "UMAP_y")
              
              bulk_pca_umap_harman_corrected <- cbind(bulk_RNA_pca.harman.corrected, bulk.umap.harman.corrected$layout)
              
              ggplot(bulk_pca_umap_harman_corrected, aes(x = UMAP_x, y = UMAP_y, color = Sequence_instrument, shape = Lesion)) +
                geom_point(size = 3) +
                theme_classic()+
                scale_shape_manual(values = c(16, 8, 0, 7, 12))
              
              bulk_RNA_umap_harman_corrected <- last_plot()    
              
              #Build Custom Color Palette:
              color.palette <- c("#080808","#F50000")
              GetPalette = colorRampPalette(color.palette)
              
              ColorCount = as.numeric(length(unique(bulk_pca_umap_harman_corrected$Lesion_severity_score)))
              custom.colors <- GetPalette(ColorCount)
              
              
             p1 <- ggplot(bulk_RNA_pca.harman.corrected, aes(x = PC1, y = PC2, color = Lesion_severity_score, shape = Lesion)) +
                      geom_jitter(size = 3, width = (max(bulk_RNA_pca.harman.corrected$PC1)-min(bulk_RNA_pca.harman.corrected$PC1))*0.05, height = (max(bulk_RNA_pca.harman.corrected$PC2)-min(bulk_RNA_pca.harman.corrected$PC2))*0.05) +
                      theme_classic()+
                      scale_shape_manual(values = c(16, 8, 0, 7, 12))+
                      scale_color_manual(values = custom.colors)
              
         
              p2 <- ggplot(bulk_pca_umap_harman_corrected, aes(x = UMAP_x, y = UMAP_y, color = Lesion_severity_score, shape = Lesion)) +
                        geom_point(size = 3) +
                        theme_classic()+
                        scale_shape_manual(values = c(16, 8, 0, 7, 12)) +
                        scale_color_manual(values = custom.colors)
                
              p1+p2
              PC_UMAP_harman_corrected_severity_colored <- last_plot()    
              
    # Bulk RNA PCA analysis using only Olink Genes:
    
            load("Data/Olink/filtered.olink.data.Rdata")    
            #filter out attempted concentrated samples since it didn't even do anything:
            filtered.olink.data <- filtered.olink.data[-which(str_detect(filtered.olink.data$sample_num, pattern = "mf-")),]
            
            all.olink.proteins <- unique(filtered.olink.data$Protein)
            
            #Reformating dataframe:
            protein.genename.correction.df <- data.frame(original_name = all.olink.proteins)
            
            protein.genenames <- str_replace(string = all.olink.proteins, pattern = " ", replacement = "")
            protein.genenames[which(str_detect(protein.genenames, pattern = "LAP"))] <- "TGFB1"
            protein.genenames <- str_replace(protein.genenames, pattern = "alpha", replacement = "A")
            protein.genenames[which(protein.genenames == "EN-RAGE")] <- "S100A12"
            protein.genenames[which(protein.genenames == "TGF-A")] <- "TGFA"
            protein.genenames[which(protein.genenames == "Flt3L")] <- "FLT3LG"
            
            
            protein.genename.correction.df$initial.corrected <- protein.genenames
            
            #aliases2entrez conversion to entrez ID:
            file <- system.file("extdata", "HGNC.txt", package = "aliases2entrez") 
            HGNC <- update_symbols()
            symbols <- protein.genename.correction.df$initial.corrected
            ids <- convert_symbols(symbols, HGNC, c=1)
            
            protein.genename.correction.df <- cbind(protein.genename.correction.df, ids)
            names(protein.genename.correction.df)[names(protein.genename.correction.df) == "Symbols"] <- "A2E_Symbols"
            
            #clusterprofiler conversion from entrez ID to RNAseq-matching symbols:
            positive_ids <- ids[which(!is.na(ids$entrezID)),]
            id.alias <- bitr(positive_ids$entrezID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
            id.alias$ENTREZID <- as.numeric(id.alias$ENTREZID)
            
            protein.genename.correction.df <- left_join(protein.genename.correction.df, id.alias, by = c("entrezID" = "ENTREZID"))
            names(protein.genename.correction.df)[names(protein.genename.correction.df) == "SYMBOL"] <- "reformated_name"
            
            olink.genes <- protein.genename.correction.df$reformated_name

            #PCA using prcomp:
            set.seed(100)
            olinkgenes_bulkRNA_pca <- prcomp(harman.corrected.df[,which(colnames(harman.corrected.df) %in% olink.genes)])
            
            olinkgenes_bulkRNA_pca <- cbind(metadata, olinkgenes_bulkRNA_pca$x)
            
            #UMAP:
            set.seed(100)
            olinkgenes_bulkRNA_umap <- umap(olinkgenes_bulkRNA_pca[,which(colnames(olinkgenes_bulkRNA_pca) %in% c(paste("PC", seq(1:20), sep = "")))])
            colnames(olinkgenes_bulkRNA_umap$layout) <- c("UMAP_x", "UMAP_y")
            
            bulk_pca_umap_olinkgenes <- cbind(olinkgenes_bulkRNA_pca, olinkgenes_bulkRNA_umap$layout)
            
            #Build Custom Color Palette:
            color.palette <- c("#080808","#F50000")
            GetPalette = colorRampPalette(color.palette)
            
            ColorCount = as.numeric(length(unique(bulk_pca_umap_olinkgenes$Lesion_severity_score)))
            custom.colors <- GetPalette(ColorCount)
            
            
            p1 <- ggplot(bulk_pca_umap_olinkgenes, aes(x = PC1, y = PC2, color = Lesion_severity_score, shape = Lesion)) +
              geom_jitter(size = 3, width = (max(bulk_pca_umap_olinkgenes$PC1)-min(bulk_pca_umap_olinkgenes$PC1))*0.05, height = (max(bulk_pca_umap_olinkgenes$PC2)-min(bulk_pca_umap_olinkgenes$PC2))*0.05) +
              theme_classic()+
              scale_shape_manual(values = c(16, 8, 0, 7, 12))+
              scale_color_manual(values = custom.colors)
            
            
            p2 <- ggplot(bulk_pca_umap_olinkgenes, aes(x = UMAP_x, y = UMAP_y, color = Lesion_severity_score, shape = Lesion)) +
              geom_point(size = 3) +
              theme_classic()+
              scale_shape_manual(values = c(16, 8, 0, 7, 12)) +
              scale_color_manual(values = custom.colors)
            
            p1+p2
            
            bulk_pca_umap_olinkgenes_plot <- last_plot()    
            

        #Bulk Dot Plots:
            ggplot(bulk_data, aes(x = Lesion, y = log2(IFNG), color = Lesion_severity_score)) +
              geom_jitter(width = 0.2, height = 0, size = 2)+
              theme_classic() +
              scale_color_manual(values = custom.colors) +
              labs(y = "Classifier Value") 
            
            ggplot(bulk_data, aes(x = Lesion, y = log2(CCL20), color = Lesion_severity_score)) +
              geom_jitter(width = 0.2, height = 0, size = 2)+
              theme_classic() +
              scale_color_manual(values = custom.colors) +
              labs(y = "Classifier Value") 
           
            ggplot(bulk_data, aes(x = Lesion, y = log2(IL4), color = Lesion_severity_score)) +
              geom_jitter(width = 0.2, height = 0, size = 2)+
              theme_classic() +
              scale_color_manual(values = custom.colors) +
              labs(y = "Classifier Value") 
            
            ggplot(bulk_data, aes(x = Lesion, y = log2(CXCL8), color = Lesion_severity_score)) +
              geom_jitter(width = 0.2, height = 0, size = 2)+
              theme_classic() +
              scale_color_manual(values = custom.colors) +
              labs(y = "Classifier Value") 
            
  # Olink Gene bulk RNA vs Protein value correlation:
            #Reformat to Merge:
            
            olinkgenes_bulkRNA_raw <- cbind(metadata, bulk_data[,which(colnames(bulk_data) %in% olink.genes)])
            olinkgenes_bulkRNA_corrected <- cbind(metadata, harman.corrected.df[,which(colnames(harman.corrected.df) %in% olink.genes)])
            
            blister.protein.data <- cbind(meta.data, blister_data_for_pca[,c(which(colnames(blister_data_for_pca) == "IL8"):length(colnames(blister_data_for_pca)))])         
            colnames(blister.protein.data)[-which(colnames(blister.protein.data) %in% colnames(meta.data))] <- protein.genename.correction.df$reformated_name
            blister.protein.data$blister_date <- as.character(blister.protein.data$blister_date)
            blister.protein.data$blister_date <- as_date(blister.protein.data$blister_date)
            
            raw_RNA_vs_protein <- olinkgenes_bulkRNA_raw %>% pivot_longer(cols = which(colnames(olinkgenes_bulkRNA_raw) %in% olink.genes), names_to = "Gene", values_to = "RNA")
            corrected_RNA_vs_protein <- olinkgenes_bulkRNA_corrected %>% pivot_longer(cols = which(colnames(olinkgenes_bulkRNA_corrected) %in% olink.genes), names_to = "Gene", values_to = "RNA")
           
            blister.protein.data <- blister.protein.data %>% pivot_longer(cols = which(colnames(blister.protein.data) %in% olink.genes), names_to = "Gene", values_to = "Protein")

            raw_RNA_vs_protein <- raw_RNA_vs_protein[,which(colnames(raw_RNA_vs_protein) %in% c("Lesion", "Patient", "Blister_date", "Lesion_severity_score", "Gene", "RNA"))]
            corrected_RNA_vs_protein <- corrected_RNA_vs_protein[,which(colnames(corrected_RNA_vs_protein) %in% c("Lesion", "Patient", "Blister_date", "Lesion_severity_score", "Gene", "RNA"))]
            
            blister.protein.data <- blister.protein.data[,which(colnames(blister.protein.data) %in% c("lesion", "patient", "blister_date", "Lesion_severity_score", "Gene", "Protein"))]
            colnames(blister.protein.data) <- c("Blister_date", "Patient", "Lesion", "Lesion_severity_score", "Gene", "Protein") 
            blister.protein.data$Lesion <- as.character(blister.protein.data$Lesion)
            blister.protein.data$Lesion <- str_replace(blister.protein.data$Lesion, pattern = " ", replacement = "_")
  
            blister.protein.data$Lesion <- factor(blister.protein.data$Lesion,
                                                  levels = levels(metadata$Lesion))
            
            #Merge:
            raw_RNA_vs_protein <- full_join(raw_RNA_vs_protein, blister.protein.data)
            corrected_RNA_vs_protein <- full_join(corrected_RNA_vs_protein, blister.protein.data)
            
            raw_RNA_vs_protein <- raw_RNA_vs_protein %>% mutate(log2_RNA = log2(RNA),
                                                                log2_Protein = log2(Protein))
            
            raw_RNA_vs_protein$log2_RNA[which(is.infinite(raw_RNA_vs_protein$log2_RNA))] <- -5
            raw_RNA_vs_protein$log2_RNA[is.na(raw_RNA_vs_protein$log2_RNA)] <- -5

            corrected_RNA_vs_protein <- corrected_RNA_vs_protein %>% mutate(log2_RNA = log2(RNA),
                                                                            log2_Protein = log2(Protein))
            
            corrected_RNA_vs_protein$log2_RNA[which(is.infinite(corrected_RNA_vs_protein$log2_RNA))] <- -10
            corrected_RNA_vs_protein$log2_RNA[is.na(corrected_RNA_vs_protein$log2_RNA)] <- -10
            

            RNA_stats <- raw_RNA_vs_protein %>% group_by(Gene) %>% summarize(max_RNA = max(log2_RNA),
                                                                            average_RNA = mean(log2_RNA))
            
           genes.to.plot <- RNA_stats %>% filter(max_RNA > 1) %>% dplyr::pull(Gene)
           
           data.for.lms <- raw_RNA_vs_protein %>% filter(log2_RNA >= -1 & Gene %in% genes.to.plot)
           
           for (g in 1:length(unique(data.for.lms$Gene))) {
             gene = unique(data.for.lms$Gene)[g]
             
             gene.data <- data.for.lms %>% filter(Gene == gene)
           
             model <- lm(log2_Protein ~ log2_RNA, data = gene.data)  # Fit linear model
             
             coef <- as.data.frame(summary(model)$coefficients)
             intercept <- coef$Estimate[1]
             coef <- coef$Estimate[2]
             
             r.squared <- summary(model)$r.squared
             
             if (g ==1) {
               rna_protein_lm_stats = data.frame(Gene = gene,
                                                 Intercept = intercept,
                                                 Coefficient = coef,
                                                 Rsquared = r.squared)
             }
             if (g > 1) {
               rna_protein_gene.stat = data.frame(Gene = gene,
                                                  Intercept = intercept,
                                                  Coefficient = coef,
                                                  Rsquared = r.squared)
               
               rna_protein_lm_stats <- bind_rows(rna_protein_lm_stats,rna_protein_gene.stat)
             }
             rm(coef)
             rm(r.squared)
          }
           
           rna_protein_lm_stats
           
           
           library(ggpmisc)
           
          correlated.proteins <- rna_protein_lm_stats %>% filter(Coefficient >= 0.2 & Rsquared >= 0.2) %>% dplyr::pull(Gene)
          uncorrelated.proteins <- rna_protein_lm_stats %>% filter(Coefficient < 0.1) %>% dplyr::pull(Gene)
          uncorrelated.proteins <- uncorrelated.proteins[-which(uncorrelated.proteins %in% c("FGF23", "IL12B"))]
          
          
           
           raw_RNA_vs_protein[which(raw_RNA_vs_protein$Gene %in% correlated.proteins),] %>%
                                          ggplot(., aes(x = log2_RNA, y = log2_Protein)) +
                                                                    geom_jitter(width = 0.5, size = 0.5, alpha = 0.8)+
                                                                    theme_classic() +
                                                                    facet_wrap(~Gene) +
                                                                    geom_smooth(data = . %>% filter(log2_RNA >= -1),
                                                                                method = "lm",
                                                                                linetype = "dashed",
                                                                                color = "red",
                                                                                se = F,
                                                                                size = 0.5) +
                                                                    scale_y_continuous(limits = c(0, 20)) +
                                                                    stat_poly_eq(data = . %>% filter(log2_RNA >= -1),
                                                                                 use_label(c("eq", "R2")),
                                                                                 size = 3)
             RNA_Protein_correlated_genes.plot <- last_plot()
            

             raw_RNA_vs_protein[which(raw_RNA_vs_protein$Gene %in% uncorrelated.proteins),] %>%
                                                                   ggplot(., aes(x = log2_RNA, y = log2_Protein)) +
                                                                   geom_jitter(width = 0.5, size = 0.5, alpha = 0.8)+
                                                                   theme_classic() +
                                                                   facet_wrap(~Gene) +
                                                                   geom_smooth(data = . %>% filter(log2_RNA >= -1),
                                                                               method = "lm",
                                                                               linetype = "dashed",
                                                                               color = "red",
                                                                               se = F,
                                                                               size = 0.5) +
                                                                   scale_y_continuous(limits = c(0, 20)) +
                                                                   stat_poly_eq(data = . %>% filter(log2_RNA >= -1),
                                                                                use_label(c("eq", "R2")),
                                                                                size = 3)
             
                                               
                                                                 
             RNA_Protein_uncorrelated_genes.plot <- last_plot()
             
             
             ggsave(RNA_Protein_correlated_genes.plot, width = 30, height = 25,units = "cm",filename = "Plots/olink/rna_protein_correlated_gene_plot.pdf", limitsize = FALSE)  
             ggsave(RNA_Protein_uncorrelated_genes.plot, width = 30, height = 25,units = "cm",filename = "Plots/olink/rna_protein_uncorrelated_gene_plot.pdf", limitsize = FALSE)  
             
             
             
             
             
             
             
             
          
            raw_RNA_vs_protein %>% filter(Gene %in% c("IFNG", "CXCL8", "LIF", "OSM")) %>%
                      ggplot(., aes(x = log2_RNA, y = log2_Protein)) +
                      geom_jitter(width = 0.5, size = 2, alpha = 0.8)+
                      theme_classic() +
                      scale_x_continuous(limits = c(-7, 15)) +
                      facet_wrap(~Gene)
            
            raw_RNA_vs_protein %>% filter(Gene %in% c("CD40", "CSF1", "IL2RB", "TGFB1")) %>%
              ggplot(., aes(x = log2_RNA, y = log2_Protein)) +
              geom_jitter(width = 0.5, size = 2, alpha = 0.8)+
              theme_classic() +
              scale_x_continuous(limits = c(-7, 15)) +
              facet_wrap(~Gene)
            
            raw_RNA_vs_protein %>% filter(Gene %in% c("TNF", "CCL19", "CX3CL1", "CXCL11")) %>%
              ggplot(., aes(x = log2_RNA, y = log2_Protein)) +
              geom_jitter(width = 0.5, size = 2, alpha = 0.8)+
              theme_classic() +
              scale_x_continuous(limits = c(-7, 15)) +
              facet_wrap(~Gene)
            