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



#load data:
setwd("/pi/john.harris-umw/frisoli/RStudio/Contact_Derm_Analysis/5_6_22_Analysis/")

load("Data/Olink/filtered.olink.data.Rdata")    
      
      #filter out attempted concentrated samples since it didn't even do anything:
      filtered.olink.data <- filtered.olink.data[-which(str_detect(filtered.olink.data$sample_num, pattern = "mf-")),]


      #Blister:
      blister_data_for_pca <- filtered.olink.data %>% filter(sample_type == "Blister") %>% pivot_wider(names_from = Protein, values_from = Concentration)
      
      #PCA using prcomp:
      meta.data <- blister_data_for_pca[,c(1:which(colnames(blister_data_for_pca) == "IL8")-1)]
      
#Continue Combination Biomarker Sensitivity Analysis:      
load("Data/Olink/olink.protein.classifier.sensitivity.Rdata")
      
    
    all.olink.proteins <- colnames(blister_data_for_pca)[c(which(colnames(blister_data_for_pca) == "IL8"):length(colnames(blister_data_for_pca)))]      
    
    single.prot.classifiers <- multi.prot.classifiers %>% filter(Protein_number == "single")
    dual.prot.classifiers <- multi.prot.classifiers %>% filter(Protein_number == "dual")
    
    ranked.single.prots <- single.prot.classifiers %>% arrange(desc(Day2sensitivity)) %>% dplyr::pull(Protein)
    best.single.prots <- ranked.single.prots[1:10]
    
    
    #Three-protein combinations as reliable classifiers:
    Three.prot.combos <- as.data.frame(t(as.data.frame(combn(all.olink.proteins, 3), cols = 3)))
    
    Three.prot.combos <- Three.prot.combos %>% filter(V1 %in% best.single.prots | V2 %in% best.single.prots | V3 %in% best.single.prots) %>% t() %>% as.data.frame()

    
    #Find which combos to still analyze and which have been complated already:
        
          combos.to.analyze <- as.data.frame(t(Three.prot.combos))
          combos.to.analyze$sorted.proteins <- NA
          
              for (c in 1:nrow(combos.to.analyze)) {
                  combo.proteins<- str_c(as.character(combos.to.analyze[c,c(1:3)]), collapse = ".")
                  combos.to.analyze$sorted.proteins[c] <- sapply(combo.proteins, FUN = function(x) str_c(sort(unlist(str_split(x, pattern = "\\."))), collapse = "."))
              }
          
          combos.to.analyze$protein1 <- sapply(combos.to.analyze$sorted.proteins, FUN = function(x) sort(unlist(str_split(x, pattern = "\\.")))[1])
          combos.to.analyze$protein2 <- sapply(combos.to.analyze$sorted.proteins, FUN = function(x) sort(unlist(str_split(x, pattern = "\\.")))[2])
          combos.to.analyze$protein3 <- sapply(combos.to.analyze$sorted.proteins, FUN = function(x) sort(unlist(str_split(x, pattern = "\\.")))[3])
          
          
          finished <- multi.prot.classifiers %>% filter(Protein_number == "triple")
        
          finished$protein1 <- sapply(finished$Protein, FUN = function(x) sort(unlist(str_split(x, pattern = "\\*|\\/")))[1])
          finished$protein2 <- sapply(finished$Protein, FUN = function(x) sort(unlist(str_split(x, pattern = "\\*|\\/")))[2])
          finished$protein3 <- sapply(finished$Protein, FUN = function(x) sort(unlist(str_split(x, pattern = "\\*|\\/")))[3])
          
          for (c in 1:nrow(finished)) {
            combo.proteins <- str_c(as.character(finished[c,c(5:7)]), collapse = ".")
            finished$sorted.proteins[c] <- sapply(combo.proteins, FUN = function(x) str_c(sort(unlist(str_split(x, pattern =  "\\."))), collapse = "."))
          }
          
          finished.stats <- finished %>% group_by(sorted.proteins) %>% summarise(count = n())
          
          if (any(finished.stats$count) <7 ) {
      
            finished.stats <- finished.stats %>% filter(count <7)
            finished <- finished[-which(finished$sorted.proteins %in% finished.stats$sorted.proteins),] 
          }

          finished.unique <- finished[-which(duplicated(finished$sorted.proteins)),]
          if (nrow(finished.unique) < nrow(combos.to.analyze)) {
            
            to.process <- combos.to.analyze[-which(combos.to.analyze$sorted.proteins %in% unique(finished$sorted.proteins)),]
            to.process <-  as.data.frame(t(to.process[,1:3]))
          }

          #Finished is completed triple.prot.classifiers df:
          triple.prot.classifiers <- finished[,1:4]
          
  #Process remaining combos:
          if (ncol(to.process) != 0) {
            Three.prot.combos <- to.process
          }
          
          if (ncol(to.process) == 0) {
            triple.prot.classifiers <- multi.prot.classifiers[0,1:3]
          }
          
          
          for (r in 1:ncol(Three.prot.combos)) {
            
            prots <- Three.prot.combos[,r]
            
            prot.data <- cbind(meta.data[,which(colnames(meta.data) %in% c("lesion", "Lesion_severity_score"))], blister_data_for_pca[,which(colnames(blister_data_for_pca) %in% prots)])
            
            prot.data <- prot.data %>% mutate(multiplicative = eval(as.symbol(prots[1]))*eval(as.symbol(prots[2]))*eval(as.symbol(prots[3])),
                                              divisional1 = eval(as.symbol(prots[1]))*eval(as.symbol(prots[2]))/eval(as.symbol(prots[3])),
                                              divisional2 = eval(as.symbol(prots[1]))*eval(as.symbol(prots[3]))/eval(as.symbol(prots[2])),
                                              divisional3 = eval(as.symbol(prots[2]))*eval(as.symbol(prots[3]))/eval(as.symbol(prots[1])),
                                              
                                              divisional4 = eval(as.symbol(prots[1]))/(eval(as.symbol(prots[2]))*eval(as.symbol(prots[3]))),
                                              divisional5 = eval(as.symbol(prots[2]))/(eval(as.symbol(prots[1]))*eval(as.symbol(prots[3]))),
                                              divisional6 = eval(as.symbol(prots[3]))/(eval(as.symbol(prots[1]))*eval(as.symbol(prots[2]))))
            
            
            colnames(prot.data)[-c(1:5)] <- c(paste(prots[1],"*",prots[2],"*", prots[3], sep = ""),
                                              paste(prots[1],"*",prots[2],"/", prots[3], sep = ""),
                                              paste(prots[1],"*",prots[3],"/", prots[2], sep = ""),
                                              paste(prots[2],"*",prots[3],"/", prots[1], sep = ""),
                                              
                                              paste(prots[1],"/(",prots[2],"*", prots[3], ")",sep = ""),
                                              paste(prots[2],"/(",prots[1],"*", prots[3], ")",sep = ""),
                                              paste(prots[3],"/(",prots[1],"*", prots[2], ")",sep = ""))
            
            
            #100 percent specificity threshold:
            max.lesion.vals <- prot.data %>% group_by(lesion) %>% summarise_at(vars(c(paste(prots[1],"*",prots[2],"*", prots[3], sep = ""),
                                                                                      paste(prots[1],"*",prots[2],"/", prots[3], sep = ""),
                                                                                      paste(prots[1],"*",prots[3],"/", prots[2], sep = ""),
                                                                                      paste(prots[2],"*",prots[3],"/", prots[1], sep = ""),
                                                                                      
                                                                                      paste(prots[1],"/(",prots[2],"*", prots[3], ")",sep = ""),
                                                                                      paste(prots[2],"/(",prots[1],"*", prots[3], ")",sep = ""),
                                                                                      paste(prots[3],"/(",prots[1],"*", prots[2], ")",sep = ""))), funs(max)) %>% 
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
                prot.classifier <- data.frame(Protein = combo,
                                              Day2sensitivity = Day2sensitivity,
                                              Day4sensitivity = Day4sensitivity)
                
                triple.prot.classifiers <- bind_rows(triple.prot.classifiers, prot.classifier)
              
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
          
          save(multi.prot.classifiers, file = "Data/Olink/olink.protein.classifier.sensitivity.Rdata")
          
