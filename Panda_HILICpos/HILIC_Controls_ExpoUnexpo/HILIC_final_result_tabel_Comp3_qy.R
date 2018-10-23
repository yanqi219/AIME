library(tidyverse)
library(zoo)
library(xlsx)

# Note: if we use mummichog version 1, then we use MTBNK+PLSDAvip2fc0.58 features, otherwise we use PLSDAvip2fc0 features

####################
# annotation
####################
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Annotation/sigfeature_annotation/")
annotation_HMDB <- read.csv(file = "HMDB/Stage5.csv", sep = ",", header = T)
annotation_KEGG <- read.csv(file = "KEGG/Stage5.csv", sep = ",", header = T)
annotation_LipidMaps <- read.csv(file = "LipidMaps/Stage5.csv", sep = ",", header = T)

####################
# feature selection - mummichog server
####################
# setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_output_PLSDA/")
# load("Res_PLSDA_result_2018-05-01_vip2fc0.58.RData")
# feature_PLSDA <- save.plsresults.allfeatures[which(save.plsresults.allfeatures$vip>=2),]
# feature_PLSDA <- feature_PLSDA[order(-feature_PLSDA$vip),]
# rm(save.cv.accuracy, save.HMDB, save.KEGG, save.LipidMaps, save.mummichog_PLSDA_VIP2, save.plsresults.allfeatures, save.plsresults.sigfeatures)

####################
# pathway analysis - mummichog server
####################
# setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/HILIC_mummichogOL_vip2fc0/tables/")
# mummichog_input <- read.table(file = "userInputData.txt", sep = "\t", header = T)
# mummichog_pathway <- read.xlsx(file = "mcg_pathwayanalysis_.xlsx", 1, header = T)
# mummichog_empirical <- read.table(file = "ListOfEmpiricalCompounds.tsv", sep = "\t", header = T)
# mummichog_link <- read.xlsx(file = "mcg_pathwayanalysis_.xlsx", 2, header = T)    # second sheet copied from the web, is the link between feature and name
# 
# # rearrange the link table
# flag <- which(grepl(pattern = "E", x = as.character(mummichog_link[,1])))
# temp <- mummichog_link[flag,]
# temp <- temp[,c(1,2,4)]
# colnames(temp) <- c("EmpiricalCompound", "CompoundName", "KEGGID")
# mummichog_link$EmpiricalCompound <- na.locf(mummichog_link$EmpiricalCompound)
# mummichog_link <- mummichog_link[-flag,]
# mummichog_link <- merge(x = mummichog_link, y = temp, by = "EmpiricalCompound", all.x = T)
# rm(temp,flag)
# 
# # rearrange pathway data from wide to long
# mummichog_pathway <- mummichog_pathway[which(mummichog_pathway$p.value<=0.05),]  #add overlap size or not
# mummichog_pathway_expand <- data.frame()
# for(i in 1:nrow(mummichog_pathway)){
#   temp <- mummichog_pathway[i,]
#   temp.compound <- as.character(temp$overlap_EmpiricalCompounds..id.)
#   temp.compound <- as.data.frame(unlist(strsplit(temp.compound, ",")))
#   temp.rest <- temp[,c(1:4)]
#   temp.merge <- merge(x = temp.rest, y = temp.compound, all.y = T)
#   mummichog_pathway_expand <- rbind(mummichog_pathway_expand, temp.merge)
# }
# colnames(mummichog_pathway_expand) <- c("pathway", "overlap_size", "pathway_size", "p.value", "EmpiricalCompound")
# rm(temp, temp.compound, temp.merge, temp.rest)
# 
# # combine pathway with link
# mummichog_pathway_complete <- merge(mummichog_link, mummichog_pathway_expand, by = "EmpiricalCompound", all = T)
# mummichog_pathway_complete <- mummichog_pathway_complete[-which(is.na(mummichog_pathway_complete$pathway)),]
# mummichog_pathway_complete <- mummichog_pathway_complete[order(mummichog_pathway_complete$pathway),]
# 
# names(mummichog_pathway_complete)[names(mummichog_pathway_complete) == "Input.m.z"] <- "mz"
# names(mummichog_pathway_complete)[names(mummichog_pathway_complete) == "Retention.time"] <- "time"
# 
# mummichog_pathway_complete$time <- round(mummichog_pathway_complete$time, 1)
# mummichog_pathway_complete$mz <- round(as.numeric(as.character(mummichog_pathway_complete$mz)), 4)

####################
# feature selection - mummichog version 1
####################
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_output_PLSDA/")
load("Res_PLSDA_result_2018-05-01_vip2fc0.58.RData")
feature_MTBNK <- read.csv(file = "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/mummichogMTBNK_input_vip2fc0.58.txt",
                          sep = "\t", header = T)
feature_MTBNK <- feature_MTBNK[which(feature_MTBNK$p.value <= 0.04), c(1:2)]
feature_PLSDA <- merge(feature_MTBNK, save.plsresults.allfeatures, by = c("mz","time"), all.X = T)
feature_PLSDA <- feature_PLSDA[order(-feature_PLSDA$vip),]
rm(save.cv.accuracy, save.HMDB, save.KEGG, save.LipidMaps, save.mummichog_PLSDA_VIP2, save.plsresults.allfeatures, save.plsresults.sigfeatures, feature_MTBNK)

####################
# pathway analysis - mummichog version 1
####################
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/HILIC_mummichog/HILIC_mummichogMTBNK_vip2fc0.58/tsv/")
mummichog_pathway <- read.xlsx(file = "mcg_pathwayanalysis_HILIC_mummichogMTBNK_vip2fc0.58.xlsx", 2, header = T)

mummichog_pathway_complete <- mummichog_pathway
names(mummichog_pathway_complete)[names(mummichog_pathway_complete) == "m.z"] <- "mz"
mummichog_pathway_complete$mz <- round(as.numeric(as.character(mummichog_pathway_complete$mz)), 4)
mummichog_pathway_complete <- mummichog_pathway_complete[,-4]

####################
# combine three datasets
####################
feature_PLSDA$time <- round(feature_PLSDA$time, 1)
annotation_HMDB$time <- round(annotation_HMDB$time, 1)
annotation_KEGG$time <- round(annotation_KEGG$time, 1)
annotation_LipidMaps$time <- round(annotation_LipidMaps$time, 1)

feature_PLSDA$mz <- round(feature_PLSDA$mz, 4)
annotation_HMDB$mz <- round(annotation_HMDB$mz, 4)
annotation_KEGG$mz <- round(annotation_KEGG$mz, 4)
annotation_LipidMaps$mz <- round(annotation_LipidMaps$mz, 4)

# merge only confidence score >= 2

# KEGG first
annotation_KEGG_highconf <- annotation_KEGG[which(annotation_KEGG$Confidence >= 1),]
feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_KEGG_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
# feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_KEGG_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
flag <- which(!is.na(feature_w_annotation$chemical_ID))
temp.KEGG <- feature_w_annotation[flag,]
feature_w_annotation <- feature_w_annotation[-flag,]   # the rest pass to HMDB

# HMDB
annotation_HMDB_highconf <- annotation_HMDB[which(annotation_HMDB$Confidence >= 1),]
# feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_HMDB_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_HMDB_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
flag <- which(!is.na(feature_w_annotation$chemical_ID))
temp.HMDB <- feature_w_annotation[flag,]
feature_w_annotation <- feature_w_annotation[-flag,]   

# LipidMaps
annotation_LipidMaps_highconf <- annotation_LipidMaps[which(annotation_LipidMaps$Confidence >= 1),]
# feature_w_annotation <- merge(feature_PLSDA[,c(1:4)], annotation_LipidMaps_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
feature_w_annotation <- merge(feature_w_annotation[,c(1:4)], annotation_LipidMaps_highconf[,c(1,2,5,6,10)], by = c("mz", "time"), all = T)
feature_w_annotation <- feature_w_annotation[-which(is.na(feature_w_annotation$vip)),]
flag <- which(!is.na(feature_w_annotation$chemical_ID))
temp.LipidMaps <- feature_w_annotation[flag,]
feature_w_annotation <- feature_w_annotation[-flag,]
# add back
feature_w_annotation <- rbind(feature_w_annotation, temp.HMDB, temp.KEGG, temp.LipidMaps)
feature_w_annotation <- feature_w_annotation[order(-feature_w_annotation$vip),]
rm(temp.HMDB, temp.KEGG, temp.LipidMaps)

# add pathway

# merge_all <- merge(feature_w_annotation, mummichog_pathway_complete, by = c("mz", "time"), all = T)  #this is for version 2!!!
# merge_all <- merge_all[,-c(8,10,16:18)]
merge_all <- merge(feature_w_annotation, mummichog_pathway_complete, by = c("mz"), all = T)
# remove those non-significant features picked by mummichog
merge_all <- merge_all[-which(is.na(merge_all$vip)),]
merge_all <- merge_all[order(-merge_all$vip),]

## Now the problem is that we need to use mummichog 1 instead of mummichog 2!!