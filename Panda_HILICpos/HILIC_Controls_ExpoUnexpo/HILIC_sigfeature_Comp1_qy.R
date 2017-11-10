##########################################
# Propose: To explore significant features derived from
#          different methods (HILICpos)
# Author: Qi Yan
# Date created: 09/12/17
# Date last modified: 09/12/17
##########################################

library(tidyverse)
library(dplyr)
library(VennDiagram)

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo")

# 1.1 significant features from PLSDA, VIP>=2, fold change>=1.5
PLSDA_VIP2 <- read.table("PANDA_output_PLSDA/Stage2/plssignalthresh0.8RSD1/pls_sigfeats_ordered_by_significance.txt", header = T, sep = "\t")
PLSDA_precombine <- PLSDA_VIP2[,c(6,7)]
PLSDA_precombine$method <- "PLSDA"
# 1.2 significant features from rfeSVM, tolerance=30
SVM_tol30 <- read.table("PANDA_output_SVM/Stage2/rfesvmsignalthresh0.8RSD1/rfesvm_sigfeats_ordered_by_significance.txt", header = T, sep = "\t")
SVM_precombine <- SVM_tol30[,c(6,7)]
SVM_precombine$method <- "rfeSVM"
# 1.3 significant features from RF, 71 features
RF_max120 <- read.table("PANDA_output_RF/Stage2/RFsignalthresh0.8RSD1/RF_sigfeats_ordered_by_significance.txt", header = T, sep = "\t")
RF_precombine <- RF_max120[,c(6,7)]
RF_precombine$method <- "RF"
# 1.4 Visualization
combined_feature <- rbind(PLSDA_precombine,SVM_precombine,RF_precombine)
ggplot(data=combined_feature,aes(x=time, y=mz)) +
  geom_point(alpha=.4, size=4, aes(color=method)) +
  ggtitle("Comparison 1: Feature Map") +
  labs(x="Time(s)", y="mz")
## Find intersect between different methods
intersect_PLSDA <- PLSDA_precombine[,c(1:2)]
intersect_SVM <- SVM_precombine[,c(1:2)]
intersect_RF <- RF_precombine[,c(1:2)]
intersect_PLSDA_SVM <- intersect(intersect_PLSDA,intersect_SVM) # N=48
intersect_PLSDA_RF <- intersect(intersect_PLSDA,intersect_RF) # N=31
intersect_RF_SVM <- intersect(intersect_RF,intersect_SVM) # N=20
intersect_all <- intersect(intersect_PLSDA_SVM,intersect_RF) # N=19

intersect_all <- merge(intersect_all,PLSDA_VIP2,by.x="mz",by.y="mz") #Retrive rank back
intersect_all <- intersect_all[,c(1,2,5)]
colnames(intersect_all)[which(names(intersect_all) == "diffexp_rank")] <- "rank_PLSDA"
intersect_all <- merge(intersect_all,SVM_tol30,by.x="mz",by.y="mz")
intersect_all <- intersect_all[,c(1,2,3,6)]
colnames(intersect_all)[which(names(intersect_all) == "diffexp_rank")] <- "rank_SVM"
intersect_all <- merge(intersect_all,RF_max120,by.x="mz",by.y="mz")
intersect_all <- intersect_all[,c(1,2,3,4,7)]
colnames(intersect_all)[which(names(intersect_all) == "diffexp_rank")] <- "rank_RF"

write.table(intersect_all,file="C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_Annotation/intersect_all.txt",sep="\t",row.names=FALSE)
## draw Venn diagram
grid.newpage()
draw.triple.venn(area1 = 115, area2 = 90, area3 = 90, n12 = 48, n23 = 20, n13 = 31, 
                 n123 = 19, category = c("PLSDA", "rfeSVM", "RF"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

############################
# 2. Prepare for mummichog
############################
# 2.1 Use PLSDA VIP>=2
PLSDA_all_VIP2 <- read.table("PANDA_output_PLSDA/Stage2/plssignalthresh0.8RSD1/plsresults_allfeatures.txt", header = T, sep = "\t")
mummichog_PLSDA_VIP2 <- PLSDA_all_VIP2[,c(6,7,4)]
mummichog_PLSDA_VIP2$"p-value" = 0.051
mummichog_PLSDA_VIP2$`p-value`[mummichog_PLSDA_VIP2$VIP>=2] <- 0.004
mummichog_PLSDA_VIP2 <- mummichog_PLSDA_VIP2[order(mummichog_PLSDA_VIP2$'p-value',-mummichog_PLSDA_VIP2$VIP),]
mummichog_PLSDA_VIP2 <- mummichog_PLSDA_VIP2[,c(1,2,4,3)]
filename <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/HILIC_mummichog_PLSDAvip2/input_vip2.txt"
write.table(mummichog_PLSDA_VIP2,file = filename, sep="\t",row.names = F,quote = F)
# 2.2 Use PLSDA VIP>=1.5
PLSDA_all_VIP2 <- read.table("PANDA_output_PLSDA/Stage2/plssignalthresh0.8RSD1/plsresults_allfeatures.txt", header = T, sep = "\t")
mummichog_PLSDA_VIP1.5 <- PLSDA_all_VIP2[,c(6,7,4)]
mummichog_PLSDA_VIP1.5$"p-value" = 0.051
mummichog_PLSDA_VIP1.5$`p-value`[mummichog_PLSDA_VIP1.5$VIP>=1.5] <- 0.004
mummichog_PLSDA_VIP1.5 <- mummichog_PLSDA_VIP1.5[order(mummichog_PLSDA_VIP1.5$'p-value',-mummichog_PLSDA_VIP1.5$VIP),]
mummichog_PLSDA_VIP1.5 <- mummichog_PLSDA_VIP1.5[,c(1,2,4,3)]
filename <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/HILIC_mummichog_PLSDAvip2/input_vip1.5.txt"
write.table(mummichog_PLSDA_VIP1.5,file = filename, sep="\t",row.names = F,quote = F)
