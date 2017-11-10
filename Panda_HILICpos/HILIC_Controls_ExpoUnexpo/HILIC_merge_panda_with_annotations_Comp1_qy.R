
library(dplyr)

#PLSDA_vip2
panda_filename<-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA/Stage2/plssignalthresh0.8RSD1/pls_sigfeats_ordered_by_significance.txt"
#SVM_tol30
panda_filename<-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_SVM/Stage2/rfesvmsignalthresh0.8RSD1/rfesvm_sigfeats_ordered_by_significance.txt"
#RF_max120
panda_filename<-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_RF/Stage2/RFsignalthresh0.8RSD1/RF_sigfeats_ordered_by_significance.txt"

annotations_filename<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_KEGG.txt"
annotations_filename_2<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_HMDB.txt"
annotations_filename_3<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_LipidMaps.txt"

parentoutput_dir="C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_Annotation"

number_significant_digits_rounding<-4

#check_original<-read.table("C:/Users/Qi/Dropbox/AIME/PNS_Ritz/Autism_PNSstudy_HILICpos_xmsPANDA_analysis/archive/merge_annotation_results/plsVIP2_w_mummichog/PANDAresults_merged_with_annotations.txt",sep="\t",header=TRUE)

########################
setwd(parentoutput_dir)
plsVIP2_w_all <- read.table("plsVIP2_w_all.txt", header = T, sep = "\t")
plsVIP2_w_HMDB <- read.table("plsVIP2_w_HMDB.txt", header = T, sep = "\t")
plsVIP2_w_KEGG <- read.table("plsVIP2_w_KEGG.txt", header = T, sep = "\t")
plsVIP2_w_LipidMaps <- read.table("plsVIP2_w_LipidMaps.txt", header = T, sep = "\t")

RF_max120_w_all <- read.table("RF_max120_w_all.txt", header = T, sep = "\t")
RF_max120_w_HMDB <- read.table("RF_max120_w_HMDB.txt", header = T, sep = "\t")
RF_max120_w_KEGG <- read.table("RF_max120_w_KEGG.txt", header = T, sep = "\t")
RF_max120_w_LipidMaps <- read.table("RF_max120_w_LipidMaps.txt", header = T, sep = "\t")

SVMtol30_w_all <- read.table("SVMtol30_w_all.txt", header = T, sep = "\t")
SVMtol30_w_HMDB <- read.table("SVMtol30_w_HMDB.txt", header = T, sep = "\t")
SVMtol30_w_KEGG <- read.table("SVMtol30_w_KEGG.txt", header = T, sep = "\t")
SVMtol30_w_LipidMaps <- read.table("SVMtol30_w_LipidMaps.txt", header = T, sep = "\t")

########################

setwd(parentoutput_dir)

d_sigfeature<-read.table(panda_filename,sep="\t",header=TRUE)
dim(d_sigfeature)

a_KEGG<-read.table(annotations_filename,sep="\t",header=TRUE)
a_HMDB<-read.table(annotations_filename_2,sep="\t",header=TRUE)
a_LipidMaps<-read.table(annotations_filename_3,sep="\t",header=TRUE)

a_KEGG$mz<-round(a_KEGG$mz,number_significant_digits_rounding)
a_KEGG<-a_KEGG[,-c(20:2032)]
a_HMDB$mz<-round(a_HMDB$mz,number_significant_digits_rounding)
a_HMDB<-a_HMDB[,-c(20:2032)]
a_LipidMaps$mz<-round(a_LipidMaps$mz,number_significant_digits_rounding)
a_LipidMaps<-a_LipidMaps[,-c(20:2032)]

d_sigfeature$mz<-round(d_sigfeature$mz,number_significant_digits_rounding)

m_KEGG<-merge(a_KEGG,d_sigfeature,by.x="mz",by.y="mz")      ## Merged data
m_KEGG<-m_KEGG[order(m_KEGG$diffexp_rank),]
m_HMDB<-merge(a_HMDB,d_sigfeature,by.x="mz",by.y="mz")
m_HMDB<-m_HMDB[order(m_HMDB$diffexp_rank),]
m_LipidMaps<-merge(a_LipidMaps,d_sigfeature,by.x="mz",by.y="mz")
m_LipidMaps<-m_LipidMaps[order(m_LipidMaps$diffexp_rank),]

mall_KEGG<-subset(m_KEGG,select=c("mz","KEGGID","Name","diffexp_rank"))
names(mall_KEGG)[names(mall_KEGG)=="Name"] <- "KEGG_name"
mall_HMDB<-subset(m_HMDB,select=c("mz","HMDBID","Name","diffexp_rank"))
names(mall_HMDB)[names(mall_HMDB)=="Name"] <- "HMDB_name"
mall_LipidMaps<-subset(m_LipidMaps,select=c("mz","LM_ID","Name","diffexp_rank"))
names(mall_LipidMaps)[names(mall_LipidMaps)=="Name"] <- "LipidMaps_name"

m_all<-full_join(mall_HMDB,mall_KEGG,by=c('mz','diffexp_rank'))
m_all<-full_join(m_all,mall_LipidMaps,by=c('mz','diffexp_rank'))

#PLSDA_vip2
write.table(m_HMDB,file="plsVIP2_w_HMDB.txt",sep="\t",row.names=FALSE)
write.table(m_KEGG,file="plsVIP2_w_KEGG.txt",sep="\t",row.names=FALSE)
write.table(m_LipidMaps,file="plsVIP2_w_LipidMaps.txt",sep="\t",row.names=FALSE)
write.table(m_all,file="plsVIP2_w_all.txt",sep="\t",row.names=FALSE)
#SVM_tol30
write.table(m_HMDB,file="SVMtol30_w_HMDB.txt",sep="\t",row.names=FALSE)
write.table(m_KEGG,file="SVMtol30_w_KEGG.txt",sep="\t",row.names=FALSE)
write.table(m_LipidMaps,file="SVMtol30_w_LipidMaps.txt",sep="\t",row.names=FALSE)
write.table(m_all,file="SVMtol30_w_all.txt",sep="\t",row.names=FALSE)
#RF_max120
write.table(m_HMDB,file="RF_max120_w_HMDB.txt",sep="\t",row.names=FALSE)
write.table(m_KEGG,file="RF_max120_w_KEGG.txt",sep="\t",row.names=FALSE)
write.table(m_LipidMaps,file="RF_max120_w_LipidMaps.txt",sep="\t",row.names=FALSE)
write.table(m_all,file="RF_max120_w_all.txt",sep="\t",row.names=FALSE)

#Check the intersect features within all three methods
intersect_all <- read.table("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_Annotation/intersect_all.txt", header = T, sep = "\t")
d_sigfeature<-intersect_all

d_sigfeature$mz<-round(d_sigfeature$mz,number_significant_digits_rounding)

m_KEGG<-merge(a_KEGG,d_sigfeature,by.x="mz",by.y="mz")   ## Merged data
m_HMDB<-merge(a_HMDB,d_sigfeature,by.x="mz",by.y="mz")
m_LipidMaps<-merge(a_LipidMaps,d_sigfeature,by.x="mz",by.y="mz")