library(dplyr)
library(plyr)
library(xmsPANDA)
library(mixOmics)

##pre-process metabolic data

class <- "C:/Users/Qi/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input/C18_classlabels_control_expo_unexpo.txt"
feature <- "C:/Users/Qi/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input/C18_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt"
outloc <-"C:/Users/Qi/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA"

ready_for_regression<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file=feature,parentoutput_dir=outloc,class_labels_file=class,num_replicates=3,feat.filt.thresh=NA,
                                      summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.5,group.missing.thresh=0.8,
                                      log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
                                      samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="halfdatamin",featselmethod=NA)

feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
na_count <-sapply(feature, function(y) sum(is.na(y)))
summary(na_count)

row.names(feature) <- c(paste("met_",1:nrow(feature),sep = ""))
after.prepro.linkid <- feature[,1:2]
after.prepro.feature <- t(feature[,-c(1:2)])
after.prepro.feature <- after.prepro.feature[ order(row.names(after.prepro.feature)), ]

setwd("C:/Users/Qi/Dropbox/AIME/PNS_Ritz/C18neg_ThermoHFQE_85to1275_mz_range")
load(file = "C18_class.rda")

# setwd("C:/Users/Qi/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input")
# feature <- read.table("C18_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt",sep="\t",header=TRUE)

##extract covariates
sampleID <- C18_class[C18_class$casecontrol..Factor.1.=="0",]
sampleID <- sampleID[sampleID$Exposure_category..Factor2.=="G"|sampleID$Exposure_category..Factor2.=="H",c(1,11:21)]
colnames(sampleID)
sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="G"] <- "Exposed"
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="H"] <- "Unexposed"
sampleID <- sampleID[!is.na(sampleID$SampleID),]

##get subjects' sample ID
# sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
# rownames(sampleID) <- c(1:489)
rownames(sampleID) <- c(1:163)

tail <- data.frame(V3=c("_2"))
tail <- data.frame(rep(tail$V3,times=163))
colnames(tail) <- "V3"

sampleID <- cbind(sampleID,tail)
sampleID$V4 <- paste(sampleID$SampleID,sampleID$V3,sep = "")

sampleID <- subset(sampleID,select=-c(V3,SampleID))

sampleID <- plyr::rename(sampleID,c('V4'='SampleID'))
sampleID <- plyr::rename(sampleID,c('Exposure_category..Factor2.'='factorcase'))
colnames(sampleID)
sampleID <- sampleID[c("SampleID","factorcase","sex","birthyear","maternal_age","maternal_raceeth","maternal_edu"
                       ,"lengthgestation","pregcompl","ttcbl","preterm","usborn")]

# Recode some covariates due tp spase data
sampleID$maternal_age[which(sampleID$maternal_age==5)] <- 4
sampleID$maternal_raceeth[which(sampleID$maternal_raceeth==4)] <- 3
sampleID$maternal_raceeth[which(sampleID$maternal_raceeth==5)] <- 3
sampleID$maternal_edu[which(sampleID$maternal_edu==2)] <- 1
sampleID$maternal_edu[which(sampleID$maternal_edu==3)] <- 2
sampleID$maternal_edu[which(sampleID$maternal_edu==4)] <- 3
sampleID$maternal_edu[which(sampleID$maternal_edu==5)] <- 4

# Check outliers
pca.datExpr = pca(after.prepro.feature, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.datExpr,main = "Variance explained by PCs using all features")
plotIndiv(pca.datExpr, group = sampleID$factorcase, ind.names = FALSE, 
          legend = TRUE, title = 'PC score using all features after preprocessing')
plotIndiv(pca.datExpr, group = sampleID$factorcase, ind.names = TRUE, 
          legend = TRUE, title = 'PC score using all features after preprocessing')
# Remove outliers
after.prepro.feature <- after.prepro.feature[-which(row.names(after.prepro.feature)=="M2249033_VS100001F16757_2"),]
sampleID <- sampleID[-which(sampleID$SampleID=="M2249033_VS100001F16757_2"),]
feature <- feature[,-which(colnames(feature)=="M2249033_VS100001F16757_2")]

# reformat feature table from wide to long, create unique metabolite id number, link with covariates table
feature <- cbind(c(1:nrow(feature)),feature)
colnames(feature)[1] <- 'metablite'
feature$metablite <- paste("met_",feature$metablite,sep = "")

save_link <- feature[,c("metablite","mz","time")]
save_feature <- subset(feature,select = -c(mz,time))

## transpose
n <- save_feature$metablite

long_save_feature <- as.data.frame(t(save_feature[,-1]))
colnames(long_save_feature) <- n
id <- row.names(long_save_feature)
long_save_feature <- cbind(id,long_save_feature)
long_save_feature <- plyr::rename(long_save_feature,c('id'='SampleID'))
row.names(long_save_feature) <- c(1:nrow(long_save_feature))

# merge covariates data wirh feature data
feature_w_cov <- merge(sampleID,long_save_feature, by = "SampleID")

# adjust for covariates on each metabolits and then get residuals
# fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,13:ncol(feature_w_cov)]) ~ as.factor(sex)+as.factor(maternal_age)+as.factor(maternal_edu)+as.factor(pregcompl)
#                  +ttcbl+as.factor(maternal_raceeth), na.action = na.exclude) # residual_3
fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,13:ncol(feature_w_cov)]) ~ as.factor(maternal_age)+as.factor(maternal_edu)+
                    as.factor(maternal_raceeth), na.action = na.exclude) # residual_3
residual_feature <- as.matrix(residuals(fit_feature),nrow = dim(feature_w_cov)[1],ncol = dim(save_feature)[1])
save_residual <- as.data.frame(residual_feature)
save_residual <- cbind(feature_w_cov$SampleID,save_residual)
save_residual <- save_residual[order(save_residual$`feature_w_cov$SampleID`),]
row.names(save_residual) <- c(1:nrow(save_residual))

# reformat residual from long to wide
n <- save_residual$`feature_w_cov$SampleID`

wide_save_residual <- as.data.frame(t(save_residual[,-1]))
colnames(wide_save_residual) <- n
id <- row.names(wide_save_residual)
wide_save_residual <- cbind(id,wide_save_residual)
names(wide_save_residual)[names(wide_save_residual) == 'id'] <- 'metablite'
row.names(wide_save_residual) <- c(1:nrow(wide_save_residual))

# link wide_save_residual with save_link, retrive m/z and time back
wide_save_residual <- merge(save_link, wide_save_residual, by = "metablite")
wide_save_residual <- wide_save_residual[,-1]
wide_save_residual <- wide_save_residual[order(wide_save_residual$mz,wide_save_residual$time),]
row.names(wide_save_residual) <- c(1:nrow(wide_save_residual))

# # replace na with 0
# wide_save_residual<-replace(wide_save_residual,is.na(wide_save_residual),0)

# remove na
wide_save_residual<-wide_save_residual[sapply(wide_save_residual, function(x) !any(is.na(x)))]
complete_sub <- row.names(t(wide_save_residual[,-c(1:2)]))
sampleID <- subset(sampleID, sampleID$SampleID %in% complete_sub)
after.prepro.feature <- subset(after.prepro.feature,row.names(after.prepro.feature) %in% complete_sub)

##save data file
setwd("C:/Users/Qi/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input")
save(sampleID, wide_save_residual, file = "C18_control_expo_unexpo_residual_nonorm_WGCNA.RData") ## for WGCNA
save_sampleID <- sampleID[,c(1:2)]
write.table(save_sampleID,file="C18_residuals_classlabels_control_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)
write.table(wide_save_residual,file="C18_residuals_control_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)
save(sampleID, after.prepro.feature,after.prepro.linkid, file = "C18_control_expo_unexpo_classification_nonorm.RData")