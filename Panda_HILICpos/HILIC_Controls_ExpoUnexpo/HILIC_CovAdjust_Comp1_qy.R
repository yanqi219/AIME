library(dplyr)
library(plyr)
library(xmsPANDA)
library(mixOmics)

##pre-process metabolic data

class <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_classlabels_control_expo_unexpo.txt"
feature <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt"
outloc <-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA_residual"

# ready_for_regression<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file=feature,parentoutput_dir=outloc,class_labels_file=class,num_replicates=3,feat.filt.thresh=NA,
#                                       summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.5,group.missing.thresh=0.8,
#                                       log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=TRUE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
#                                       samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="zeros",featselmethod=NA)
ready_for_regression<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file=feature,parentoutput_dir=outloc,class_labels_file=class,num_replicates=3,feat.filt.thresh=NA,
                                      summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.5,group.missing.thresh=0.8,
                                      log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
                                      samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="halfdatamin",featselmethod=NA)

feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
na_count <-sapply(feature, function(y) sum(is.na(y)))
summary(na_count)

setwd("C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range")
load(file = "HILIC_class.rda")

# setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")
# feature <- read.table("HILIC_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt",sep="\t",header=TRUE)

##extract covariates
sampleID <- HILIC_class[HILIC_class$casecontrol..Factor.1.=="0",]
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

tail <- data.frame(V3=c("_1"))
# tail <- data.frame(V3=c("_1","_3","_5"))
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
# fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,13:ncol(feature_w_cov)]) ~ as.factor(sex), na.action = na.exclude) # residual_1
# fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,13:ncol(feature_w_cov)]) ~ as.factor(sex)+as.factor(maternal_edu), na.action = na.exclude) # residual_2
# fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,13:ncol(feature_w_cov)]) ~ as.factor(sex)+as.factor(maternal_edu)+as.factor(pregcompl)
#                   +ttcbl+as.factor(maternal_raceeth), na.action = na.exclude) # residual_3
fit_feature <- lm(data = feature_w_cov, as.matrix(feature_w_cov[,13:ncol(feature_w_cov)]) ~ as.factor(sex)+as.factor(maternal_age)+as.factor(maternal_raceeth)+
                    as.factor(maternal_edu)+lengthgestation+as.factor(pregcompl)+ttcbl+as.factor(preterm)+as.factor(usborn), na.action = na.exclude)
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

##save data file
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")
save_sampleID <- sampleID[,c(1:2)]
write.table(save_sampleID,file="HILIC_residuals_classlabels_control_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)
write.table(wide_save_residual,file="HILIC_residuals_control_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)

######################### Get class table with covariates for univariates regression ####################
##extract covariates
sampleID <- HILIC_class[HILIC_class$casecontrol..Factor.1.=="0",]
sampleID <- sampleID[sampleID$Exposure_category..Factor2.=="G"|sampleID$Exposure_category..Factor2.=="H",c(1,11:21)]
colnames(sampleID)
sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="G"] <- "Exposed"
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="H"] <- "Unexposed"
sampleID <- sampleID[!is.na(sampleID$SampleID),]

##get subjects' sample ID
sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
rownames(sampleID) <- c(1:489)

tail <- data.frame(V3=c("_1","_3","_5"))
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

class <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_classlabels_control_expo_unexpo.txt"
class <- read.table(class,sep="\t",header=TRUE)
ordered_sampleID <- inner_join(class,sampleID,by=c("SampleID","factorcase"))

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")
write.table(ordered_sampleID,file="HILIC_withcov_classlabels_control_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)

####################################################### Check covariates ########################################################

## Distribution of covariates
count(feature_w_cov,var = c("factorcase","sex"))
count(feature_w_cov,var = c("factorcase","maternal_age"))
count(feature_w_cov,var = c("factorcase","maternal_raceeth"))
count(feature_w_cov,var = c("factorcase","maternal_edu"))
count(feature_w_cov,var = c("factorcase","preterm"))
count(feature_w_cov,var = c("factorcase","usborn"))

## Exposure status ~ cov
regression_data <- feature_w_cov
regression_data$factorcase <- ifelse(regression_data$factorcase=='Exposed',1,ifelse(regression_data$factorcase=='Unexposed',2,99))
regression_data$factorcase <- as.factor(regression_data$factorcase)
expo_cov_logit <- glm(factorcase ~ as.factor(sex)+as.factor(birthyear)+as.factor(maternal_age)+as.factor(maternal_raceeth)+as.factor(maternal_edu)+
                        lengthgestation+as.factor(pregcompl)+ttcbl+as.factor(preterm)+as.factor(usborn),data=regression_data,family=binomial())
summary(expo_cov_logit) # birthyear, maternal race, maternal edu, pregcompl, ttcbl significant

## Metabolites ~ exposure + cov
fit <- lm(data = regression_data, met_1 ~ as.factor(factorcase)+as.factor(sex)+as.factor(birthyear)+as.factor(maternal_age)+as.factor(maternal_raceeth)+
            lengthgestation+as.factor(maternal_edu)+as.factor(pregcompl)+ttcbl+as.factor(preterm)+as.factor(usborn), na.action = na.exclude)
summary(fit) # maternal edu significant
p_value <- save_link
p_value$p <- NA
for (i in 13:ncol(regression_data)){
met_expo_lm <- lm(data = regression_data, regression_data[,i] ~ as.factor(factorcase)+as.factor(sex)+as.factor(birthyear)+as.factor(maternal_age)+as.factor(maternal_raceeth)+
                    lengthgestation+as.factor(maternal_edu)+as.factor(pregcompl)+ttcbl+as.factor(preterm)+as.factor(usborn), na.action = na.exclude)
p_value[i-12,4]<-summary(met_expo_lm)$coefficients[2,4]
}
p_value <- p_value[order(p_value$p),]
p_value$Bonferroni <- p.adjust(p_value$p, method = "bonferroni")
p_value$BH <- p.adjust(p_value$p, method = "BH")
