library(dplyr)

setwd("C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range")

class_example <- read.table("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_input/classlabels_case_control_noexposure.txt",header = T)
feature_example <- read.table("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_input/ftrsmzcalib_combat_ordered_case_control_noexposure.txt",header = T)

load(file = "HILIC_feature_combat.rda")
load(file = "HILIC_feature_raw.rda")
load(file = "HILIC_class.rda")
load(file = "HILIC_sample_link.rda")
load(file = "HILIC_match_kegg.rda")
load(file = "HILIC_match_hmdb.rda")
load(file = "HILIC_match_lipidmaps.rda")
load(file = "HILIC_replicates_corr.rda")

#############################
#Unexposed cases and controls
#############################

##extract cases and controls ID from class dataset
sampleID <- HILIC_class[HILIC_class$Exposure_category..Factor2.=="H",c(1,10)]
colnames(sampleID)
sampleID$casecontrol..Factor.1.[sampleID$casecontrol..Factor.1.==1] <- "case"
sampleID$casecontrol..Factor.1.[sampleID$casecontrol..Factor.1.==0] <- "control"
sampleID <- sampleID[!is.na(sampleID$SampleID),]

##get subjects' 3 replicates sample ID

sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
rownames(sampleID) <- c(1:351)

tail <- data.frame(V3=c("_2","_4","_6"))
tail <- data.frame(rep(tail$V3,times=117))
colnames(tail) <- "V3"

sampleID <- cbind(sampleID,tail)
sampleID$V4 <- paste(sampleID$SampleID,sampleID$V3,sep = "")

sampleID <- sampleID[,c("V4","casecontrol..Factor.1.")]

colnames(sampleID) <- c("SampleID","factorcase")

##get the feature table

link <- sampleID
colnames(link) <- c("Sample.ID","factorcase")
link <- merge(link,HILIC_sample_link,by="Sample.ID")
link <- link[,1:3]
link$V4 <- ".mzXML"
link$File.Name <- paste(link$File.Name,link$V4,sep = "")
addon <- data.frame(Sample.ID=c("mz","time"),factorcase=c("mz","time"),File.Name=c("mz","time"),V4=c("mz","time"))
link <- rbind(addon,link)

featurename <- as.vector(link$File.Name)
twomore <- c("mz","time")
featurename <- c(twomore,featurename)

# names(HILIC_feature_combat) <- lapply(HILIC_feature_combat[1, ], as.character)

feature <- HILIC_feature_combat[,colnames(HILIC_feature_combat) %in% featurename]    ##keep only wanted columns

colname_feature <- colnames(feature)

long_feature <- as.data.frame(colname_feature)
colnames(long_feature) <- "File.Name"
long_feature <- merge(long_feature,link,by="File.Name")

class <- long_feature[3:nrow(long_feature),2:3]
colnames(class) <- c("SampleID","factorcase")

long_feature <- long_feature[,1:2]
wide_feature <- as.data.frame(t(long_feature))

names(feature) <- lapply(wide_feature[2, ], as.character)


# wide_feature <- as.data.frame(t(colname_feature))
# 
# colnames(wide_feature) <- "V3"
# wide_feature <- merge(wide_feature,link,by="V3")

# wide_feature <- wide_feature[,1:2]
# long_feature <- as.data.frame(t(wide_feature))
# 
# long_feature <- long_feature[2,]
# 
# colnames(feature) <- colnames(long_feature)
# 
# feature <- rbind(long_feature,feature)
# feature <- feature[-c(2),]
# rownames(feature) <- 1:nrow(feature)    ##change first row in feature data
# 
# addon <- data.frame(V1="SampleID",V2="factorcase")
# sampleID <- rbind(addon,sampleID)

##save data file

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_input")

write.table(class,file="HILIC_classlabels_case_control_noexposure.txt",sep="\t",row.names = F,quote = F)
write.table(feature,file="HILIC_ftrsmzcalib_combat_ordered_case_control_noexposure.txt",sep = "\t",row.names = F,quote = F)

#############################
#Air pollution exposed cases and controls
#############################

##extract cases and controls ID from class dataset
sampleID <- HILIC_class[HILIC_class$Exposure_category..Factor2.=="G",c(1,10)]
colnames(sampleID)
sampleID$casecontrol..Factor.1.[sampleID$casecontrol..Factor.1.==1] <- "case"
sampleID$casecontrol..Factor.1.[sampleID$casecontrol..Factor.1.==0] <- "control"
sampleID <- sampleID[!is.na(sampleID$SampleID),]

##get subjects' 3 replicates sample ID

sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
rownames(sampleID) <- c(1:660)

tail <- data.frame(V3=c("_1","_3","_5"))
tail <- data.frame(rep(tail$V3,times=220))
colnames(tail) <- "V3"

sampleID <- cbind(sampleID,tail)
sampleID$V4 <- paste(sampleID$SampleID,sampleID$V3,sep = "")

sampleID <- sampleID[,c("V4","casecontrol..Factor.1.")]

colnames(sampleID) <- c("SampleID","factorcase")

##get the feature table

link <- sampleID
colnames(link) <- c("Sample.ID","factorcase")
link <- merge(link,HILIC_sample_link,by="Sample.ID")
link <- link[,1:3]
link$V4 <- ".mzXML"
link$File.Name <- paste(link$File.Name,link$V4,sep = "")
addon <- data.frame(Sample.ID=c("mz","time"),factorcase=c("mz","time"),File.Name=c("mz","time"),V4=c("mz","time"))
link <- rbind(addon,link)

featurename <- as.vector(link$File.Name)
twomore <- c("mz","time")
featurename <- c(twomore,featurename)

# names(HILIC_feature_combat) <- lapply(HILIC_feature_combat[1, ], as.character)

feature <- HILIC_feature_combat[,colnames(HILIC_feature_combat) %in% featurename]    ##keep only wanted columns

colname_feature <- colnames(feature)

long_feature <- as.data.frame(colname_feature)
colnames(long_feature) <- "File.Name"
long_feature <- merge(long_feature,link,by="File.Name")

class <- long_feature[3:nrow(long_feature),2:3]
colnames(class) <- c("SampleID","factorcase")

long_feature <- long_feature[,1:2]
wide_feature <- as.data.frame(t(long_feature))

names(feature) <- lapply(wide_feature[2, ], as.character)

##save data file

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Exposed_CasesControls/PANDA_input")

write.table(class,file="HILIC_classlabels_case_control_exposure.txt",sep="\t",row.names = F,quote = F)
write.table(feature,file="HILIC_ftrsmzcalib_combat_ordered_case_control_exposure.txt",sep = "\t",row.names = F,quote = F)

#############################
#Air pollution exposed controls and unexposed controls
#############################

##extract cases and controls ID from class dataset
sampleID <- HILIC_class[HILIC_class$casecontrol..Factor.1.=="0",]
sampleID <- sampleID[sampleID$Exposure_category..Factor2.=="G"|sampleID$Exposure_category..Factor2.=="H",c(1,11)]
colnames(sampleID)
sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="G"] <- "Exposed"
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="H"] <- "Unexposed"
sampleID <- sampleID[!is.na(sampleID$SampleID),]

##get subjects' 3 replicates sample ID

sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
rownames(sampleID) <- c(1:489)

tail <- data.frame(V3=c("_1","_3","_5"))
tail <- data.frame(rep(tail$V3,times=163))
colnames(tail) <- "V3"

sampleID <- cbind(sampleID,tail)
sampleID$V4 <- paste(sampleID$SampleID,sampleID$V3,sep = "")

sampleID <- sampleID[,c("V4","Exposure_category..Factor2.")]

colnames(sampleID) <- c("SampleID","factorcase")

##get the feature table

link <- sampleID
colnames(link) <- c("Sample.ID","factorcase")
link <- merge(link,HILIC_sample_link,by="Sample.ID")
link <- link[,1:3]
link$V4 <- ".mzXML"
link$File.Name <- paste(link$File.Name,link$V4,sep = "")
addon <- data.frame(Sample.ID=c("mz","time"),factorcase=c("mz","time"),File.Name=c("mz","time"),V4=c("mz","time"))
link <- rbind(addon,link)

featurename <- as.vector(link$File.Name)
twomore <- c("mz","time")
featurename <- c(twomore,featurename)

# names(HILIC_feature_combat) <- lapply(HILIC_feature_combat[1, ], as.character)

feature <- HILIC_feature_combat[,colnames(HILIC_feature_combat) %in% featurename]    ##keep only wanted columns

colname_feature <- colnames(feature)

long_feature <- as.data.frame(colname_feature)
colnames(long_feature) <- "File.Name"
long_feature <- merge(long_feature,link,by="File.Name")

class <- long_feature[3:nrow(long_feature),2:3]
colnames(class) <- c("SampleID","factorcase")

long_feature <- long_feature[,1:2]
wide_feature <- as.data.frame(t(long_feature))

names(feature) <- lapply(wide_feature[2, ], as.character)

##save data file

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")

write.table(class,file="HILIC_classlabels_control_expo_unexpo.txt",sep="\t",row.names = F,quote = F)
write.table(feature,file="HILIC_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)

#############################
#Air pollution exposed cases and unexposed cases
#############################

##extract cases and controls ID from class dataset
sampleID <- HILIC_class[HILIC_class$casecontrol..Factor.1.=="1",]
sampleID <- sampleID[sampleID$Exposure_category..Factor2.=="G"|sampleID$Exposure_category..Factor2.=="H",c(1,11)]
colnames(sampleID)
sampleID$Exposure_category..Factor2.<-as.character(sampleID$Exposure_category..Factor2.)
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="G"] <- "Exposed"
sampleID$Exposure_category..Factor2.[sampleID$Exposure_category..Factor2.=="H"] <- "Unexposed"
sampleID <- sampleID[!is.na(sampleID$SampleID),]

##get subjects' 3 replicates sample ID

sampleID <- sampleID[rep(row.names(sampleID),3),]
sampleID <- sampleID[order(sampleID$SampleID),]
rownames(sampleID) <- c(1:522)

tail <- data.frame(V3=c("_1","_3","_5"))
tail <- data.frame(rep(tail$V3,times=174))
colnames(tail) <- "V3"

sampleID <- cbind(sampleID,tail)
sampleID$V4 <- paste(sampleID$SampleID,sampleID$V3,sep = "")

sampleID <- sampleID[,c("V4","Exposure_category..Factor2.")]

colnames(sampleID) <- c("SampleID","factorcase")

##get the feature table

link <- sampleID
colnames(link) <- c("Sample.ID","factorcase")
link <- merge(link,HILIC_sample_link,by="Sample.ID")
link <- link[,1:3]
link$V4 <- ".mzXML"
link$File.Name <- paste(link$File.Name,link$V4,sep = "")
addon <- data.frame(Sample.ID=c("mz","time"),factorcase=c("mz","time"),File.Name=c("mz","time"),V4=c("mz","time"))
link <- rbind(addon,link)

featurename <- as.vector(link$File.Name)
twomore <- c("mz","time")
featurename <- c(twomore,featurename)

# names(HILIC_feature_combat) <- lapply(HILIC_feature_combat[1, ], as.character)

feature <- HILIC_feature_combat[,colnames(HILIC_feature_combat) %in% featurename]    ##keep only wanted columns

colname_feature <- colnames(feature)

long_feature <- as.data.frame(colname_feature)
colnames(long_feature) <- "File.Name"
long_feature <- merge(long_feature,link,by="File.Name")

class <- long_feature[3:nrow(long_feature),2:3]
colnames(class) <- c("SampleID","factorcase")

long_feature <- long_feature[,1:2]
wide_feature <- as.data.frame(t(long_feature))

names(feature) <- lapply(wide_feature[2, ], as.character)

##save data file

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Cases_ExpoUnexpo/PANDA_input")

write.table(class,file="HILIC_classlabels_case_expo_unexpo.txt",sep="\t",row.names = F,quote = F)
write.table(feature,file="HILIC_ftrsmzcalib_combat_ordered_case_expo_unexpo.txt",sep = "\t",row.names = F,quote = F)


