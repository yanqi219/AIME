# library(tidyverse)
library(mixOmics)
library(e1071)
library(ROCR)
library(RColorBrewer)
# library(caret)

time.start <- proc.time()

# Set log and graph direction
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA")
logfile <- paste("Log_",Sys.Date(),".log",sep = "")
logfile <- file(logfile)
sink(logfile) # sink()

pdf_file<-paste("Output_",Sys.Date(),".pdf",sep="")
pdf(file = pdf_file,width=10,height=10)

# Load data
print("Good Luck, Have Fun!")
print("Read metabolomics data")
print("After preprocess")

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")
load(file = "HILIC_control_expo_unexpo_classification.RData")

X <- after.prepro.feature
linkid <- after.prepro.linkid
row.names(sampleID) <- sampleID$SampleID
sampleID$factorcase <- as.numeric(ifelse(sampleID$factorcase == "Exposed", 1, ifelse(sampleID$factorcase == "Unexposed", 0, 99)))
sampleID$factorcase <- as.factor(sampleID$factorcase)
Y <- sampleID$factorcase
levels(Y) = c("control","case")

X = X[ ,apply(X[,1:ncol(X)],2,var) != 0] # remove 0 variance variables

X_caret <- as.data.frame(cbind(X,Y))
X_caret$Y <- as.factor(X_caret$Y)
levels(X_caret$Y) = c("control","case")

rm(after.prepro.feature)
rm(after.prepro.linkid)

# exp.X_caret <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/arcene/ARCENE/arcene_train.data", sep = " ",
#                      colClasses = c(rep("numeric", 10000), "NULL"))
# exp.X_caret$Y <- factor(scan("https://archive.ics.uci.edu/ml/machine-learning-databases/arcene/ARCENE/arcene_train.labels", sep = "\t"))
# levels(exp.X_caret$Y) = c("control","case")
# 
# exp.X_caret = exp.X_caret[ ,apply(exp.X_caret[,1:10000],2,var) != 0] # remove 0 variance variables
# exp.Y <- exp.X_caret$Y
# exp.X <- as.matrix(exp.X_caret[,-9921])

# data(srbct)
# X = srbct$gene
# dim(X)
# Y = srbct$class
# summary(Y)

optimal_comp = FALSE
components = 5
vip_threshold = 2
foldchange_threshold = 0.58
is.log = TRUE
fold_cv = 10
nrepeat = 50
pred.eval.method = "CV"
cpu = 4


###################################################
### I. PCA
###################################################

# PCA
print("PCA analysis")
pca.datExpr = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.datExpr,main = "Variance explained by PCs using all features")
plotIndiv(pca.datExpr, group = Y, ind.names = FALSE, 
          legend = TRUE, title = 'PC score using all features after preprocessing')
print("PCA done")

###################################################
### II. PLSDA
###################################################

# PLSDA
print("PLS-DA analysis")
datExpr.plsda <- plsda(X, Y, ncomp = 10) # set ncomp to 10 for performance assessment later
plotIndiv(datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLS score using all features after preprocessing')

# Plot PLS1 vs. PLS2 with background
background = background.predict(datExpr.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, title = "PLS score using all features after preprocessing (Maximum distance)", ellipse = TRUE,
          legend = TRUE,  background = background)

# Assess the preformance of PLSDA, tuning, select number of component
print("Optimize number of component if necessary")
if (optimal_comp == TRUE){
set.seed(1106)
system.time(
perf.plsda.datExpr <- perf(datExpr.plsda, validation = "Mfold", folds = fold_cv,   ## 5 fold CV
                  progressBar = TRUE, auc = TRUE, nrepeat = nrepeat, cpus = cpu)
)
plot(perf.plsda.datExpr, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

error_rate <- as.data.frame(perf.plsda.datExpr$error.rate$BER)
num_comp <- rownames(error_rate)[which.min(apply(error_rate,MARGIN=1,min))]
dist_method <- as.numeric(which.min(apply(error_rate,MARGIN=2,min)))
# if (optimal_comp == TRUE){
  opt_comp = as.numeric(perf.plsda.datExpr$choice.ncomp[2,dist_method])
}else{
  opt_comp = components
}

# auc.plsda = auroc(datExpr.plsda, roc.comp = 5)

# USing caret
# Compile cross-validation settings
set.seed(1106)
myfolds <- caret::createMultiFolds(X_caret$Y, k = fold_cv, times = nrepeat)              ## 5 fold CV
control <- caret::trainControl("repeatedcv", index = myfolds, selectionFunction = "oneSE")

# Train PLS model
system.time(
  mod1 <- caret::train(Y ~ ., data = X_caret,
                       method = "pls",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = control,
                       preProc = c("zv","center","scale"))
)

# Check CV profile
plot(mod1,main = paste(fold_cv,"-fold CV accuracy over components",sep = ""))
# plot(varImp(mod1), 10, main = "PLS-DA")
pkg <- "package:caret"
detach(pkg, character.only = TRUE)

# select important variables
print("Select important variables")
vip.plsda.datExpr <- vip(datExpr.plsda)
good_feats<-{} 
for(c1 in 1:opt_comp){
  good_feats<-c(good_feats,which(vip.plsda.datExpr[,c1]>=vip_threshold))
}
good_feats<-unique(good_feats)

# fold change
print("Calculate fold change")
fc.Y <- as.vector(as.numeric(Y))-1
fc.X <- t(X)
colnames(fc.X) <- c(1:ncol(fc.X))
f <- DEDS::comp.FC(L=fc.Y, is.log = is.log)
fc <- f(fc.X)
vip_w_fc <- cbind(vip.plsda.datExpr,fc)

# another way to calculate fold change and direction
fc.updown <- as.data.frame(t(aggregate(X_caret[,1:ncol(X_caret)-1],list(X_caret$Y), mean)))
colnames(fc.updown) <- c("control","case")
fc.updown <- fc.updown[-1,]
fc.updown <- as.data.frame(apply(fc.updown,2,function(x) as.numeric(x)))
if(is.log == TRUE){
  fc.updown$foldchange = fc.updown$case-fc.updown$control
}else{
  fc.updown$foldchange = fc.updown$case/fc.updown$control
}

# get vip score for all expr features and plot manhattan plots
print("Calculate vip score and plot manhattan plots")
vip.for.selection<-as.data.frame(apply(vip.plsda.datExpr[,c(1:opt_comp)],1,max))
vip.for.selection <- cbind(vip.for.selection,fc.updown$foldchange)
colnames(vip.for.selection) = c("vip","foldchange")

vip.for.selection <- cbind(linkid,vip.for.selection)
vip.for.selection$group[vip.for.selection$vip < vip_threshold] <- 0
vip.for.selection$group[vip.for.selection$vip >= vip_threshold & vip.for.selection$foldchange <= 0] <- 1 ## down
vip.for.selection$group[vip.for.selection$vip >= vip_threshold & vip.for.selection$foldchange > 0] <- 2  ## up

ggplot2::ggplot(vip.for.selection,aes(x=mz,y=vip)) +
  geom_point(aes(colour = cut(group, c(-Inf,0,1,2,Inf))),size=1,show.legend = FALSE) + 
  scale_fill_hue(c=20, l=20) + 
  scale_color_manual(values = c("#999999","#66CC99","#CC6666")) +
  geom_hline(aes(yintercept = 2),color = "red",size = 1,linetype = "dashed") +
  ggtitle("Type 1 manhattan plot (VIP vs mz)
red: lower in class case & green: higher in class case")

ggplot2::ggplot(vip.for.selection,aes(x=time,y=vip)) +
  geom_point(aes(colour = cut(group, c(-Inf,0,1,2,Inf))),size=1,show.legend = FALSE) + 
  scale_fill_hue(c=20, l=20) + 
  scale_color_manual(values = c("#999999","#66CC99","#CC6666")) +
  geom_hline(aes(yintercept = 2),color = "red",size = 1,linetype = "dashed") +
  ggtitle("Type 2 manhattan plot (VIP vs retention time)
          green: lower in class case & red: higher in class case")
print("PLS-DA phase 1 done")

# Using significant features
print(paste("Assess performance of PLS-DA using significant features, threshold vip score >=",vip_threshold,sep = " "))
good_feats_ordered <- vip.for.selection[order(-vip.for.selection$vip,-abs(vip.for.selection$foldchange)),]
good_feats_ordered <- subset(good_feats_ordered,vip>=vip_threshold & abs(foldchange)>=foldchange_threshold)
good_feats_ordered.name <- {}
for (i in 1:nrow(good_feats_ordered)) {
  good_feats_ordered.name <- c(good_feats_ordered.name,which(colnames(X)==row.names(good_feats_ordered[i,])))
}
sig.X <- X[,c(good_feats_ordered.name)]

#pca using significant features
print("PCA using significant features")
sig.pca.datExpr = pca(sig.X, ncomp = 10, center = TRUE, scale = TRUE)
plot(sig.pca.datExpr,main = "Variance explained by PCs using significant features")
plotIndiv(sig.pca.datExpr, group = Y, ind.names = FALSE, 
          legend = TRUE, title = 'PC score using significant features')

# PLSDA using significant features
print("PLS-DA using significant features")
sig.datExpr.plsda <- plsda(sig.X, Y, ncomp = 10) # set ncomp to 10 for performance assessment later
plotIndiv(sig.datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLS score using significant features')

# Plot PLS1 vs. PLS2 with background using significant features
background = background.predict(sig.datExpr.plsda, comp.predicted=2, dist = "max.dist")
plotIndiv(sig.datExpr.plsda, comp = 1:2,
          group = Y, ind.names = FALSE, title = "PLS score using significant features (Maximum distance)", ellipse = TRUE,
          legend = TRUE,  background = background)

# Assess the preformance of PLSDA using significant features
print("Assess the preformance")
set.seed(1106)
system.time(
  sig.perf.plsda.datExpr <- perf(sig.datExpr.plsda, validation = "Mfold", folds = fold_cv,   ## 5 fold CV
                           progressBar = TRUE, auc = TRUE, nrepeat = nrepeat, cpus = cpu)
)

# get ROC curve using significant features
print("Generating ROC curve using top features on training set")
source("C:/Users/QiYan/Dropbox/AIME/Archive/get_roc.R")
roc.dataA <- t(sig.X)
get_roc(dataA=roc.dataA,classlabels=Y,classifier="svm",kname="radial",
        rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,mainlabel="Training set ROC curve using top features")
print("ROC done")

# get CV accuracy using significant features via SVM
print("get k-fold CV accuracy")
source("C:/Users/QiYan/Dropbox/AIME/Archive/svm_cv.R")
xvec<-{}
yvec<-{}
best_acc<-0
for(i in 2:nrow(roc.dataA)){
  subdata<-t(roc.dataA[1:i,])
  svm_model<-try(svm_cv(v=fold_cv,x=subdata,y=Y,kname="radial",
                        errortype=pred.eval.method,conflevel=95))
  if(is(svm_model,"try-error")){
    svm_model<-NA
    }else{
      xvec<-c(xvec,i)
      yvec<-c(yvec,svm_model$avg_acc)
      if(svm_model$avg_acc>best_acc){
      best_acc<-svm_model$avg_acc
      best_subset<-seq(1,i)
      }
      if(svm_model$avg_acc<best_acc){
        diff_acc<-best_acc-svm_model$avg_acc
        if(diff_acc>50){
          break;
        }
      }
    }
  }

if(pred.eval.method=="CV"){
  ylab_text=paste(pred.eval.method," accuracy (%)",sep="")
  }else{
    if(pred.eval.method=="BER"){
      ylab_text=paste("Balanced accuracy"," (%)",sep="")
      }else{
        ylab_text=paste("AUC"," (%)",sep="")
      }
    }
print(length(yvec))
print(length(xvec))
print(xvec)
print(yvec)
if(length(yvec)>0){
  plot(x=xvec,y=yvec,main="k-fold CV classification accuracy based on forward selection of top features",xlab="Feature index",ylab=ylab_text,type="b",col="brown")
  
  cv_mat<-cbind(xvec,yvec)
  colnames(cv_mat)<-c("Feature Index",ylab_text)
  
  # write.table(cv_mat,file="kfold_cv_mat.txt",sep="\t")
}

cv.acc.sigfeats <- svm_model$avg_acc
print(paste(paste("K-fold CV accuracy is",cv.acc.sigfeats,sep = " "),"%",sep = ""))

#permutation test
print("Calculating permuted CV accuracy")
cv.acc.permut<-{}
subdata <- t(roc.dataA)
cv.acc.permut<-lapply(1:100,function(j){
  rand_order<-sample(1:dim(as.data.frame(Y))[1],size=dim(as.data.frame(Y))[1])
  classlabels_permut<-as.data.frame(Y)[rand_order,]
  classlabels_permut<-as.data.frame(classlabels_permut)
  svm_permut_res<-svm_cv(v=fold_cv,x=subdata,y=classlabels_permut,kname="radial",errortype=pred.eval.method,conflevel=95)
  return(svm_permut_res$avg_acc)
})

cv.acc.permut<-unlist(cv.acc.permut)
cv.acc.permut<-mean(cv.acc.permut,na.rm=TRUE)
cv.acc.permut<-round(cv.acc.permut,2)

print(paste(paste("mean Permuted accuracy is:",cv.acc.permut,sep = " "),"%",sep = ""))

# Plot var plot for permutation test
bar <- as.data.frame(c(cv.acc.sigfeats,cv.acc.permut))
bar$name <- c("Actual","Permutated")
colnames(bar) <- c("value","name")
bar$value <- round(bar$value,2)
ggplot2::ggplot(bar, aes(x=name, y=value)) +
  geom_bar(stat = "identity",aes(fill = as.factor(name)),width = 0.5,show.legend = FALSE) +
  geom_text(aes(label=value),vjust=-0.3) +
  labs(title = "K-fold CV accuracy and permutation test",x = element_blank(),y = "Accuracy (%)")

print("PLS-DA phase 2 done")

# Two-way HCA
print("Two-way HCA")
source("C:/Users/QiYan/Dropbox/AIME/Archive/get_hca.R")
get_hca(data_m = sig.X,classlabels = Y,is.data.znorm = FALSE,clu.method = "bicor")

# Save files
save.plsresults.allfeatures <- cbind(vip.for.selection,t(X))
good_feats_ordered$rank <- 1:nrow(good_feats_ordered)
save.plsresults.sigfeatures <- cbind(good_feats_ordered,roc.dataA)
save.cv.accuracy <- as.data.frame(cv_mat)

###################################################
### III. Prepare for mummichog
###################################################

save.mummichog_PLSDA_VIP2 <- save.plsresults.allfeatures[,1:3]
save.mummichog_PLSDA_VIP2$"p-value" = 0.051
save.mummichog_PLSDA_VIP2$`p-value`[save.mummichog_PLSDA_VIP2$vip>=2] <- 0.004
save.mummichog_PLSDA_VIP2 <- save.mummichog_PLSDA_VIP2[order(save.mummichog_PLSDA_VIP2$`p-value`,-save.mummichog_PLSDA_VIP2$vip),]
save.mummichog_PLSDA_VIP2 <- save.mummichog_PLSDA_VIP2[,c(1,2,4,3)]

###################################################
### IV. Annotation
###################################################

#PLSDA_vip2
annotations_filename<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_KEGG.txt"
annotations_filename_2<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_HMDB.txt"
annotations_filename_3<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/HILICpos_ThermoHFQE_85to1275_mz_range/DBmatches_LipidMaps.txt"

number_significant_digits_rounding<-4

d_sigfeature<-save.plsresults.sigfeatures
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
m_KEGG<-m_KEGG[order(m_KEGG$rank),]
m_HMDB<-merge(a_HMDB,d_sigfeature,by.x="mz",by.y="mz")
m_HMDB<-m_HMDB[order(m_HMDB$rank),]
m_LipidMaps<-merge(a_LipidMaps,d_sigfeature,by.x="mz",by.y="mz")
m_LipidMaps<-m_LipidMaps[order(m_LipidMaps$rank),]

save.KEGG<-subset(m_KEGG,select=c("mz","KEGGID","Name","rank"))
names(save.KEGG)[names(save.KEGG)=="Name"] <- "KEGG_name"
save.HMDB<-subset(m_HMDB,select=c("mz","HMDBID","Name","rank"))
names(save.HMDB)[names(save.HMDB)=="Name"] <- "HMDB_name"
save.LipidMaps<-subset(m_LipidMaps,select=c("mz","LM_ID","Name","rank"))
names(save.LipidMaps)[names(save.LipidMaps)=="Name"] <- "LipidMaps_name"

# save.all<-dplyr::full_join(save.HMDB,save.KEGG,by=c('mz','rank'))
# save.all<-dplyr::full_join(save.all,save.LipidMaps,by=c('mz','rank'))

rm(a_HMDB,a_KEGG,a_LipidMaps,m_HMDB,m_KEGG,m_LipidMaps)

###################################################
### V. Save files
###################################################

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA")
filename <- paste("PLSDA_result_",Sys.Date(),".RData",sep="")
save(save.cv.accuracy,save.HMDB,save.KEGG,save.LipidMaps,save.mummichog_PLSDA_VIP2,save.plsresults.allfeatures,
     save.plsresults.sigfeatures,file = filename)

time.end <- proc.time()
time.run <- time.end - time.start

print("The running time is")
print(time.run)

dev.off()
sink()

###################################################
### code chunk number 11: PLSDA-analysis.Rnw:267-285 (eval = FALSE)
###################################################
## # run internally and saved
## 
## set.seed(2543) # for reproducibility
## # grid of possible keepX values that will be tested for each component
## list.keepX <- c(1:10,  seq(20, 300, 10))  
## 
## t1 = proc.time()  
## tune.splsda.datExpr <- tune.splsda(X, Y, ncomp = 6, validation = 'Mfold', folds = 5, 
##                            progressBar = FALSE, dist = 'max.dist', measure = "BER",
##                           test.keepX = list.keepX, nrepeat = 10)#, cpus = 2)
## t2 = proc.time()
## running_time = t2 - t1; running_time # running time
## 
## error <- tune.splsda.datExpr$error.rate
## ncomp <- tune.splsda.datExpr$choice.ncomp$ncomp
## select.keepX <- tune.splsda.datExpr$choice.keepX[1:ncomp]
## 
## save(tune.splsda.datExpr, ncomp, list.keepX, error, select.keepX, running_time,  file = 'RData/result-datExpr-sPLSDA.RData')


###################################################
### code chunk number 12: PLSDA-analysis.Rnw:288-289
###################################################
load('RData/result-datExpr-sPLSDA.Rdata')


###################################################
### code chunk number 13: PLSDA-analysis.Rnw:292-302 (eval = FALSE)
###################################################
## set.seed(2543) # for reproducibility
## # grid of possible keepX values that will be tested for each component
## list.keepX <- c(1:10,  seq(20, 300, 10))
## 
## t1 = proc.time()  
## tune.splsda.datExpr <- tune.splsda(X, Y, ncomp = 6, validation = 'Mfold', folds = 5, 
##                            progressBar = TRUE, dist = 'max.dist', measure = "BER",
##                           test.keepX = list.keepX, nrepeat = 10)
## t2 = proc.time()
## running_time = t2 - t1; running_time # running time


###################################################
### code chunk number 14: PLSDA-analysis.Rnw:305-306
###################################################
running_time


###################################################
### code chunk number 15: PLSDA-analysis.Rnw:309-310 (eval = FALSE)
###################################################
## error <- tune.splsda.datExpr$error.rate  # error rate per component for the keepX grid


###################################################
### code chunk number 16: PLSDA-analysis.Rnw:314-316 (eval = FALSE)
###################################################
## ncomp <- tune.splsda.datExpr$choice.ncomp$ncomp # optimal number of components based on t-tests
## ncomp


###################################################
### code chunk number 17: PLSDA-analysis.Rnw:319-320
###################################################
ncomp


###################################################
### code chunk number 18: PLSDA-analysis.Rnw:324-326 (eval = FALSE)
###################################################
## select.keepX <- tune.splsda.datExpr$choice.keepX[1:ncomp]  # optimal number of variables to select
## select.keepX


###################################################
### code chunk number 19: PLSDA-analysis.Rnw:328-329
###################################################
select.keepX


###################################################
### code chunk number 20: PLSDA-analysis.Rnw:337-338
###################################################
plot(tune.splsda.datExpr, col = color.jet(6))


###################################################
### code chunk number 21: PLSDA-analysis.Rnw:347-348
###################################################
splsda.datExpr <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX) 


###################################################
### code chunk number 22: PLSDA-analysis.Rnw:354-358
###################################################
plotIndiv(splsda.datExpr, comp = c(1,2),
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on datExpr, comp 1 & 2')


###################################################
### code chunk number 23: PLSDA-analysis.Rnw:365-369
###################################################
plotIndiv(splsda.datExpr, comp = c(1,3),
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on datExpr, comp 1 & 3')


###################################################
### code chunk number 24: PLSDA-analysis.Rnw:377-378
###################################################
auc.splsda = auroc(splsda.datExpr, roc.comp = 2)


###################################################
### code chunk number 25: PLSDA-analysis.Rnw:383-384
###################################################
auc.splsda = auroc(splsda.datExpr, roc.comp = ncomp)


###################################################
### code chunk number 26: PLSDA-analysis.Rnw:390-398
###################################################
set.seed(40) # for reproducibility, only when the `cpus' argument is not used

t1 = proc.time()  
perf.datExpr <- perf(splsda.datExpr, validation = "Mfold", folds = 5,
                   dist = 'max.dist', nrepeat = 10,
                   progressBar = FALSE) 
t2 = proc.time()
running_time = t2 - t1; running_time


###################################################
### code chunk number 27: PLSDA-analysis.Rnw:402-405
###################################################
# perf.datExpr  # lists the different outputs
perf.datExpr$error.rate
plot(perf.datExpr, col = color.mixo(5))


###################################################
### code chunk number 28: PLSDA-analysis.Rnw:412-420
###################################################
par(mfrow=c(1,3))
plot(perf.datExpr$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.datExpr$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)
plot(perf.datExpr$features$stable[[3]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 3', las =2)
par(mfrow=c(1,1))


###################################################
### code chunk number 29: PLSDA-analysis.Rnw:426-433
###################################################
# here we match the selected variables to the stable features
ind.match = match(selectVar(splsda.datExpr, comp = 1)$name, 
                  names(perf.datExpr$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.datExpr$features$stable[[1]][ind.match])

data.frame(selectVar(splsda.datExpr, comp = 1)$value, Freq)


###################################################
### code chunk number 30: PLSDA-analysis.Rnw:440-442
###################################################
plotLoadings(splsda.datExpr, comp = 1, title = 'Loadings on comp 1', 
             contrib = 'max', method = 'mean')


###################################################
### code chunk number 31: PLSDA-analysis.Rnw:446-448
###################################################
plotLoadings(splsda.datExpr, comp = 2, title = 'Loadings on comp 2', 
             contrib = 'max', method = 'mean')


###################################################
### code chunk number 32: PLSDA-analysis.Rnw:451-453
###################################################
plotLoadings(splsda.datExpr, comp = 3, title = 'Loadings on comp 3', 
             contrib = 'max', method = 'mean')


###################################################
### code chunk number 33: PLSDA-analysis.Rnw:459-461
###################################################
plotLoadings(splsda.datExpr, comp = 1, title = 'Loadings on comp 1', 
             contrib = 'min', method = 'mean')


###################################################
### code chunk number 34: PLSDA-analysis.Rnw:467-468
###################################################
cim(splsda.datExpr)


###################################################
### code chunk number 35: PLSDA-analysis.Rnw:471-472
###################################################
cim(splsda.datExpr, comp=1, title ="Component 1")


###################################################
### code chunk number 36: PLSDA-analysis.Rnw:478-479
###################################################
plotArrow(splsda.datExpr, legend=T)


###################################################
### code chunk number 37: PLSDA-analysis.Rnw:488-489
###################################################
sessionInfo()


###################################################
### code chunk number 38: PLSDA-analysis.Rnw:493-495
###################################################
# the end of that chapter. Command not to be run
Stangle('PLSDA-analysis.Rnw', encoding = 'utf8')


