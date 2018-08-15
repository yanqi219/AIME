# library(tidyverse)
library(mixOmics)
library(e1071)
library(ROCR)
library(RColorBrewer)
# library(caret)

#################################
# Set universal parameters
#################################
is.residual = TRUE

optimal_comp = FALSE
components = 5
vip_threshold = 2
foldchange_threshold = 0
is.log = TRUE
fold_cv = 10
nrepeat = 50
pred.eval.method = "BER"
cluster.method = "dist" #"dist" or "bicor" #Residual can only use "dist"
cpu = 4

dir.folder <- "C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/"

#################################
# Start!
#################################

time.start <- proc.time()

# Set log and graph direction
dir.temp <- paste(dir.folder,"PANDA_output_PLSDA",sep = "")
setwd(dir.temp)

if(is.residual == TRUE){
  logfile <- paste("Res_Log_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".log",sep = "")
  pdf_file<-paste("Res_Output_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".pdf",sep="")
}else{
  logfile <- paste("Log_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".log",sep = "")
  pdf_file<-paste("Output_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".pdf",sep="")
}

logfile <- file(logfile)
sink(logfile) # sink()

pdf(file = pdf_file,width=10,height=10)

# Load data
print("Good Luck, Have Fun!")
print("Read metabolomics data")
print("After preprocess")

dir.temp <- paste(dir.folder,"PANDA_input",sep = "")
setwd(dir.temp)

if(is.residual == FALSE){
  load(file = "C18_control_expo_unexpo_classification_nonorm.RData")
  
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
}else{
  # For fold change
  
  load(file = "C18_control_expo_unexpo_classification_nonorm.RData")
  X <- after.prepro.feature
  linkid <- after.prepro.linkid
  row.names(sampleID) <- sampleID$SampleID
  sampleID$factorcase <- as.numeric(ifelse(sampleID$factorcase == "Exposed", 1, ifelse(sampleID$factorcase == "Unexposed", 0, 99)))
  sampleID$factorcase <- as.factor(sampleID$factorcase)
  Y <- sampleID$factorcase
  levels(Y) = c("control","case")
  
  X = X[ ,apply(X[,1:ncol(X)],2,var) != 0] # remove 0 variance variables
  
  fc.X_caret <- as.data.frame(cbind(X,Y))
  fc.X_caret$Y <- as.factor(fc.X_caret$Y)
  levels(fc.X_caret$Y) = c("control","case")
  rm(after.prepro.feature, after.prepro.linkid)
  
  load(file = "C18_control_expo_unexpo_residual_nonorm_WGCNA.RData")
  
  row.names(wide_save_residual) <- c(paste("met_",1:nrow(wide_save_residual),sep = ""))
  linkid <- wide_save_residual[,1:2]
  X <- t(as.matrix(wide_save_residual[,-c(1:2)]))
  row.names(sampleID) <- sampleID$SampleID
  sampleID$factorcase <- as.numeric(ifelse(sampleID$factorcase == "Exposed", 1, ifelse(sampleID$factorcase == "Unexposed", 0, 99)))
  sampleID$factorcase <- as.factor(sampleID$factorcase)
  Y <- sampleID$factorcase
  levels(Y) = c("control","case")
  
  X = X[ ,apply(X[,1:ncol(X)],2,var) != 0] # remove 0 variance variables
  
  X_caret <- as.data.frame(cbind(X,Y))
  X_caret$Y <- as.factor(X_caret$Y)
  levels(X_caret$Y) = c("control","case")
  
  rm(wide_save_residual)
}

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
fc.updown <- as.data.frame(t(aggregate(fc.X_caret[,1:ncol(fc.X_caret)-1],list(fc.X_caret$Y), mean)))
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
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="M/Z",y="VIP") +
  ggtitle("Type 1 manhattan plot (VIP vs mz)
          red: lower in class case & green: higher in class case")

ggplot2::ggplot(vip.for.selection,aes(x=time,y=vip)) +
  geom_point(aes(colour = cut(group, c(-Inf,0,1,2,Inf))),size=1,show.legend = FALSE) + 
  scale_fill_hue(c=20, l=20) + 
  scale_color_manual(values = c("#999999","#66CC99","#CC6666")) +
  geom_hline(aes(yintercept = 2),color = "red",size = 1,linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Retention Time",y="VIP") +
  ggtitle("Type 2 manhattan plot (VIP vs retention time)
          green: lower in class case & red: higher in class case")


volcano <- vip.for.selection[,c(1:4)]
volcano$group <- ifelse(volcano$vip>=vip_threshold&volcano$foldchange>=foldchange_threshold,2,
                        ifelse(volcano$vip>=vip_threshold&volcano$foldchange<=-foldchange_threshold,1,0))
volcano$size <- ifelse(volcano$group == 0,0,1)

ggplot2::ggplot(volcano,aes(x=foldchange,y=vip)) +
  geom_point(aes(colour=cut(group, c(-Inf,0,1,2,Inf)),size=size),show.legend = FALSE) + 
  xlim(-3.8,3.8) +
  scale_fill_hue(c=20, l=20) + 
  scale_color_manual(values = c("#999999","#66CC99","#CC6666")) + 
  scale_size_continuous(range = c(1,3)) +
  geom_hline(aes(yintercept = 2),color = "red",size = 0.8,linetype = "dashed") +
  geom_vline(aes(xintercept = 0),color = "red",size = 0.8,linetype = "dashed") +
  labs(x="Log2(Fold Change)",y="VIP") +
  # theme_minimal()
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Volcano Plot")

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
get_hca(data_m = sig.X,classlabels = Y,is.data.znorm = FALSE,clu.method = cluster.method)

# Save files
save.plsresults.allfeatures <- cbind(vip.for.selection,t(X))
good_feats_ordered$rank <- 1:nrow(good_feats_ordered)
save.plsresults.sigfeatures <- cbind(good_feats_ordered,roc.dataA)
save.cv.accuracy <- as.data.frame(cv_mat)

print(paste("Number of significant features:",nrow(save.plsresults.sigfeatures),sep = " "))

###################################################
### III. Prepare for mummichog
###################################################

save.mummichog_PLSDA_VIP2 <- save.plsresults.allfeatures[,1:4]
save.mummichog_PLSDA_VIP2$"p-value" = 0.051
save.mummichog_PLSDA_VIP2$`p-value`[save.mummichog_PLSDA_VIP2$vip>=vip_threshold&abs(save.mummichog_PLSDA_VIP2$foldchange)>=foldchange_threshold] <- 0.004
save.mummichog_PLSDA_VIP2 <- save.mummichog_PLSDA_VIP2[order(save.mummichog_PLSDA_VIP2$`p-value`,-save.mummichog_PLSDA_VIP2$vip),]
save.mummichog_PLSDA_VIP2 <- save.mummichog_PLSDA_VIP2[,c(1,2,5,4)]

###################################################
### IV. Annotation
###################################################

#PLSDA_vip2
annotations_filename<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/C18neg_ThermoHFQE_85to1275_mz_range/DBmatches_KEGG.txt"
annotations_filename_2<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/C18neg_ThermoHFQE_85to1275_mz_range/DBmatches_HMDB.txt"
annotations_filename_3<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/C18neg_ThermoHFQE_85to1275_mz_range/DBmatches_LipidMaps.txt"

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

a_KEGG$time<-round(a_KEGG$time,2)
a_KEGG<-a_KEGG[,-c(20:2032)]
a_HMDB$time<-round(a_HMDB$time,2)
a_HMDB<-a_HMDB[,-c(20:2032)]
a_LipidMaps$time<-round(a_LipidMaps$time,2)
a_LipidMaps<-a_LipidMaps[,-c(20:2032)]

d_sigfeature$mz<-round(d_sigfeature$mz,number_significant_digits_rounding)
d_sigfeature$time<-round(d_sigfeature$time,2)

m_KEGG<-merge(a_KEGG,d_sigfeature,by.x=c("mz","time"),by.y=c("mz","time"))      ## Merged data
m_KEGG<-m_KEGG[order(m_KEGG$rank),]
m_HMDB<-merge(a_HMDB,d_sigfeature,by.x=c("mz","time"),by.y=c("mz","time"))
m_HMDB<-m_HMDB[order(m_HMDB$rank),]
m_LipidMaps<-merge(a_LipidMaps,d_sigfeature,by.x=c("mz","time"),by.y=c("mz","time"))
m_LipidMaps<-m_LipidMaps[order(m_LipidMaps$rank),]

save.KEGG<-subset(m_KEGG,select=c("mz","time","KEGGID","Name","rank"))
names(save.KEGG)[names(save.KEGG)=="Name"] <- "KEGG_name"
save.HMDB<-subset(m_HMDB,select=c("mz","time","HMDBID","Name","rank"))
names(save.HMDB)[names(save.HMDB)=="Name"] <- "HMDB_name"
save.LipidMaps<-subset(m_LipidMaps,select=c("mz","time","LM_ID","Name","rank"))
names(save.LipidMaps)[names(save.LipidMaps)=="Name"] <- "LipidMaps_name"

# save.all<-dplyr::full_join(save.HMDB,save.KEGG,by=c('mz','rank'))
# save.all<-dplyr::full_join(save.all,save.LipidMaps,by=c('mz','rank'))

rm(a_HMDB,a_KEGG,a_LipidMaps,m_HMDB,m_KEGG,m_LipidMaps)

###################################################
### V. Save files
###################################################

dir.temp <- paste(dir.folder,"PANDA_output_PLSDA",sep = "")
setwd(dir.temp)

if(is.residual == TRUE){
  filename <- paste("Res_PLSDA_result_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".RData",sep="")
}else{
  filename <- paste("PLSDA_result_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".RData",sep="")
}
save(save.cv.accuracy,save.HMDB,save.KEGG,save.LipidMaps,save.mummichog_PLSDA_VIP2,save.plsresults.allfeatures,
     save.plsresults.sigfeatures,file = filename)

time.end <- proc.time()
time.run <- time.end - time.start

print("The running time is")
print(time.run)

dev.off()
sink()

# Create mummichog files

load(filename)
mummichog_name <- paste("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/mummichog_input","_vip",vip_threshold,"fc",foldchange_threshold,"_",Sys.Date(),".txt",sep="")
write.table(save.mummichog_PLSDA_VIP2,file=mummichog_name,sep = "\t",row.names = F,quote = F)

###################################################
### Get correlated list with MetabNet
###################################################
#Get target list and feature table for MetabNet
dir.temp <- paste(dir.folder,"PANDA_output_PLSDA",sep = "")
setwd(dir.temp)
if(is.residual == TRUE){
  filename <- paste("Res_PLSDA_result_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".RData",sep="")
}else{
  filename <- paste("PLSDA_result_",Sys.Date(),"_vip",vip_threshold,"fc",foldchange_threshold,".RData",sep="")
}
load(filename)
target_list <- save.plsresults.sigfeatures[,1:2]
colnames(target_list) <- c("mz","time")

dir.temp <- paste(dir.folder,"PANDA_input",sep = "")
setwd(dir.temp)
load(file = "C18_control_expo_unexpo_residual_nonorm_WGCNA.RData")

target_list_name <- paste("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input/C18_metabnet_target","_vip",vip_threshold,"fc",foldchange_threshold,".txt",sep="")
write.table(target_list,file=target_list_name,sep = "\t",row.names = F,quote = F)
feature_name <- paste("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input/C18_metabnet_feature","_vip",vip_threshold,"fc",foldchange_threshold,".txt",sep="")
write.table(wide_save_residual,file=feature_name,sep = "\t",row.names = F,quote = F)

#Conduct MetabNet
library(MetabNet)

#output location
outloc<-"C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/MetabNet/vip2fc0.58"

allowWGCNAThreads()
net_res<-metabnet(feature_table_file=feature_name,
                  target.metab.file=target_list_name, sig.metab.file=NA,
                  parentoutput_dir=outloc,
                  class_labels_file=NA,cor.method="spearman",abs.cor.thresh=0.6,cor.fdrthresh=0.2,
                  cor.fdrmethod="BH",
                  target.mzmatch.diff=0.01,
                  target.rtmatch.diff=0.01,max.cor.num=150,num_replicates=1,summarize.replicates=FALSE,
                  all.missing.thresh=0, rep.max.missing.thresh = 0,
                  group.missing.thresh=NA,
                  log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=NA,
                  networktype="complete", summary.na.replacement="bpca",samplermindex=NA,
                  net_node_colors=c("yellow", "green"), net_node_shapes=c("rectangle","circle"),net_edge_colors=c("red","blue"),net_legend=FALSE,
                  netrandseed =1106,num_nodes=6)

#Retrive corelated features and save as mummichog file
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/MetabNet/vip2fc0.58/Stage2")
MetabNet.sigcorr.mztime <- read.table(file = "significant_correlations_targeted_matrix_mzlabels.txt", header = T, sep = "\t")
MetabNet.FDR <- read.table(file = "correlation_FDR.txt", header = T, sep = "\t")
MetabNet.list <- {}
for(i in 1:nrow(MetabNet.sigcorr.mztime)){
  temp <- colnames(MetabNet.sigcorr.mztime)[which(MetabNet.sigcorr.mztime[i,]!=0)]
  MetabNet.list <- append(MetabNet.list, temp)
}
MetabNet.list <- substring(MetabNet.list,2)
MetabNet.list <- append(MetabNet.list,row.names(MetabNet.FDR))
MetabNet.list <- unique(MetabNet.list)
MetabNet.list <- as.data.frame(MetabNet.list)
for(i in 1:nrow(MetabNet.list)){
  MetabNet.list[i,2] <- strsplit(as.character(MetabNet.list[i,1]),"_")[[1]][1]
  MetabNet.list[i,3] <- strsplit(as.character(MetabNet.list[i,1]),"_")[[1]][2]
}
MetabNet.list <- MetabNet.list[,2:3]
colnames(MetabNet.list) <- c("mz","time")
MetabNet.list$flag <- 1
MetabNet.list$mz <- as.numeric(MetabNet.list$mz)
MetabNet.list$time <- as.numeric(MetabNet.list$time)

save.mummichog_MetabNet<-merge(x=save.mummichog_PLSDA_VIP2,y=MetabNet.list,by=c("mz","time"),all.x=TRUE)
save.mummichog_MetabNet$`p-value`[which(save.mummichog_MetabNet$flag==1)] <- 0.04
save.mummichog_MetabNet <- save.mummichog_MetabNet[order(save.mummichog_MetabNet$`p-value`),-5]

mummichog_name <- paste("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/mummichogMTBNK_input","_vip",vip_threshold,"fc",foldchange_threshold,"_",Sys.Date(),".txt",sep="")
write.table(save.mummichog_MetabNet,file=mummichog_name,sep = "\t",row.names = F,quote = F)

# save pathway_metaboanalyst file
pathway_metaboanalyst <- save.mummichog_MetabNet[,c(1,3,4)]
colnames(pathway_metaboanalyst) <- c("m.z","p.value","t.score")
write.table(pathway_metaboanalyst,file="C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/pathway_metaboanalyst.txt",sep = "\t",row.names = F,quote = F)

###################################################
### VI. Create box plots
###################################################
library(xlsx)

# load comp3: unexpo case vs. control
load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Non_Exposed_CasesControls/PANDA_input/C18_case_control_noexposure_residual3_WGCNA.RData")
comp3.wide_save_residual <- wide_save_residual
comp3.sampleID <- sampleID
comp3.wide_save_residual$mz <- round(comp3.wide_save_residual$mz,4)
comp3.wide_save_residual$time <- round(comp3.wide_save_residual$time,2)
rm(sampleID)
rm(wide_save_residual)
# load comp4: unexpo vs. expo case
load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Cases_ExpoUnexpo/PANDA_input/C18_case_expo_unexpo_residual3_WGCNA.RData")
comp4.wide_save_residual <- wide_save_residual
comp4.sampleID <- sampleID
comp4.wide_save_residual$mz <- round(comp4.wide_save_residual$mz,4)
comp4.wide_save_residual$time <- round(comp4.wide_save_residual$time,2)
rm(sampleID)
rm(wide_save_residual)

load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input/C18_control_expo_unexpo_residual3_WGCNA.RData")
load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/Res_PLSDA_result_2018-01-23.RData")
wide_save_residual$mz <- round(wide_save_residual$mz,4)
wide_save_residual$time <- round(wide_save_residual$time,2)
save.plsresults.sigfeatures$mz <- round(save.plsresults.sigfeatures$mz,4)
save.plsresults.sigfeatures$time <- round(save.plsresults.sigfeatures$time,2)

# annotations_filename<-"C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/C18neg_ThermoHFQE_85to1275_mz_range/DBmatches_KEGG.txt"
# a_KEGG<-read.table(annotations_filename,sep="\t",header=TRUE)

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichog_0113/Res_C18_PLSDAvip2_model4/tsv")
mummichog.pathway <- read.xlsx("mcg_pathwayanalysis_C18_PLSDAvip2.xlsx",sheetName = "Sheet1")
mummichog.MatchMetab <- read.xlsx("_tentative_featurematch_C18_PLSDAvip2.xlsx",sheetName = "_tentative_featurematch_")

mummichog.sig.pathway <-as.data.frame(mummichog.pathway[which(as.numeric(as.character(mummichog.pathway$p.value)) <= 0.05),])
mummichog.metabID <- {}
mummichog.metabID <- strsplit(as.character(mummichog.sig.pathway$overlap_features..id.),";")
mummichog.metabmz <- {}
mummichog.metabmz <- strsplit(as.character(mummichog.sig.pathway$used_input_mzs),";")

mummichog.sig.pathway$pathway <- as.character(mummichog.sig.pathway$pathway)
mummichog.sig.pathway$pathway[mummichog.sig.pathway$pathway=="Urea cycle/amino group metabolism"] <- "Urea cycle metabolism" # Change some pathway name since folder name cannot have /

# Begin box plot
for(i in 1:nrow(mummichog.sig.pathway)){
  pathway.name <- as.character(mummichog.sig.pathway[i,1])
  print(pathway.name)     # Obtain the name of pathways
  
  setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/Pathway Box Plots")
  pdf_file<-paste(pathway.name,".pdf",sep="")
  pdf(file = pdf_file,width=10,height=10)
  
  for(j in 1:length(mummichog.metabID[[i]])){
    par(mfrow=c(1,3))
    metab.mz <- mummichog.metabmz[[i]][j]
    metab.mz <- as.vector(strsplit(as.character(metab.mz),","))[[1]]
    metab.mz <- as.numeric(metab.mz)
    metab.name <- mummichog.MatchMetab[which(as.character(mummichog.MatchMetab$id)==mummichog.metabID[[i]][j]),]
    metab.name <- as.character(metab.name[1,5])
    if(length(metab.name)==0){
      next}
    # metab.time <- a_KEGG[which(a_KEGG$KEGGID==metab.id),2]
    for(l in 1:length(metab.mz)){
      matchmz <- metab.mz[l]
      metab.expr <- save.plsresults.sigfeatures[which(save.plsresults.sigfeatures$mz==matchmz),-c(3:6)]    # link mz with residual data, get the residual record of that metablite
      if(nrow(metab.expr)==0){
        metab.expr <- wide_save_residual[which(wide_save_residual$mz==matchmz),]
      }
      for(k in 1:nrow(metab.expr)){                                           # For each metablite, it can be linked with multiple records
        box.mz <- metab.expr[k,1]
        box.time <- metab.expr[k,2]
        
        box.expr <- t(metab.expr[k,-c(1,2)])
        box.expr <- as.data.frame(cbind(box.expr,sampleID$factorcase))
        colnames(box.expr) <- c("Intensity","Group")
        box.expr$Intensity <- as.numeric(as.character(box.expr$Intensity))
        box.expr$Group <- as.factor(box.expr[,2])
        # title <- paste("metab: ",metab.name,"mz:",box.mz,"time:",box.time,sep = " ")
        title <- metab.name
        boxplot(Intensity~Group,data = box.expr,outline = FALSE,main=title)
        
        # Draw second comp
        comp3.metab.expr <- comp3.wide_save_residual[which(comp3.wide_save_residual$mz==box.mz & comp3.wide_save_residual$time==box.time),]
        comp3.box.expr <- t(comp3.metab.expr[1,-c(1,2)])
        comp3.box.expr <- as.data.frame(cbind(comp3.box.expr,comp3.sampleID$factorcase))
        colnames(comp3.box.expr) <- c("Intensity","Group")
        comp3.box.expr$Intensity <- as.numeric(as.character(comp3.box.expr$Intensity))
        comp3.box.expr$Group <- as.factor(comp3.box.expr[,2])
        if(nrow(comp3.metab.expr)==0){
          comp3.box.expr$Intensity=0
        }
        boxplot(Intensity~Group,data = comp3.box.expr,outline = FALSE)
        
        # Draw third comp
        comp4.metab.expr <- comp4.wide_save_residual[which(comp4.wide_save_residual$mz==box.mz & comp4.wide_save_residual$time==box.time),]
        comp4.box.expr <- t(comp4.metab.expr[1,-c(1,2)])
        comp4.box.expr <- as.data.frame(cbind(comp4.box.expr,comp4.sampleID$factorcase))
        colnames(comp4.box.expr) <- c("Intensity","Group")
        comp4.box.expr$Intensity <- as.numeric(as.character(comp4.box.expr$Intensity))
        comp4.box.expr$Group <- as.factor(comp4.box.expr[,2])
        if(nrow(comp4.metab.expr)==0){
          comp4.box.expr$Intensity=0
        }
        boxplot(Intensity~Group,data = comp4.box.expr,outline = FALSE)
      }
    }
  }
  dev.off()
}

###################################################
### VI. Create box plots II
###################################################
library(xlsx)

load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_all_residual.RData")
load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/Res_PLSDA_result_2018-04-30_vip2fc0.58.RData")
wide_save_residual$mz <- round(wide_save_residual$mz,4)
wide_save_residual$time <- round(wide_save_residual$time,2)
save.plsresults.sigfeatures$mz <- round(save.plsresults.sigfeatures$mz,4)
save.plsresults.sigfeatures$time <- round(save.plsresults.sigfeatures$time,2)

sampleID <- class
rm(class)
sampleID$factorcase <- paste(sampleID$casecontrol..Factor.1.,sampleID$Exposure_category..Factor2.,sep = "")

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichogMTBNK_vip2fc0.58/tsv")
mummichog.pathway <- read.xlsx("mcg_pathwayanalysis_C18_mummichogMTBNK_vip2fc0.58.xlsx",sheetName = "Sheet1")
mummichog.MatchMetab <- read.xlsx("_tentative_featurematch_C18_mummichogMTBNK_vip2fc0.58.xlsx",sheetName = "_tentative_featurematch_")

mummichog.sig.pathway <-as.data.frame(mummichog.pathway[which(as.numeric(as.character(mummichog.pathway$p.value)) <= 0.06),]) ## not enough pathways, so losen the p value
mummichog.metabID <- {}
mummichog.metabID <- strsplit(as.character(mummichog.sig.pathway$overlap_features..id.),";")
mummichog.metabmz <- {}
mummichog.metabmz <- strsplit(as.character(mummichog.sig.pathway$used_input_mzs),";")

mummichog.sig.pathway$pathway <- as.character(mummichog.sig.pathway$pathway)
mummichog.sig.pathway$pathway[mummichog.sig.pathway$pathway=="Urea cycle/amino group metabolism"] <- "Urea cycle metabolism" # Change some pathway name since folder name cannot have /

# Begin box plot
for(i in 1:nrow(mummichog.sig.pathway)){
  pathway.name <- as.character(mummichog.sig.pathway[i,1])
  print(pathway.name)     # Obtain the name of pathways
  
  setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/Pathway Box Plots")
  pdf_file<-paste(pathway.name,".pdf",sep="")
  pdf(file = pdf_file,width=10,height=10)
  
  for(j in 1:length(mummichog.metabID[[i]])){
    metab.mz <- mummichog.metabmz[[i]][j]
    metab.mz <- as.vector(strsplit(as.character(metab.mz),","))[[1]]
    metab.mz <- as.numeric(metab.mz)
    metab.name <- mummichog.MatchMetab[which(as.character(mummichog.MatchMetab$id)==mummichog.metabID[[i]][j]),]
    metab.name <- as.character(metab.name[1,5])
    if(length(metab.name)==0){
      next}
    for(l in 1:length(metab.mz)){
      matchmz <- metab.mz[l]
      metab.expr <- save.plsresults.sigfeatures[which(save.plsresults.sigfeatures$mz==matchmz),-c(3:6)]    # link mz with residual data, get the residual record of that metablite
      if(nrow(metab.expr)==0){
        metab.expr <- wide_save_residual[which(wide_save_residual$mz==matchmz),]
      }
      for(k in 1:nrow(metab.expr)){                                           # For each metablite, it can be linked with multiple records
        box.mz <- metab.expr[k,1]
        box.time <- metab.expr[k,2]
        box.expr <- wide_save_residual[which(wide_save_residual$mz==box.mz & wide_save_residual$time==box.time),]
        if(nrow(box.expr)==0){
          next
        }                                                 # Some metabolites can be found in saave.sigfeatures but not wide.save.residula for all subjects, that's because it might be missing > 80% in other groups and discard during preprocess
        box.expr <- t(box.expr[1,-c(1,2)])
        box.expr <- as.data.frame(cbind(box.expr,sampleID$factorcase))
        colnames(box.expr) <- c("Intensity","Group")
        box.expr$Intensity <- as.numeric(as.character(box.expr$Intensity))
        box.expr$Group <- as.factor(box.expr[,2])
        title <- paste(metab.name,"\nmz:",box.mz,sep = " ")
        # title <- metab.name
        boxplot(Intensity~Group,data = box.expr,outline = FALSE,main=title)
      }
    }
  }
  dev.off()
}

#################
# Plot modules
#################
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichogMTBNK_vip2fc0.58/tsv")
mummichog.sig.module <- read.xlsx("mcg_modularanalysis_C18_mummichogMTBNK_vip2fc0.58.xlsx",sheetName = "Sheet1")
mummichog.MatchMetab <- read.xlsx("_tentative_featurematch_C18_mummichogMTBNK_vip2fc0.58.xlsx",sheetName = "_tentative_featurematch_")

mummichog.metabID <- {}
mummichog.metabID <- strsplit(as.character(mummichog.sig.module$members..id.),",")

mummichog.sig.module$MODULE <- as.character(mummichog.sig.module$MODULE)

# Begin box plot
for(i in 1:nrow(mummichog.sig.module)){
  module.name <- as.character(mummichog.sig.module[i,1])
  print(module.name)     # Obtain the name of modules
  
  setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/pathway Box Plots")
  pdf_file<-paste(module.name,".pdf",sep="")
  pdf(file = pdf_file,width=10,height=10)
  
  for(j in 1:length(mummichog.metabID[[i]])){
    metab.id <- mummichog.metabID[[i]][j]
    metab.mz <- mummichog.MatchMetab[which(as.character(mummichog.MatchMetab$id)==metab.id),1]
    # metab.mz <- as.vector(strsplit(as.character(metab.mz),","))[[1]]
    # metab.mz <- as.numeric(metab.mz)
    metab.name <- mummichog.MatchMetab[which(as.character(mummichog.MatchMetab$id)==metab.id),]
    metab.name <- as.character(metab.name[1,5])
    if(length(metab.name)==0){
      next}
    for(l in 1:length(metab.mz)){
      matchmz <- metab.mz[l]
      metab.expr <- save.plsresults.sigfeatures[which(save.plsresults.sigfeatures$mz==matchmz),-c(3:6)]    # link mz with residual data, get the residual record of that metablite
      if(nrow(metab.expr)==0){
        metab.expr <- wide_save_residual[which(wide_save_residual$mz==matchmz),]
      }
      for(k in 1:nrow(metab.expr)){                                           # For each metablite, it can be linked with multiple records
        box.mz <- metab.expr[k,1]
        box.time <- metab.expr[k,2]
        box.expr <- wide_save_residual[which(wide_save_residual$mz==box.mz & wide_save_residual$time==box.time),]
        if(nrow(box.expr)==0){
          next
        }                                                 # Some metabolites can be found in save.sigfeatures but not wide.save.residula for all subjects, that's because it might be missing > 80% in other groups and discard during preprocess
        box.expr <- t(box.expr[1,-c(1,2)])
        box.expr <- as.data.frame(cbind(box.expr,sampleID$factorcase))
        colnames(box.expr) <- c("Intensity","Group")
        box.expr$Intensity <- as.numeric(as.character(box.expr$Intensity))
        box.expr$Group <- as.factor(box.expr[,2])
        title <- paste(metab.name,"\nmz:",box.mz,sep = " ")
        # title <- metab.name
        boxplot(Intensity~Group,data = box.expr,outline = FALSE,main=title)
      }
    }
  }
  dev.off()
}


###################################################
### VII. Metaboanalyst
###################################################

library(xlsx)

load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input/C18_control_expo_unexpo_classification_nonorm.RData")
load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_input/C18_control_expo_unexpo_residual_nonorm_WGCNA.RData")

row.names(wide_save_residual) <- paste("met_",1:nrow(wide_save_residual),sep = "")
feature <- cbind(sampleID[,c(1:2)],t(wide_save_residual[,-c(1:2)]))
write.table(feature,file="C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/metaboanalyst.txt",sep = "\t",row.names = F,quote = F)

metaboanalyst.sig <- read.csv("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/imp_features_cv.csv")
after.prepro.linkid$X <- row.names(after.prepro.linkid)
metaboanalyst.sig_cov <- merge(after.prepro.linkid,metaboanalyst.sig, by = "X")
metaboanalyst.sig_cov <- metaboanalyst.sig_cov[order(-metaboanalyst.sig_cov$Importance),]

load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/Res_PLSDA_result_2018-01-23.RData")
temp <- save.plsresults.sigfeatures[,c(1,6)]
check.cons <- merge(temp,metaboanalyst.sig_cov,by = "mz")
check.cons <- check.cons[order(check.cons$rank),]

###################################################
### VIII. Activity network for cytoscape/KEGG mapper
###################################################

# Network

load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/Res_PLSDA_result_2018-04-30_vip2fc0.58.RData")
save.plsresults.allfeatures$mz <- round(save.plsresults.allfeatures$mz,4)
save.plsresults.sigfeatures$mz <- round(save.plsresults.sigfeatures$mz,4)

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichogMTBNK_vip2fc0.58/tsv")
mummichog.NetAct <- read.xlsx("InspectedNodes_ActivityNetwork.xlsx",sheetName = "Sheet1")
mummichog.NetAct$number <- c(1:nrow(mummichog.NetAct))
mummichog.NetAct <- merge(mummichog.NetAct,save.plsresults.allfeatures,by = "mz", all.x = TRUE)
mummichog.NetAct <- mummichog.NetAct[order(mummichog.NetAct$number),1:6]

mummichog.NetAct <- merge(mummichog.NetAct,save.plsresults.sigfeatures,by = "mz", all.x = TRUE)
mummichog.NetAct <- mummichog.NetAct[order(mummichog.NetAct$number),c(1:6,11)]

mummichog.cor <- read.table("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/mummichogMTBNK_input_vip2fc0.58.txt",sep = "\t",header = TRUE)
mummichog.cor <- mummichog.cor[which(mummichog.cor$p.value == 0.04),]
mummichog.cor$mz <- round(mummichog.cor$mz,4)
mummichog.NetAct <- merge(mummichog.NetAct,mummichog.cor,by = "mz", all.x = TRUE)
mummichog.NetAct <- mummichog.NetAct[order(mummichog.NetAct$number),c(1:7,9)]

# Pathways

load("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/PANDA_output_PLSDA/Res_PLSDA_result_2018-05-31_vip2fc0.RData")
save.plsresults.allfeatures$mz <- round(save.plsresults.allfeatures$mz,4)
save.plsresults.sigfeatures$mz <- round(save.plsresults.sigfeatures$mz,4)

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichog_vip2fc0_0531/tsv")
mummichog.PatAct <- read.xlsx("mcg_pathwayanalysis_C18_mummichog_vip2fc0_0531.xlsx",sheetName = "Sheet1",header = TRUE)
# mummichog.PatAct <- mummichog.PatAct[,c(5,1,2)]
# colnames(mummichog.PatAct) <- c("name","mz","id")
colnames(mummichog.PatAct) <- c("mz","id","match.form","mz_difference","name","pathway")
mummichog.PatAct$number <- c(1:nrow(mummichog.PatAct))
mummichog.PatAct <- merge(mummichog.PatAct,save.plsresults.allfeatures,by = "mz", all.x = TRUE)
mummichog.PatAct <- mummichog.PatAct[order(mummichog.PatAct$number),1:10]

mummichog.PatAct <- merge(mummichog.PatAct,save.plsresults.sigfeatures,by = "mz", all.x = TRUE)
mummichog.PatAct <- mummichog.PatAct[order(mummichog.PatAct$number),c(1:7,15)]

mummichog.cor <- read.table("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/mummichog_input_vip2fc0_2018-05-31.txt",sep = "\t",header = TRUE)
mummichog.cor <- mummichog.cor[which(mummichog.cor$p.value == 0.004),]
mummichog.cor$mz <- round(mummichog.cor$mz,4)
mummichog.PatAct <- merge(mummichog.PatAct,mummichog.cor,by = "mz", all.x = TRUE)
mummichog.PatAct <- mummichog.PatAct[order(mummichog.PatAct$number),c(1:8,11)]

mummichog.PatAct <- mummichog.PatAct[-which(is.na(mummichog.PatAct$rank)),]
mummichog.PatAct <- mummichog.PatAct[-which(duplicated(mummichog.PatAct$id)),]
mummichog.PatAct$color[mummichog.PatAct$foldchange > 0] <- "red"
mummichog.PatAct$color[mummichog.PatAct$foldchange < 0] <- "green"
for.kegg <- mummichog.PatAct[,c(2,10)]

# add HILIC

load("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA/Res_PLSDA_result_2018-05-31_vip2fc0.RData")
save.plsresults.allfeatures$mz <- round(save.plsresults.allfeatures$mz,4)
save.plsresults.sigfeatures$mz <- round(save.plsresults.sigfeatures$mz,4)

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/HILIC_mummichog_vip2fc0_0531/tsv")
HILIC.mummichog.PatAct <- read.xlsx("mcg_pathwayanalysis_HILIC_mummichog_vip2fc0_0531.xlsx",sheetName = "Sheet1",header = TRUE)
# HILIC.mummichog.PatAct <- HILIC.mummichog.PatAct[,c(5,1,2)]
# colnames(HILIC.mummichog.PatAct) <- c("name","mz","id")
colnames(HILIC.mummichog.PatAct) <- c("mz","id","match.form","mz_difference","name","pathway")
HILIC.mummichog.PatAct$number <- c(1:nrow(HILIC.mummichog.PatAct))
HILIC.mummichog.PatAct <- merge(HILIC.mummichog.PatAct,save.plsresults.allfeatures,by = "mz", all.x = TRUE)
HILIC.mummichog.PatAct <- HILIC.mummichog.PatAct[order(HILIC.mummichog.PatAct$number),1:10]

HILIC.mummichog.PatAct <- merge(HILIC.mummichog.PatAct,save.plsresults.sigfeatures,by = "mz", all.x = TRUE)
HILIC.mummichog.PatAct <- HILIC.mummichog.PatAct[order(HILIC.mummichog.PatAct$number),c(1:7,15)]

HILIC.mummichog.cor <- read.table("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_mummichog/mummichog_input_vip2fc0_2018-05-31.txt",sep = "\t",header = TRUE)
HILIC.mummichog.cor <- HILIC.mummichog.cor[which(HILIC.mummichog.cor$p.value == 0.004),]
HILIC.mummichog.cor$mz <- round(HILIC.mummichog.cor$mz,4)
HILIC.mummichog.PatAct <- merge(HILIC.mummichog.PatAct,HILIC.mummichog.cor,by = "mz", all.x = TRUE)
HILIC.mummichog.PatAct <- HILIC.mummichog.PatAct[order(HILIC.mummichog.PatAct$number),c(1:8,11)]

HILIC.mummichog.PatAct <- HILIC.mummichog.PatAct[-which(is.na(HILIC.mummichog.PatAct$rank)),]
HILIC.mummichog.PatAct <- HILIC.mummichog.PatAct[-which(duplicated(HILIC.mummichog.PatAct$id)),]
HILIC.mummichog.PatAct$color[HILIC.mummichog.PatAct$foldchange > 0] <- "red"
HILIC.mummichog.PatAct$color[HILIC.mummichog.PatAct$foldchange < 0] <- "green"

# Combined HILIC and C18 pathways together
mummichog.PatAct <- rbind(mummichog.PatAct, HILIC.mummichog.PatAct)

# Draw pathway plots
source("C:/Users/QiYan/Dropbox/AIME/Archive/pathway_Plotly.R")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Fatty acid activation")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "De novo fatty acid biosynthesis")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Fatty Acid Metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Butanoate metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Linoleate metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Lysine metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Methionine and cysteine metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Glycolysis and Gluconeogenesis")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Prostaglandin formation from arachidonate")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Vitamin E metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Glycerophospholipid metabolism")

pathway_Plotly(data = mummichog.PatAct, pathway_name = "Urea cycle/amino group metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Glycosphingolipid metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Pyrimidine metabolism")
pathway_Plotly(data = mummichog.PatAct, pathway_name = "Leukotriene metabolism")

# Use pathview to draw KEGG pathways
cmp.list <- mummichog.PatAct$foldchange
names(cmp.list) <- mummichog.PatAct$id

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_C18neg/C18_Controls_ExpoUnexpo/C18_mummichog/C18_mummichog_vip2fc0_0531")
pv.out <- pathview(cpd.data = cmp.list,
                   pathway.id = '01100', species = "hsa", keys.align = "y", kegg.native = TRUE, res = 1000, bins = 80,
                   plot.col.key = FALSE, new.signature = FALSE,
                   low = list(gene = "blue", cpd = "green"), high = list(gene = "yellow", cpd ="red"))
