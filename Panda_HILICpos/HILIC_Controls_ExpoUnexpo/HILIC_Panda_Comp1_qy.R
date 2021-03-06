library(tidyverse)
library(xmsPANDA)
# library(versions)
# install.versions('mixOmics','5.2.0')
library(mixOmics)

#xmsPANDA

class <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_residuals_classlabels_control_expo_unexpo.txt"
# feature <- "C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt"
# residual <- "C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_residuals_control_expo_unexpo.txt"
# 
# feature <- read.table(feature,sep="\t",header=TRUE)
# # class <- read.table(class,sep="\t",header=TRUE)
# residual <- read.table(residual,sep="\t",header=TRUE)

# PLS-DA_residual
# residual: All
residual <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_residuals_control_expo_unexpo.txt"
outloc <-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA_residual/residual"

ptm <- proc.time()

demetabs_res<-diffexp(feature_table_file=residual,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 1,
                      feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="median",summary.na.replacement="none",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE,
                      quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(-9999999), pairedanalysis = FALSE, featselmethod="pls",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor=TRUE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=NA,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2,
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=-9999999,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=NA,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=30)

sink(file=NULL)
proc.time() - ptm
# residual_1: sex
residual_1 <- "C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_residuals_1_control_expo_unexpo.txt"
outloc <-"C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA_residual/residual_1"

ptm <- proc.time()

demetabs_res<-diffexp(feature_table_file=residual_1,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 1,
                      feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="median",summary.na.replacement="halfdatamin",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE,
                      quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="pls",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor  =TRUE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=NA,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2,
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=30)

sink(file=NULL)
proc.time() - ptm
# residual_2: sex, maternal edu
residual_2 <- "C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_residuals_2_control_expo_unexpo.txt"
outloc <-"C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA_residual/residual_2"

ptm <- proc.time()

demetabs_res<-diffexp(feature_table_file=residual_2,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 1,
                      feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="median",summary.na.replacement="halfdatamin",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = FALSE, medcenter=FALSE, znormtransform = TRUE,
                      quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="pls",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor  =TRUE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2,
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=30)

sink(file=NULL)
proc.time() - ptm
# residual_3: sex, maternal edu, pregcompl, ttcbl, maternal raceeth
residual_3 <- "C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_residuals_3_control_expo_unexpo.txt"
outloc <-"C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA_residual/residual_3"

ptm <- proc.time()

demetabs_res<-diffexp(feature_table_file=residual_3,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 1,
                      feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="median",summary.na.replacement="halfdatamin",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = FALSE, medcenter=FALSE, znormtransform = TRUE,
                      quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="pls",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor  =TRUE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2,
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=30)

sink(file=NULL)
proc.time() - ptm

#SVM
outloc <-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_SVM"

ptm <- proc.time()
demetabs_res<-diffexp(feature_table_file=feature,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 3,
                      feat.filt.thresh =NA, summarize.replicates =TRUE, summary.method="median",summary.na.replacement="halfdatamin",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE,
                      quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="rfesvm",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor=TRUE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 6,
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=5)

sink(file=NULL)
proc.time() - ptm

#SVM_residual
outloc <-"C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_SVM_residual"

ptm <- proc.time()
demetabs_res<-diffexp(feature_table_file=residual,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 1,
                      feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="median",summary.na.replacement="halfdatamin",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = FALSE, medcenter=FALSE, znormtransform = TRUE,
                      quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="rfesvm",
                      fdrthresh = 0.05, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor=TRUE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 6,
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=5)

sink(file=NULL)
proc.time() - ptm

#LIMMA
# outloc <-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Exposed_CasesControls/PANDA_output_LIMMA"
# 
# ptm <- proc.time()
# demetabs_res<-diffexp(feature_table_file=feature,
#                       parentoutput_dir=outloc,
#                       class_labels_file=class,
#                       num_replicates = 3,
#                       feat.filt.thresh =NA, summarize.replicates =TRUE, summary.method="median",summary.na.replacement="halfdatamin",
#                       rep.max.missing.thresh=0.5,
#                       all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
#                       log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE,
#                       quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
#                       rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="limma",
#                       fdrthresh = 0.05, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
#                       kfold=10,feat_weight=1,globalcor=TRUE,target.metab.file=NA,rf_selmethod = "rawVIMsig",
#                       target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
#                       samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
#                       net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2,
#                       max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
#                       rocfeatincrement=TRUE,
#                       rocclassifier="svm",
#                       foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
#                       optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
#                       pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=30)
# 
# sink(file=NULL)
# proc.time() - ptm

#RF
outloc <-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_RF_residual"

ptm <- proc.time()
demetabs_res<-diffexp(feature_table_file=feature,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 3,
                      feat.filt.thresh =NA, summarize.replicates =TRUE, summary.method="median",summary.na.replacement="halfdatamin",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE,
                      quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="RF",
                      fdrthresh = 0.2, fdrmethod="none",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor=TRUE,target.metab.file=NA,rf_selmethod = "rankbased",
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=2000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 2,
                      max_varsel = 150, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",
                      foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=30)

sink(file=NULL)
proc.time() - ptm

#lmreg
class <- "C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_withcov_classlabels_control_expo_unexpo.txt"
feature <- "C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input/HILIC_ftrsmzcalib_combat_ordered_control_expo_unexpo.txt"
outloc <-"C:/Users/Qi/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_LogitReg"

ptm <- proc.time()
demetabs_res<-diffexp(feature_table_file=feature,
                      parentoutput_dir=outloc,
                      class_labels_file=class,
                      num_replicates = 3,
                      feat.filt.thresh =NA, summarize.replicates =TRUE, summary.method="median",summary.na.replacement="halfdatamin",
                      rep.max.missing.thresh=0.5,
                      all.missing.thresh=0.5, group.missing.thresh=0.8, input.intensity.scale="raw",
                      log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE,
                      quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
                      rsd.filt.list = c(1), pairedanalysis = FALSE, featselmethod="lmreg",
                      fdrthresh = 0.3, fdrmethod="BH",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
                      kfold=10,feat_weight=1,globalcor=TRUE,target.metab.file=NA,
                      target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,missing.val=0,networktype="complete",
                      samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
                      net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3, pls_vip_thresh = 2, num_nodes = 6,
                      max_varsel = 100, pls_ncomp = 5,pcacenter=TRUE,pcascale=TRUE,pred.eval.method="CV",rocfeatlist=seq(2,10,1),
                      rocfeatincrement=TRUE,
                      rocclassifier="svm",foldchangethresh=0.58,wgcnarsdthresh=30,WGCNAmodules=FALSE,
                      optselect=TRUE,max_comp_sel=2,saveRda=FALSE,pca.cex.val=4,pls.permut.count=100,
                      pca.ellipse=FALSE,ellipse.conf.level=0.95,legendlocation="bottomleft",svm.acc.tolerance=5)

sink(file=NULL)
proc.time() - ptm
