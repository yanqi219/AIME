library(ggplot2)
# library(xmsPANDA)

#Load data

load("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_input/HILIC_case_control_noexposure_classification_nonorm.RData")
load("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_output_PLSDA/Res_PLSDA_result_2018-08-16_vip2fc0.RData")

# comp.mode = "Unexpo_CaseControl" # "Control_ExpoUnexpo","Unexpo_CaseControl","Case_ExpoUnexpo","Expo_CaseControl"
# 
# if(comp.mode == "Control_ExpoUnexpo"){
#   
#   dir.folder <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/"
#   filename.ins <- "_control_expo_unexpo"
#   
# }else if(comp.mode == "Unexpo_CaseControl"){
#   
#   dir.folder <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/"
#   filename.ins <- "_case_control_noexposure"
#   
# }else if(comp.mode =="Case_ExpoUnexpo"){
#   
#   dir.folder <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Cases_ExpoUnexpo/"
#   filename.ins <- "_case_expo_unexpo"
#   
# }else if(comp.mode =="Expo_CaseControl"){
#   
#   dir.folder <- "C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Exposed_CasesControls/"
#   filename.ins <- "_case_control_exposure"
#   
# }
# 
# # TO plot box plot, we need raw intensity data (after filtering)
# 
# class <- paste(dir.folder,paste("PANDA_input/HILIC_classlabels",filename.ins,".txt",sep = ""),sep = "")
# feature <- paste(dir.folder,paste("PANDA_input/HILIC_ftrsmzcalib_combat_ordered",filename.ins,".txt",sep = ""),sep = "")
# outloc <- paste(dir.folder,"PANDA_output_PLSDA",sep = "")
# 
# ready_for_regression<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file=feature,parentoutput_dir=outloc,class_labels_file=class,num_replicates=3,feat.filt.thresh=NA,
#                                       summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.5,group.missing.thresh=0.8,
#                                       log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,
#                                       samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="halfdatamin",featselmethod=NA)
# 
# feature <- as.data.frame(ready_for_regression$data_matrix_afternorm_scaling)
# na_count <-sapply(feature, function(y) sum(is.na(y)))
# summary(na_count)
# 
# # remove na
# col.name <- c("mz","time",sampleID$SampleID)
# feature <- feature[,colnames(feature) %in% col.name]
# feature.temp <- feature[,c(1:2)]
# feature.order <- as.data.frame(t(feature[,-c(1:2)]))
# feature.order <- feature.order[order(row.names(feature.order)),]
# feature.order <- as.data.frame(t(feature.order))
# feature <- cbind(feature.temp,feature.order)
###############################################################
# feature <- as.data.frame(t(after.prepro.feature))
# feature <- cbind(after.prepro.linkid, feature)

load("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Non_Exposed_CasesControls/PANDA_input/HILIC_case_control_noexposure_raw.RData")

feature$mz <- round(feature$mz,4)
feature$time <- round(feature$time,1)

feature_w_fc <- cbind(feature[,c(1,2)], save.plsresults.allfeatures[,c(3,4)], feature[,c(3:ncol(feature))])

targeted_box_plot <-
  function(Expr=feature_w_fc,case=sampleID,mz,time,metablite_name,adductform,pathway){
    box.mz <- mz
    box.time <- time
    
    metab.expr <- Expr[which(Expr$mz==box.mz & Expr$time==box.time),]
    
    box.vip <- round(metab.expr$vip,2)
    box.foldchange <- round(metab.expr$foldchange,2)
    
    box.expr <- t(metab.expr[1,-c(1:4)])
    box.expr <- as.data.frame(cbind(box.expr,sampleID$factorcase))
    colnames(box.expr) <- c("Intensity","Group")
    box.expr$Intensity <- as.numeric(as.character(box.expr$Intensity))
    box.expr$Group <- as.factor(box.expr[,2])
    title <- paste(metablite_name,"\nM/Z:",box.mz,"Time:",box.time,"\nLog(fc):",box.foldchange,"Adduct:",adductform,"Pathway:",pathway,sep = " ")
    # title <- metablite_name
    # ggplot(data = box.expr, aes(x = Group,y = Intensity)) + 
    #   stat_boxplot(geom ='errorbar') + 
    #   geom_boxplot(outlier.shape = NA)
    boxplot(Intensity~Group,data = box.expr,outline = FALSE,main=title, range = 1.5, width = c(1,1), ylab = "Log(Intensity)")
  }

# targeted_box_plot(mz = 608.08805,time = 25.81, metablite_name = "UDP-N-acetyl-D-galactosamine", adductform = "M+H")
# targeted_box_plot(mz = 113.03552,time = 217.62, metablite_name = "Uracil", adductform = "M+H")
# targeted_box_plot(mz = 155.04269,time = 98.32, metablite_name = "L-Asparagine", adductform = "M+H")
# targeted_box_plot(mz = 197.10093,time = 103.31, metablite_name = "L-Arginine", adductform = "M+H")
# targeted_box_plot(mz = 974.28712,time = 57.75, metablite_name = "(2R,6R,10R)-trimethyl-hendecanoyl-CoA", adductform = "M+H")
# targeted_box_plot(mz = 281.24742,time = 33.64, metablite_name = "Linoleate", adductform = "M+H")
# targeted_box_plot(mz = 169.03562,time = 100.23, metablite_name = "Urate", adductform = "M+H")

#HILICpos
targeted_box_plot(mz = 281.2474,time = 33.6, metablite_name = "Linoleate", adductform = "M+H[1+]", pathway = "Fatty acid activation")
targeted_box_plot(mz = 251.2006,time = 27.4, metablite_name = "Tetradecanoic acid", adductform = "M+Na[1+]", pathway = "Fatty acid activation")
targeted_box_plot(mz = 297.2422,time = 239.1, metablite_name = "13(S)-HODE", adductform = "M+H[1+]", pathway = "Linoleate metabolism")
targeted_box_plot(mz = 279.2318,time = 33.6, metablite_name = "Gamma-Linolenic acid", adductform = "M+H[1+]", pathway = "Linoleate metabolism")
targeted_box_plot(mz = 385.3460,time = 33.4, metablite_name = "Cholestenone", adductform = "M+H[1+]", pathway = "sterol lipid")
targeted_box_plot(mz = 585.2699,time = 29.1, metablite_name = "Bilirubin", adductform = "M+H[1+]", pathway = "Bile acid")
#HILICpos
targeted_box_plot(mz = 335.2227,time = 119.9, metablite_name = "LTB4", adductform = "M-H[1-]", pathway = "Leukotriene metabolism")
targeted_box_plot(mz = 327.2330,time = 234.4, metablite_name = "Docosahexaenoic acid", adductform = "M-H[1-]", pathway = "W3 fatty acid")
targeted_box_plot(mz = 317.2119,time = 179.8, metablite_name = "Prostaglandin B1", adductform = "M-H2O-H[1-]", pathway = "Arachidonic acid metabolism")
targeted_box_plot(mz = 359.2227,time = 114.5, metablite_name = "Progesterone", adductform = "M+HCOO[1-]", pathway = "Steroid hormones")
targeted_box_plot(mz = 361.2382,time = 146.5, metablite_name = "5beta-Pregnane-3,20-dione", adductform = "M+HCOO[1-]", pathway = "Steroid hormones")



# Draw bar plots

pathway <- read.table(file = "C:/Users/QiYan/Dropbox/AIME/Results/Air pollution/Enrichment pathways.txt", header = T, sep = "\t")
pathway <- pathway[-15,]
pathway$new.name <- paste(pathway$Pathways,"(",pathway$overlap.size,"/",pathway$pathway.size,")",sep = "")

p <- ggplot(data=pathway, aes(x=reorder(new.name, -log(p.value)), y=-log(p.value))) +
  geom_bar(aes(fill=-log(p.value)),stat="identity", position=position_dodge()) +
  scale_fill_gradient(low = "#edf8b1", high = "#2c7fb8", space = "Lab",
                      na.value = "grey50", guide = FALSE) +
  geom_hline(yintercept=-log(0.05), linetype="dashed", color = "deeppink3", size = 1.2) +
  # theme_minimal() + 
  theme(axis.text.y=element_text(size=rel(3)), axis.title.y=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "-Log(P value)") +
  coord_flip()
p

# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")
