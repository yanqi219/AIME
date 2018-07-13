###################################################
### VIII. Re-annotate significant features with loose adducts
###################################################

# Install packages

# source("http://bioconductor.org/biocLite.R")
# 
# install.packages("snow",repos="http://cran.r-project.org")
# install.packages("doSNOW",repos="http://cran.r-project.org")
# install.packages("parallel",repos="http://cran.r-project.org")
# install.packages("e1071",repos="http://cran.r-project.org")
# install.packages("XML",repos="http://cran.r-project.org")
# install.packages("R2HTML",repos="http://cran.r-project.org")
# install.packages("RCurl",repos="http://cran.r-project.org")
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rdisop",suppressUpdates=TRUE)
# biocLite("KEGGREST",suppressUpdates=TRUE)
# biocLite("pcaMethods",suppressUpdates=TRUE)
# install.packages("flashClust",repos="http://cran.r-project.org")
# install.packages("plyr",repos="http://cran.r-project.org")
# install.packages("png",repos="http://cran.r-project.org")
# install.packages("rjson",repos="http://cran.r-project.org")
# 
# install.packages("/u/home/q/qyan/AIME/r_package/SSOAP_0.9-0.tar.gz", repos = NULL, type = "source")
# install.packages("/u/home/q/qyan/AIME/r_package/XMLSchema_0.7-0.tar.gz", repos = NULL, type = "source")
# install.packages("/u/home/q/qyan/AIME/r_package/xMSannotator_1.3.2.tar.gz", repos = NULL, type = "source")

library(xMSannotator)

# setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_output_PLSDA")
# load(file = "Res_PLSDA_result_2018-05-31_vip2fc0.RData")
setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")
# setwd("/u/home/q/qyan/AIME/Annotation")
load(file = "HILIC_control_expo_unexpo_residual_nonorm_WGCNA.RData")

data(adduct_table)
data(adduct_weights)

normalization <- function(x){
  return((x-mean(x))/(max(x)-min(x)))
}

###########Parameters to change##############

# input_data <- save.plsresults.sigfeatures[,-c(3:6)]
# input_data.norm <- cbind(input_data[,c(1,2)],apply(input_data[,-c(1,2)],2,normalization))
input_data <- wide_save_residual
dataA<-input_data

outloc<-"C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/HILIC_Annotation/sigfeature_annotation"
# outloc<-"/u/home/q/qyan/AIME/Annotation/HILIC"

max.mz.diff<-10  #mass search tolerance for DB matching in ppm
max.rt.diff<-10 #retention time tolerance between adducts/isotopes
corthresh<-0.5 #correlation threshold between adducts/isotopes
max_isp=5 #maximum number of isotopes to search for
mass_defect_window=0.01 #mass defect window for isotope search

num_nodes<-4   #number of cores to be used; 2 is recommended for desktop computers due to high memory consumption

db_name="HMDB" #other options: HMDB,Custom,KEGG, LipidMaps, T3DB
status="Detected and Quantified" #other options: "Detected", NA, "Detected and Quantified", "Expected and Not Quantified"
num_sets<-300 #number of sets into which the total number of database entries should be split into;

mode<-"pos" #ionization mode
# queryadductlist=c("M+H","M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O") #other options: c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H"); c("positive"); c("negative"); c("all");see data(adduct_table) for complete list
queryadductlist=c("all")

#provide list of database IDs (depending upon selected database) for annotating only specific metabolites
customIDs<-NA #c("HMDB15352","HMDB60549","HMDB00159","HMDB00222"); read.csv("/Users/mzmatch_95stdmx_HMDBIDs.csv")
customDB<-NA 
  
#########################

dataA<-unique(dataA)
print(dim(dataA))
print(format(Sys.time(), "%a %b %d %X %Y"))

system.time(annotres<-multilevelannotation(dataA=dataA,max.mz.diff=max.mz.diff,max.rt.diff=max.rt.diff,cormethod="spearman",
                                           num_nodes=num_nodes,queryadductlist=queryadductlist,
                                           mode=mode,outloc=outloc,db_name=db_name, adduct_weights=adduct_weights,num_sets=num_sets,allsteps=TRUE,
                                           corthresh=corthresh,NOPS_check=TRUE,customIDs=customIDs,missing.value=NA,
                                           deepsplit=4,networktype="unsigned",minclustsize=5,module.merge.dissimilarity=0.5,filter.by=c("M+H"),
                                           biofluid.location="Blood",origin=NA,status=status,boostIDs=NA,max_isp=max_isp,
                                           customDB=customDB,MplusH.abundance.ratio.check = FALSE, min_ions_perchem = 1,
                                           HMDBselect="all",mass_defect_window=mass_defect_window,pathwaycheckmode="pm",mass_defect_mode="both")
)


print(format(Sys.time(), "%a %b %d %X %Y"))

# pkg <- "package:xMSannotator"
# detach(pkg, character.only = TRUE)