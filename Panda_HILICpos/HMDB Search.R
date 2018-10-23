library('rvest')
library(doSNOW)
library(foreach)

cl<-makeCluster(4)
registerDoSNOW(cl)

HILIClist <- read.csv("C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/IROA_input_file_hilicpos_vjune282018.txt", sep = "\t", header = T)
C18list <- read.csv("C:/Users/QiYan/Dropbox/AIME/PNS_Ritz/IROA_input_file_c18neg_vjune282018.txt", sep = "\t", header = T)

HILIC_HMDBid <- as.character(HILIClist$HMDBID)
C18_HMDBid <- as.character(C18list$HMDBID)

#HMDB
HILIC_name <- vector('character')
for(i in 1: length(HILIC_HMDBid)){

  scraping <- read_html(paste("http://www.hmdb.ca/metabolites/",HILIC_HMDBid[i],sep = ""))
  
  METname <- scraping %>%
    html_nodes("strong") %>%
    html_text()
  METname <- METname[1]
  HILIC_name[i] <- METname
}

#C18
C18_name <- foreach(i=1:length(C18_HMDBid),.packages='rvest',.combine = rbind) %dopar% {
  scraping <- read_html(paste("http://www.hmdb.ca/metabolites/",C18_HMDBid[i],sep = ""))
  METname <- scraping %>%
    html_nodes("strong") %>%
    html_text()
  METname <- METname[1]
  C18_name[i] <- METname
}
stopCluster(cl)

HILIClist <- cbind(HILIClist,HILIC_name)
HILIClist$HILIC_name[which(HILIClist$HMDBID=="")] <- NA
C18list <- cbind(C18list,C18_name)
C18list$C18_name[which(C18list$HMDBID=="")] <- NA
