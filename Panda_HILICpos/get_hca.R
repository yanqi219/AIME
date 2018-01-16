get_hca <-
function(data_m,classlabels,cor.method="spearman",is.data.znorm=FALSE, clu.method="bicor",
plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3)
  {
    # data_m<-t(data_m)
    col_samples<-TRUE
    
    class_labels_levels<-levels(classlabels)
    ordered_labels<-classlabels
    
    class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="")
    
    col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
    
    sampleclass<-{}
    patientcolors<-{}
    classlabels<-as.data.frame(classlabels)
    f<-factor(classlabels[,1])
    for(c in 1:length(class_labels_levels)){
      num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
      sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
      patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group_cur))
    }

    heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    heatmap_cols<-rev(heatmap_cols)
    
    if(clu.method == "bicor"){
      hc <- try(hclust(as.dist(1-WGCNA::bicor(data_m,use="pairwise.complete.obs")))) #metabolites
      hr <- try(hclust(as.dist(1-WGCNA::bicor(t(data_m),use="pairwise.complete.obs")))) #sample
    }else{
      hr <- try(hclust(dist(data_m))) #metabolites
      hc <- try(hclust(dist(t(data_m)))) #sample
    }
    
    # hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE)
    # 
    # hr <- try(hclust(as.dist(1-WGCNA::cor(t(data_m),method=cor.method,use="pairwise.complete.obs"))),silent=TRUE)
    
    if(is(hr,"try-error") || is(hc,"try-error")){
								
                                print("Hierarchical clustering can not be performed. ")
    }else{
      
        # heatmap_file<-paste("heatmap2.tiff",sep="")
        # 
        # tiff(heatmap_file,width=plots.width,height=plots.height,res=plots.res)
        
        if(is.data.znorm==FALSE){
                
                heatmap.2(t(data_m), Rowv=as.dendrogram(hc), Colv=as.dendrogram(hr),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="", ColSideColors=patientcolors)
            }else{
                heatmap.2(t(data_m), Rowv=as.dendrogram(hc), Colv=as.dendrogram(hr),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="", ColSideColors=patientcolors)
            }
        }
       # dev.off()
    }
