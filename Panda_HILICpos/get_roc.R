get_roc <-
function(dataA,classlabels,classifier="svm",kname="radial",rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=NA,testclasslabels=NA,mainlabel=NA){

    d1<-dataA
    rm(dataA)
    
    d1<-as.data.frame(d1)
    
    #print(dim(d1))
    class_inf<-classlabels
    
    d2<-t(d1)
    
    if(is.na(testset)==TRUE){
        testset<-d1
        testclasslabels<-classlabels
    }
    testset<-t(testset)
    featlist<-rocfeatlist
    featincrement<-rocfeatincrement
    
    
    print(featlist)
    
    d3<-cbind(classlabels,d2)
   
    d3<-as.data.frame(d3)
    
    d3$classlabels <- as.factor(d3$classlabels)
    
    colnames(d3)[1] <- "class"
    
    featlist<-unique(featlist)
   
    mod.lab<-{}
    #  col_vec<-c(heat.colors(256),topo.colors(256))
    
    #col_lab<-col_vec[sample(x=1:length(col_vec),size=length(featlist))] #c("purple","blue","darkgreen", "green","yellow","orange","red")
    
    col_lab<-palette(rainbow(length(featlist)))
    
    if(featincrement==TRUE){
        for(n in 1:length(featlist)){
            
           
            
            num_select<-featlist[n]
            
            if(num_select>dim(d1)[1]){
                
                num_select<-dim(d1)[1]
            }
            
            if(is.na(mainlabel)==TRUE){
                
                
                if(classifier=="logitreg"){
                    mainlab<-"ROC curves using logistic regression"
                }else{
                    
                    mainlab<-"ROC curves using SVM"
                }
                
            }else{
                
                mainlab=mainlabel
            }

            model1 <- svm(class~., data=d3[,c(1:num_select)], type="C",probability=TRUE,kernel=kname)
                
                
            pred<-predict(model1,testset,probability=TRUE,decision.values=TRUE,type="prob")
            pred1 <- prediction(attributes(pred)$probabilities[,2], testclasslabels)
            
            stats1a <- performance(pred1, 'tpr', 'fpr')
            x1<-seq(0,1,0.01)
            y1<-x1
            p1<-performance(pred1,"auc")
            mod.lab <-c(mod.lab,paste('using top ',(num_select),' m/z features: AUC ',round(p1@y.values[[1]],2),sep=""))
            if(n==1){
                plot(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', ylab=stats1a@y.name, xlab=stats1a@x.name, col=col_lab[n], lty=2, main=mainlab,cex.main=1,lwd=2)
            }else{
                lines(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', col=col_lab[n], lty=2,lwd=2)
            }
            
        }
        lines(x1, y1, col="black", type="l",lwd=2)
        legend('bottomright', c(mod.lab), col=col_lab, lwd=c(2,1), lty=1:3, cex=0.8, bty='n')
    }else{

        model1 <- svm(class~., data=d3[,c(featlist)], type="C",probability=TRUE,kernel=kname)
            
        pred<-predict(model1,testset,probability=TRUE,decision.values=TRUE,type="prob")
        pred1 <- prediction(attributes(pred)$probabilities[,2], testclasslabels)

        stats1a <- performance(pred1, 'tpr', 'fpr')
        x1<-seq(0,1,0.01)
        y1<-x1
        p1<-performance(pred1,"auc")
        mod.lab <-c(mod.lab,paste('using top ',(num_select),' m/z features: AUC ',round(p1@y.values[[1]],2),sep=""))
        if(n==1){
            plot(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', ylab=stats1a@y.name, xlab=stats1a@x.name, col=col_lab[n], lty=2, main=mainlab,cex.main=0.8,lwd=2)
        }else{
            lines(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', col=col_lab[n], lty=2,lwd=2)
        }
        
    }
    lines(x1, y1, col="black", type="l",lwd=2)
    legend('bottomright', c(mod.lab), col=col_lab, lwd=c(2,1), lty=1:3, cex=0.8, bty='n')
}
