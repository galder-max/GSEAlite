computeFDR<-function(r)
  {
    ##max
    mES<-sapply(1:length(r$r0),function(y) sapply(r$rR,function(x) x[[y]]$essmax))
    meanES<-colMeans(mES)
    mNES<-t(t(1/mES)*meanES)
    ES<-sapply(r$r0,function(x) x$essmax)
    NES<-meanES/ES
    resUp<-compute.Stats(as.vector(mNES),NES,length(NES),nrow(mNES))
    rownames(resUp)<-names(r$r0)
    ##min
    mES<-sapply(1:length(r$r0),function(y) sapply(r$rR,function(x) x[[y]]$essmin))
    meanES<-colMeans(mES)
    mNES<-t(t(1/mES)*meanES)
    ES<-sapply(r$r0,function(x) x$essmin)
    NES<-meanES/ES
    resDown<-compute.Stats(as.vector(mNES),NES,length(NES),nrow(mNES))
    rownames(resDown)<-names(r$r0)
    resTot<-cbind(resUp,resDown)
    resTot<-resTot[,c(1,2,9,3,10,4,5,11,12)]
    colnames(resTot)<-c("genesetNB",
                        "scores.up",
                        "scores.down",
                        "rank.up",
                        "rank.down","pvals.up","qvals.up","pvals.do","qvals.do")        
    resUp<-resUp[order(resUp$scores,decreasing=T),]
    resDown<-resDown[order(resDown$scores,decreasing=T),]    
    return(list(resTot=resTot,resUp=resUp,resDo=resDown))
  }