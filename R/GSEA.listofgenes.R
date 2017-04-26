GSEA.listofgenes<-function(m,lG,o=NULL,nperms=1000, metric=c("SNR","FC"))
  {    
    r0<-computeES.lists(m,o,lG,metric[1])
    rR<-lapply(1:nperms,function(x)
               {
                 if(x%%10==0) print(paste(x,"/",nperms))
                 computeES.lists(m=if(is.null(o)) apply(m,2,function(x) if(rnorm(1)<0) -x else x) else m,
                                  o=if(is.null(o)) NULL else sample(o),lG=lG,metric[1])
               })
    lRes<-t(sapply(1:length(r0),function(x)
                 {
                   essmaxp<-sum(r0[[x]]$essmax>sapply(rR,function(y) y[[x]]$essmax))/nperms
                   essminp<-sum(r0[[x]]$essmin<sapply(rR,function(y) y[[x]]$essmin))/nperms
                   return(c(Up.pval=essmaxp,Down.pval=essminp))
                 }))
    rownames(lRes)<-names(lG)
    return(list(lRes=lRes,r0=r0,rR=rR))
  }