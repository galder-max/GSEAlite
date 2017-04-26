computeES.lists<-function(m, o, lG, metric)
  {
    lS<-get(paste("snr.",metric,sep=""))(m,o)
    res<-lapply(lG,function(x)
                  {                   
                    lleading.edge(rownames(m)%in%x,lS)
                  })
    return(res)
  }