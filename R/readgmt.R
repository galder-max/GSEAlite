readgmt<-function(path.gmt)
  {
    t<-read.csv(path.gmt,sep=" ")
    names.<-sapply(1:nrow(t),function(x)
                   {
                     return(unlist(strsplit(as.character(t[x,]),split="\t"))[1])
                   })
    t<-lapply(1:nrow(t),function(x)
              {
                return(unlist(strsplit(as.character(t[x,]),split="\t"))[-c(1,2)])
              })    
    names(t)<-names.
    return(t)
  }
