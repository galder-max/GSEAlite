snr.FC <-
function(m,o=NULL)
  {
    if(is.null(o))
      return(rowMeans(m))
    else
      {
        isl1<-o==levels(as.factor(o))[1]
        isl2<-o==levels(as.factor(o))[2] 
        return(rowMeans(m[,isl2])-rowMeans(m[,isl1]))
      }
  }
