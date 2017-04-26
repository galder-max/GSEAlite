snr.SNR <- function(m,o=NULL)
{
  if(is.null(o))
     {
       means<-rowMeans(m,na.rm=T)
       N<-ncol(m)
       sd<-rowSums((m-means)^2/N,na.rm=T)
       ## s2n
       s2n <- means/sd
       s2n[sd==0]<-0
       return(s2n)
     }
  islevel1<-as.factor(o) == levels(as.factor(o))[1]
  islevel2<-!islevel1
  n1<-sum(islevel1)
  n2<-sum(islevel2)
  means1<-rowMeans(m[,islevel1],na.rm=T)
  means2<-rowMeans(m[,islevel2],na.rm=T)
  term1<-1/n1+1/n2
  term2<-n1+n2-2
  term<-term1/term2
  n1m1<-n1-1
  n2m1<-n2-1
  sd1<-rowSums(( m[,islevel1]-means1)^2/n1,na.rm=T)
  sd2<-rowSums(( m[,islevel2]-means2)^2/n2,na.rm=T)
  ## s2n
  signals <- means2-means1
  noises<-sqrt(term*(sd1*n1m1+sd2*n2m1))
  s2n <- signals/noises
  s2n[noises==0]<-0
  return(s2n)
}
