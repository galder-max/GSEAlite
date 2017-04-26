lleading.edge <-function(k, #boolean vector: TRUE if genes in genelist
                         s2n)#scores 
  {
    forder<-order(s2n)
    Ngenes<-length(s2n)
    running.sumt <- vector(length=length(s2n))
    emplacement<- k[forder]
    sum.k.signal<-sum(abs(s2n[k]))
    running.sumt<-s2n[forder]
    running.sumt[!emplacement]<--abs(1/(Ngenes-sum(k)))
    running.sumt[emplacement]<-abs(running.sumt[emplacement] )/sum.k.signal
    running.sum<-cumsum(running.sumt)
    return(list(essmin=min(running.sum),
                essmax=max(running.sum)))
  }
