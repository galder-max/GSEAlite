computeFDRmedian<-function(r)
{
  .compute.Stats.median<-function(null.distrib,final.scores, nb.genes, nperms)
    {
      ord.finalscore<-order(final.scores,decreasing=TRUE)
      final.scores.ordered<-final.scores[ord.finalscore]
      ords<-order(c(final.scores.ordered,null.distrib),decreasing=TRUE)
      nb.better<-(match(1:nb.genes,ords)-1:nb.genes)
      pvals<-nb.better
      revert.ind<-match(1:nb.genes,ord.finalscore)    
      pvals<-pvals[revert.ind]
      return(pvals)
    }
  
  .compute.Statsmedian<-function(mrand,final.scores, nb.genes, nperms)
    {
      mBetters<-t(apply(mrand,1,.compute.Stats.median,final.scores,nb.genes,nperms))
      ord<-order(final.scores,decreasing=TRUE)
      ##Up
      expected.nb.better<-apply(mBetters,2,median)
      observed.nb.better<-rank(1/final.scores)
      fdr<-expected.nb.better/observed.nb.better  
      fdr.<-fdr[ord]
      fdr.<-correctlocalQval(fdr.)
      res<-data.frame(score=final.scores[ord],fdr=fdr.)
      return(res)
      revert<-match(1:nb.genes,ord)
      fdr<-fdr.[revert]
      fdr[fdr>1]<-1      
      res<-data.frame(score=final.scores,fdr=fdr)
      return(res)
    }
  ##max
  mES<-sapply(1:length(r$r0),function(y) sapply(r$rR,function(x) x[[y]]$essmax))
  meanES<-colMeans(mES)
  mNES<-t(t(1/mES)*meanES)
  ES<-sapply(r$r0,function(x) x$essmax)
  NES<-meanES/ES
  resUp<-.compute.Statsmedian(mNES,NES,length(NES),nrow(mNES))
  ##rownames(resUp)<-names(r$r0)
  ##min
  mES<-sapply(1:length(r$r0),function(y) sapply(r$rR,function(x) x[[y]]$essmin))
  meanES<-colMeans(mES)
  mNES<-t(t(1/mES)*meanES)
  ES<-sapply(r$r0,function(x) x$essmin)
  NES<-meanES/ES
  resDown<-.compute.Statsmedian(mNES,NES,length(NES),nrow(mNES))
  ##rownames(resDown)<-names(r$r0)
  return(list(resUp=resUp,resDo=resDown))
}
