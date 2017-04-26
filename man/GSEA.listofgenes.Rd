\name{GSEA.listofgenes}
\alias{GSEA.listofgenes}

\title{
GSEA for a list of genesets
}
\description{
GSEA.listofgenes actually takes a list of genesets as input and performs
the GSEA analysis for all of the genesets. Useful for collection of
genesets, e.g. from the Molecular Signatures Database
}
\usage{
GSEA.listofgenes(m, lG, o = NULL, nperms = 1000, metric = c("SNR", "FC"))
}

\arguments{
  \item{m}{
matrix of expressions. columns represent samples, rows are genes.
}
  \item{lG}{
list of gene sets. genes must be in rownames(m).
}
  \item{o}{
2-class factor where length(o) equals ncol(m). If left to NULL, m is
  considered as a paired difference between two classes.
}
  \item{nperms}{
number of permutations to use. Default:1000
}
  \item{metric}{
one of SNR or FC to respectively use snr.SNR or snr.FC function as
  metric in the GSEA analysis
}
}
\details{

}
\value{
returns a list with 3 elements : lRes is a matrix of two columns and
length(lG) rows. For each geneset, column 1 gives the Up-regulation
p-value, column 2 gives the Down-regulation p-value.
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}
\seealso{
GSEA.plot; readgmt;
}
\examples{

##creates matrix of expression with rownames = gene names
m<-matrix(rnorm(10*10000),10000,10)
rownames(m)<-paste("r",1:10000)

## creates 2-level factor
o<-as.factor(c(rep(1,5),rep(2,5)))


## creates list of 10 genesets of size 100
lG<-lapply(1:10,function(x) paste("r",sample(1:10000,100,rep=F)))
names(lG)<-paste("GS",1:10)


##perform GSEA analysis; may take ~1 minute for 1000 permutations;
##depending on your machine of course
res<-GSEA.listofgenes(m,lG,o,1000)
fdr<-computeFDRmedian(res)
str(fdr)



}
