\name{GSEA.plot}
\alias{GSEA.plot}

\title{
GSEA.plot performs the GSEA and plots the running sum for one geneset
}
\description{

}
\usage{
GSEA.plot(m, geneset, o = NULL, ga = NULL, nperms = 1000, metric = c("SNR", "FC"), col = rgb(0.3, 0.73, 0.3))
}

\arguments{
  \item{m}{
matrix of expressions. columns represent samples, rows are genes.
}
  \item{geneset}{
vector of genes symbols or ID present in rownames(m)
make sure the geneset is not too small (typically 15) or too big (typically 500)
}
\item{o}{
  2-class factor where length(o) equals ncol(m). If left to NULL, m is
  considered as a paired difference between two classes.
}
  \item{nperms}{
number of permutations; the higher, the better the NES estimation;
  default to 1000
}
  \item{metric}{
metric used SNR to use Signal-to-Noise (snr.SNR function) FC to use Fold
  changes (snr.FC function)
}
  \item{col}{
color of the running sum curve. Matter of taste...
}
}
\details{

}
\value{
returns NULL value since it is primarily used for the plot
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}



\seealso{
GSEA.listofgenes
}
\examples{
## creating synthetic inputs
m<-matrix(rnorm(10*10000),10000,10)
rownames(m)<-paste("r",1:10000)
o<-as.factor(c(rep(1,5),rep(2,5)))

## creating random geneset of size 100
geneset<-paste("r",sample(1:10000,100,rep=F))

## plotting the GSEA running sum
GSEA.plot(m,geneset,o,nperms=100)

## adding some signal to the synthetic data and taking geneset of size
## 100 among the 1000 genes with signal
geneset<-paste("r",sample(1:1000,100,rep=F))
m[1:1000,1:5]<-m[1:1000,1:5]+rnorm(1000*5,mean=2)

## performing the GSEA plot again
##note that the enrichment is negative since class2 is less expressed
##than class 1
GSEA.plot(m,geneset,o,nperms=100)


##PAIRED TEST
##if no factor is provided m is considered as a paired difference
##between class1 & class2
GSEA.plot(m[,1:5]-m[,6:10],geneset,nperms=100)




}
\keyword{ ~kwd1 }
