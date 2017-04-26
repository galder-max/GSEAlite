\name{lleading.edge}
\alias{lleading.edge}

\title{
lleading edge analysis
}
\description{
performs the leading edge analysis involving the running sum described
in the GSEA methodology.
}
\usage{
lleading.edge(k, s2n)
}

\arguments{
  \item{k}{
boolean vector of length = number of genes. TRUE if the nth gene/row of
  the expression matrix is found in the gene set to be tested.
}
  \item{s2n}{
scores between two classes : see snr.SNR or snr.FC
}
}
\details{

}
\value{
return a list with two elements : essmin is the minimum enrichment score
in the running sum; essmax is the maximum enrichment score in the
running score.
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}



\seealso{
  snr.SNR; snr.FC; GSEA.plot
}

\examples{

## generates random scores
scores<-rnorm(10000)

## generate boolean vector with 50 genes present in the 500 highest
   scores
k<-sample(which(rank(scores)<=500),50,rep=F)

## print the leading edge result of 50 random genes vs. 50 genes in top 500
lleading.edge(k,scores) ##top 500
lleading.edge(sample(k),scores) ##random

}
\keyword{ ~kwd1 }
