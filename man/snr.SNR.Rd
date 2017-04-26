\name{snr.SNR}
\alias{snr.SNR}

\title{
snr.SNR - signal-to-noise ratios
}
\description{
snr.SNR computes the signal-to-noise ratios in a matrix of expression. A
2-level factor is used to compare two classes (unpaired). If no factor
is given as input the matrix is considered as a paired difference in
log2 between two groups.
}
\usage{
snr.SNR(m, o = NULL)
}

\arguments{
  \item{m}{
matrix of expressions. columns represent samples, rows are genes.
}
  \item{o}{
2-class factor where length(o) equals ncol(m). If left to NULL, m is
  considered as a paired difference between two classes.
}
}
\details{
This version of the SNR actually computes a t-like statistics. It was
  named after the SNR used by the original version of GSEA.
}
\value{
returns a vector of size nrow(m) correspondig to the SNR metric for each gene/row.
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{
}

\seealso{
  snr.FC
}
\examples{
##creates a matrix of 10000 rows/genes and 10 columns/samples
m<-matrix(rnorm(10*10000),10000,10)


##creates a factor vector of length nrow(m)=10
o<-as.factor(c(rep(1,5),rep(2,5)))

##print results
head(snr.SNR(m,o))
head(snr.SNR(m))

##compare the SNR metric to the fold changes
plot(snr.SNR(m,o),snr.FC(m,o),xlab="SNR",ylab="FC",main="Comparison
between the two metrics",pch=".")

##same as above but paired
plot(snr.SNR(m,o=NULL),snr.FC(m,o=NULL),xlab="SNR",ylab="FC",main="Comparison
between the two metrics",pch=".")

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }