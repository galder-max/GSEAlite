\name{snr.FC}
\alias{snr.FC}

\title{
snr.FC - Fold Changes
}
\description{
snr.FC computes mean fold changes between two classes in a matrix of
expression. A 2-level factor is used to compare two classes (unpaired). If no factor
is given as input the matrix is considered as a paired difference in
log2 between two groups.
}
\usage{
snr.FC(m, o = NULL)
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

}
\value{
  returns a vector of size nrow(m) correspondig to the fold changes for each gene/row.
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}



\seealso{
snr.SNR
}
\examples{
##creates a matrix of 10000 rows/genes and 10 columns/samples
m<-matrix(rnorm(10*10000),10000,10)


##creates a factor vector of length nrow(m)=10
o<-as.factor(c(rep(1,5),rep(2,5)))

##print results
head(snr.FC(m,o))
head(snr.FC(m))

##compare the SNR metric to the fold changes
plot(snr.SNR(m,o),snr.FC(m,o),xlab="SNR",ylab="FC",main="Comparison
between the two metrics",pch=".")

##same as above but paired
plot(snr.SNR(m,o=NULL),snr.FC(m,o=NULL),xlab="SNR",ylab="FC",main="Comparison
between the two metrics",pch=".")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
