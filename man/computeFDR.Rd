\name{computeFDR}
\alias{computeFDR}

\title{
  computes False Discovery Rate (FDR)
}
\description{
computes FDR from scores and permutation scores
}
\usage{
computeFDR(r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
result list from GSEA analysis with r0 scores; rR random scores 
}
}
\details{
computes FDR from scores and permutation scores
}
\value{
returns a list with two dataframes; one for Up regulation and one for
  Down regulation
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}



\seealso{

}
\examples{
     ##creates matrix of expression with rownames = gene names
     m<-matrix(rnorm(10*10000),10000,10)
     rownames(m)<-paste("r",1:10000)    
     
     ## creates list of 10 genesets of size 100
     lG<-lapply(1:10,function(x) paste("r",sample(1:10000,100,rep=F)))
     names(lG)<-paste("GS",1:10)
     
     
     ##performs GSEA analysis
     r<-GSEA.listofgenes(m,lG,nperms=100)

     ##computesFDR for the results
     computeFDR(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
