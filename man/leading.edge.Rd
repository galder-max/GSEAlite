\name{leading.edge}
\alias{leading.edge}

\title{
lleading edge analysis
}
\description{
performs the leading edge analysis involving the running sum described
in the GSEA methodology. It IS a mess...
}
\usage{
leading.edge<-function(k, 
                       forder,
                       s2n,
                       meanES,
                       Ngenes,
                       doplot=F,
                       label="",
                       alless=NA,
                       col="black",
                       add=F,
                       pas=10,
                       ...)
}

\arguments{
  \item{k}{
boolean vector of length = number of genes. TRUE if the nth gene/row of
  the expression matrix is found in the gene set to be tested.
}
  \item{forder}{
order of the scores between two classes 
}
  \item{s2n}{
scores between two classes : see snr.SNR or snr.FC
}
  \item{meanES}{

}
  \item{Ngenes}{

}
  \item{doplot}{

}
  \item{label}{

}
  \item{alless}{

}
  \item{col}{

}
  \item{add}{

}
  \item{pas}{

}
  \item{...}{

}


}
\details{
used internally for the GSEA plot
}
\value{

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


}
\keyword{ ~kwd1 }
