\name{computeES.lists}
\alias{computeES.lists}

\title{
computeES.lists
}
\description{
computes enrichment score for geneset list
}
\usage{
computeES.lists(m,o,lG,metric)
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
  \item{metric}{
one of SNR or FC to respectively use snr.SNR or snr.FC function as
  metric in the GSEA analysis
}
}
\details{

}
\value{
returns a list with length(lG) elements, each of which is a list of two elements returned by the lleading.edge function
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}
\seealso{
GSEA.plot
}
\examples{
##used internally
}
