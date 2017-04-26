\name{compute.Stats}
\alias{compute.Stats}

\title{
compute.Stats
}
\description{
computes statistics from a vector of null distribution values (expected) and a vector of observed values
}
\usage{
compute.Stats<-function(null.distrib,final.scores, nb.genes, nperms)
}

\arguments{
  \item{null.distrib}{
null distribution = expected vector of values
}
  \item{final.scores}{
scores = observed vector of values
}
  \item{nb.genes}{
number of observed values = length(final.scores)
}
  \item{nperms}{
number of nperms = length(null.distrib)/nb.genes
}
}
\details{

}
\value{
returns a dataframe with 7 columns : 1. the gene's row 2. observed score 3. rank of the observed score 4. pval.up 5. qval.up 6. pval.do 7. qval.do
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
