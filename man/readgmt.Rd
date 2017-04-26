\name{readgmt}
\alias{readgmt}

\title{
read.gmt
}
\description{
reads a gmt file on your disk to the memory.
}
\usage{
readgmt(path.gmt)
}

\arguments{
  \item{path.gmt}{
    the path to the gmt file on your disk.
}
}
\details{
}
\value{
named list of genesets. Each element is a gene set, i.e. a vector of
  characters representing official gene symbols read from the gmt file
}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{
}

\seealso{
GSEA.listofgenes;
}
\examples{

## collections at: http://www.broadinstitute.org/gsea/msigdb/collections.jsp
path<-"/Users/yourUserName/Desktop/c5.gmt"
lG<-readgmt(path)

##use function GSEA.listofgenes to run GSEA analysis on the collection.
res<-GSEA.listofgenes(m,lG,o,1000)
fdr<-computeFDRmedian(res)
str(fdr)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }