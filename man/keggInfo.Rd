\name{keggInfo}
\alias{keggInfo}
\alias{info}
\title{
Displays the current statistics of a given database
}
\description{
Displays statistics of a given database, such as number of
entries, version, release date, and source. 
}
\usage{
keggInfo(database)
}
\arguments{
  \item{database}{
  Either a KEGG database (list available via \code{\link{listDatabases}()}),
  a KEGG organism code (list available by calling \code{\link{keggList}()})
  with the \code{organism} argument), or a T number (list available by
  calling \code{\link{keggList}()} with the \code{genome} argument.)
  
}
}
\value{
A character vector containing statistics about \code{database}.
}
\references{
  \url{http://www.kegg.jp/kegg/docs/keggapi.html}
}
\author{
Dan Tenenbaum
}
\examples{
keggInfo("kegg") ## displays the current statistics of the KEGG database
keggInfo("pathway") ## displays the number pathway entries including both
                    ## the reference and organism-specific pathways
keggInfo("hsa") ## displays the number of gene entries for the KEGG organism
                ## Homo sapiens
}
\keyword{ info }
\keyword{ metadata }
