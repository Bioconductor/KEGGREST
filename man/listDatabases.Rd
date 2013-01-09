\name{listDatabases}
\alias{listDatabases}
\title{
Lists the KEGG databases which may be searched.
}
\description{
Lists the KEGG databases which may be searched. In most cases,
you can also use a KEGG organism name or T number (genome identifier)
as a database name.
}
\usage{
listDatabases()
}
\value{
A character vector of database names.
}
\references{
  \url{http://www.kegg.jp/kegg/docs/keggapi.html}
}
\author{
Dan Tenenbaum
}
\seealso{
\code{\link{keggList}}
}
\examples{
listDatabases()
keggList("organism") ## list all organisms
keggList("hsa") ## list all human genes
keggList("genome") ## list all genome identifiers
keggList("T01001") ## list all human genes
}
\keyword{ database }
\keyword{ databases }
