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
  \url{https://www.kegg.jp/kegg/docs/keggapi.html}
}
\author{
Dan Tenenbaum
}
\seealso{
\code{\link{keggList}}
}
\examples{
listDatabases()
res <- keggList("organism") ## list all organisms
nrow(res)
head(res)
res <- keggList("hsa") ## list all human genes
length(res)
head(res)
## keggList("T01001") ## list all human genes
res <- keggList("genome") ## list all genome identifiers
length(res)
head(res)
}
\keyword{ database }
\keyword{ databases }
