\name{keggList}
\alias{keggList}
\title{
Returns a list of entry identifiers and associated definition for a given
database or a given set of database entries.
%%  ~~function to do ... ~~
}
\description{
Returns a list of entry identifiers and associated definition for a given
database or a given set of database entries.
}
\usage{
keggList(database, organism)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database}{
%%     ~~Describe \code{x} here~~
Either a KEGG database (list available via \code{\link{listDatabases}()}),
a KEGG organism code (list available via \code{\link{keggList}()} with the
\code{organism} argument,  a T number (list available via
\code{\link{keggList}()} with the \code{genome} argument), or a character
vector of KEGG identifiers.
}
  \item{organism}{
  \code{Optional.} A KEGG organism identifier (list available via
  \code{\link{keggList}()} with the \code{organism} argument).
}
}
\value{
A named character vector containing entry identifiers and
associated definition.
}
\references{
  \url{https://www.kegg.jp/kegg/docs/keggapi.html}
}
\author{
Dan Tenenbaum
}
\examples{
res <- keggList("pathway") ## returns the list of reference pathways
length(res)
head(res) 
res <- keggList("pathway", "hsa") ## returns the list of human pathways
length(res)
head(res)
res <- keggList("organism") ## returns the list of KEGG organisms with
                     ## taxonomic classification
nrow(res)
head(res)
res <- keggList("hsa")  ## returns the entire list of human genes
length(res)
head(res)
## keggList("T01001") ## same as above
keggList(c("hsa:10458", "ece:Z5100")) ## returns the list of a human gene
                                      ## and an E.coli O157 gene
keggList(c("cpd:C01290","gl:G00092")) ## returns the list of a compound entry
                                      ## and a glycan entry
keggList(c("C01290+G00092")) ## same as above (prefixes are not necessary)
}
\keyword{ list }
