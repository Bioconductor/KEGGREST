\name{keggLink}
\alias{keggLink}
\alias{link}
\title{
Find related entries by using database cross-references.
}
\description{
Find related entries by using database cross-references.
}
\usage{
keggLink(target, source)
}
\arguments{
  \item{target}{
  Either the name of a single KEGG database (list available via
  \code{\link{listDatabases}()}, a "T number" genome identifier,
  or a KEGG organism code (lists of both available via
  \code{keggList("organism")}).
}
  \item{source}{
  The same as \code{target}, but may also be one or more
  KEGG identifiers.
}
}
\details{
Many of the old KEGGSOAP functions whose names
started with 'get', such as \code{get.pathways.by.genes} and
\code{get.pathways.by.reactions},
are replaced by using \code{keggLink} (see examples).



}
\value{
A named character vector.
}
\references{
  \url{https://www.kegg.jp/kegg/docs/keggapi.html}
}
\author{
Dan Tenenbaum
}
\examples{
res <- keggLink("pathway", "hsa") ## KEGG pathways linked from each of
          ## the human genes equivalent to 'get.genes.by.pathway' in KEGGSOAP
length(res)
head(res)
res <- keggLink("hsa", "pathway") ## human genes linked from each of the
          ## KEGG pathways equivalent to 'get.pathways.by.genes' in KEGGSOAP
keggLink("pathway", c("hsa:10458", "ece:Z5100")) ## KEGG pathways
          ## linked from a human gene and an E. coli O157 gene
res <- keggLink("hsa:126") ## LinkDB search shows all KEGG
          ## resources related to hsa:126
head(res)
}
\keyword{ link }
