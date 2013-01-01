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
}
  \item{source}{
}
}
\details{
Many of the old KEGGSOAP functions whose names
started with 'get', such as \code{getFIXME} and \code{getFIXME},
are replaced by using \code{keggLink} (see examples).
}
\value{
A named character vector.
}
\references{
  \url{http://www.kegg.jp/kegg/docs/keggapi.html}
}
\author{
Dan Tenenbaum
}
\examples{
keggLink("pathway", "hsa") ## KEGG pathways linked from each of the human genes
                           ## equivalent to 'get.genes.by.pathway' in KEGGSOAP
keggLink("hsa", "pathway") ## human genes linked from each of the KEGG pathways
                           ## equivalent to 'get.pathways.by.genes' in KEGGSOAP
keggLink("pathway", c("hsa:10458", "ece:Z5100")) ## KEGG pathways linked from a
                                                 ## human gene and an E. coli
                                                 ## O157 gene
                                                 


}
\keyword{ link }
