\name{keggConv}
\alias{keggConv}
\alias{conv}
\alias{bconv}
\title{
Convert KEGG identifiers to/from outside identifiers
}
\description{
Convert KEGG identifiers to/from outside identifiers.
}
\usage{
keggConv(target, source)
}
\arguments{
  \item{target}{
  A KEGG organism code (), T number, or one of the external
  databases \code{ncbi-gi}, \code{ncbi-geneid}, \code{uniprot}, or
  (for chemical substance identifiers) 
  \code{drug}, \code{compound}, or \code{glycan}, \code{pubchem},
  or \code{chebi}.
}

  \item{source}{
  Same as \code{target}, but may also be a list of KEGG identifers
  representing internal or external names.
}

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
keggConv("eco", "ncbi-geneid") ## conversion from NCBI GeneID to KEGG ID for
                               ## E. coli genes
keggConv("ncbi-geneid", "eco") ## opposite direction
keggConv("ncbi-gi", c("hsa:10458", "ece:Z5100")) ## conversion from KEGG ID
                                                 ## to NCBI GI

}
\keyword{ conv }
