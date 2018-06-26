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
keggConv(target, source, querySize = 100)
}
\arguments{
  \item{target}{
  A KEGG organism code (), T number, or one of the external
  databases \code{ncbi-gi}, \code{ncbi-geneid}, \code{ncbi-proteinid},
  \code{uniprot}, or
  (for chemical substance identifiers) 
  \code{drug}, \code{compound}, or \code{glycan}, \code{pubchem},
  or \code{chebi}.
}

  \item{source}{
  Same as \code{target}, but may also be a list of KEGG identifers
  representing internal or external names.
}

  \item{querySize}{
  Empirically, KEGG limits queries to 100 source identifiers per query.
  This argument enables larger queries by dividing \code{source} into
  sub-queries of no more than \code{querySize} identifiers.
}

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
## conversion from NCBI GeneID to KEGG ID for E. coli genes
head(keggConv("eco", "ncbi-geneid"))
head(keggConv("ncbi-geneid", "eco")) ## opposite direction

## conversion from KEGG ID to NCBI GI
head(keggConv("ncbi-proteinid", c("hsa:10458", "ece:Z5100")))

## conversion from NCBI GI to KEGG ID when the organism code is not known:
head(keggConv("genes", "ncbi-geneid:3113320"))
}
\keyword{ conv }
