\name{bconv}
\alias{bconv}

\title{Client-side interface to obtain the KEGG ids for
  external gene IDs} 
\description{
  Given a gene identifier, the functions queries KEGG to retrieve the
  appropriate KEGG ID.
}
\usage{
bconv(id.list, organism)
}

\arguments{
  \item{id.list}{a character vector containing the IDs that you wish to
  convert to KEGG IDs.  These IDs must have the appropriate prefix!}
  \item{organism}{Three-letter KEGG organism code, which can be 
  obtained with \code{list.organisms}.}


}
\details{
  Depending on the kind of ID you wish to convert, you must use the
  appropriate prefix followed by a colo and then the correct ID.

  Prefixes supported by KEGG:

  External database  Database prefix
  -----------------  ---------------
  NCBI GI            ncbi-gi:
  NCBI GeneID        ncbi-geneid:
  GenBank            genbank:
  UniGene            unigene:
  UniProt            uniprot:
  
}
\value{
  The functions return a named vector with your initial IDs as the names
  and the appropriate KEGG IDs as the value.
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Marc Carlson}

\examples{
        bconv("ncbi-geneid:10", "hsa")
        bconv(c("ncbi-geneid:100008586", "ncbi-geneid:10"), "hsa")
}
\keyword{ datasets }

