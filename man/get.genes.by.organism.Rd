\name{get.genes.by.organism}
\alias{get.genes.by.organism}

\title{Client-side interface to obtain the KEGG ids for all the genes of
  a given organism}
\description{
  Given a KEGG organism id, the function searches the KEGG GENES
  database for all the genes of the organism
}
\usage{
get.genes.by.organism(org)
}

\arguments{
  \item{org}{\code{org} a character string for the id used by KEGG for
    organisms. The organism ids are normally three-letter codes with the
    first letter being the first letter of the genus name and the rest
    being the first two letters of the species name of the scientic name
    of the organism of concern}
}
\details{
  The gene ids returned by the query normally consist of three letters
  followed by a colon and then numbers or a combination of letters and
  numbers. The three letters are from the first letter of the genus name
  and the first two letters of the species name of the scientific name
  of the organism of concern (e. g. hsa:111 for Homo Sapiens) 
}
\value{
  The function returns a vector of character strings of ids used by KEGG
  to represent genes
}
\references{\url{http://www.kegg.jp/kegg/docs/keggapi.html}}
\author{Jianhua Zhang}

\examples{
        genes <- get.genes.by.organism("hsa")
}
\keyword{ datasets }

