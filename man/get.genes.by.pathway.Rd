\name{get.genes.by.pathway}
\alias{get.genes.by.pathway}
\alias{get.enzymes.by.pathway}
\alias{get.compounds.by.pathway}
\alias{get.reactions.by.pathway}

\title{Client-side interface to obtain the KEGG ids for
  genes/enzymes/compounds/reactions that are involved in the
  interactions in a given pathway} 
\description{
  Given a KEGG pathway identifier, the functions query the KEGG PATHWAY
  database for all the genes/enzymes/compounds/reactions that that are
  involved in the interactions in the specified pathway.
}
\usage{
get.genes.by.pathway(pathway.id)
get.enzymes.by.pathway(pathway.id)
get.compounds.by.pathway(pathway.id)
get.reactions.by.pathway(pathway.id) 
}

\arguments{
  \item{pathway.id}{\code{pathway.id} a character string for a KEGG
    pathway id. KEGG pathway ids consist of the string path followed by
    a colon, a three-letter code for the organism of concern, and then
    a number (e. g. "path:eco00020"). The three-letter organism code
    consists of the first letter of the genus name and the first two
    letters of the species name of the scientific name of the organism
    of concern} 
}
\details{
  KEGG pathway identifiers for a given organism can be obtained using
  function \code{\link{list.pathways}}
}
\value{
  The functions return a vector of KEGG gene/enzyme/compound/reation ids
  found in the pathway 
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Jianhua Zhang}

\seealso{\code{\link{list.pathways}}}
\examples{
        genes <- get.genes.by.pathway("path:eco00020")
        enzymes <- get.enzymes.by.pathway("path:eco00020")
        compounds <- get.compounds.by.pathway("path:eco00020")
        reactions <- get.reactions.by.pathway("path:eco00020")
}
\keyword{ datasets }

