\name{get.pathways.by.genes}
\alias{get.pathways.by.genes}
\alias{get.pathways.by.enzymes}
\alias{get.pathways.by.compounds}
\alias{get.pathways.by.reactions}

\title{Client-side interface to obtain the KEGG pathway ids}
\description{
  Given a set of KEGG gene/enzyme/compound/reation identifiers, the
  functions query the KEGG PATHWAY database for all the pathways in
  which items represented by the given set of identifers are involved
}
\usage{ 
get.pathways.by.genes(genes.id.list)
get.pathways.by.enzymes(enzyme.id.list)
get.pathways.by.compounds(compound.id.list)
get.pathways.by.reactions(reaction.id.list)
}

\arguments{
  \item{genes.id.list}{\code{genes.id.list} a vector of character
    strings for the ids used by  KEGG to represent genes. An id
    normally consists of three letters followed by a colon and then
    several numbers. The three letters are from the first letter
    of the genus name and the first two letters of the species name of
    the scientific name of the organism of concern (e. g. hsa:111 for
    Homo Sapiens)} 
  \item{enzyme.id.list}{\code{enzyme.id.list} a vector of character
    strings for enzyme commission numbers}
  \item{compound.id.list}{\code{compound.id.list} a vector of character
    strings for the ids used by  KEGG to represent compounds. A compound id
    begins with cpd: followed by a combination of letters and numbers
    (e. g. cpd:C00579)}
  \item{reaction.id.list}{\code{reaction.id.list} a vector of character
    strings for the ids used by  KEGG to represent reactions. A reaction id
    begins with rn: followed by a combination of letters and numbers
    (e. g. rn:R00268)}
}

\value{
  The functions return a vector of KEGG pathway ids
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Jianhua Zhang}

\seealso{\code{\link{get.genes.by.pathway}},
  \code{\link{get.enzymes.by.pathway}},
  \code{\link{get.compounds.by.pathway}},
  \code{\link{get.reactions.by.pathway}}} 
\examples{
        # There seems to be some problem at the server side. Use try
        pathways <- try(get.pathways.by.genes(c("eco:b0077", "eco:b0078")))
        pathways <- try(get.pathways.by.enzymes("ec:1.3.99.1"))
        pathways <- try(get.pathways.by.compounds(c("cpd:C00033", "cpd:C00158"))) 
        pathways <- try(get.pathways.by.reactions(c("rn:R00959",
                                 "rn:R02740", "rn:R00960", "rn:R01786")))
}
\keyword{ datasets }

