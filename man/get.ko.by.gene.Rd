\name{get.ko.by.gene}
\alias{get.ko.by.gene}
\alias{get.ko.by.ko.class}
\alias{get.genes.by.ko.class}
\alias{get.genes.by.ko}
\alias{get.kos.by.pathway}
\alias{get.pathways.by.kos}

\title{Client-side interfaces to obtain the KEGG ko ids for a pathway
  and vice versa}
\description{
  Given a KEGG pathway ko identifier, the functions query the KEGG
  PATHWAY database for all the pathway id or vice versa.
}
\usage{
get.ko.by.gene(genes.id)
get.ko.by.ko.class(ko.class.id)
get.genes.by.ko.class(ko.class.id, org , offset, limit)
get.genes.by.ko(ko.id, org)
get.kos.by.pathway(pathway.id)
get.pathways.by.kos(ko.id.list, org)
}

\arguments{
  \item{genes.id}{a vector of gene IDs}
  \item{ko.id}{a vector of ko IDs}
  \item{ko.class.id}{a vector of ko class IDs}
  \item{pathway.id}{\code{pathway.id} a character string for a KEGG
    pathway id. KEGG pathway ids consist of the string path followed by
    a colon, a three-letter code for the organism of concern, and then
    a number (e. g. "path:eco00020"). The three-letter organism code
    consists of the first letter of the genus name and the first two
    letters of the species name of the scientific name of the organism
    of concern}
  \item{ko.id.list}{\code{pathway.id} a vector of KEGG
    ko IDs.}
  \item{org}{\code{pathway.id} a string containing the three letter
    KEGG prefix to use in looking up the IDs}
  \item{offset}{an offset}
  \item{limit}{how many}
}
\value{
  The functions return a vector or a named list of values depending on
  what the function is supposed to retrieve.
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Marc Carlson}

\seealso{\code{\link{list.pathways}}}
\examples{
	ko <- get.ko.by.gene("eco:b0002")
	ko <- get.ko.by.ko.class("00524")
	genes <- get.genes.by.ko.class("00903", "hsa" , 1, 100)
	genes <- get.genes.by.ko("ko:K12524", "eco")
        kos <- get.kos.by.pathway("path:hsa00010")
        pathways <- get.pathways.by.kos(c("ko:K00016","ko:K00382"), "hsa")
}
\keyword{ datasets }