\name{get.genes.by.motifs}
\alias{get.genes.by.motifs}

\title{Client-side interface to obtain the name of genes that contain
  the motifs represented by a set of motif ids}
\description{
  Given a set of motif ids, the function searches the databases
  implied by the motif ids for genes containing the motifs specified
 by the motif ids.
}
\usage{
get.genes.by.motifs(motif.id.list, start, max.results)
}

\arguments{
  \item{motif.id.list}{\code{motif.id.list} a vector of character strings
    for the ids of the motifs that are conserved by genes across
    organisms}
  \item{start}{\code{start}an integer to indicate the location of the
    entry in the query results from which the results will be
    extracted and returned}
  \item{max.results}{\code{max.results} an integer to indicate the
    maximum number of entries that will be extracted from the query
    results and returned}
}
\details{
  KEGG seems to have two ways of defining the ids for motifs. One is the
  motif ids obtained through \code{\link{get.motifs.by.gene}}, where
  pfam, tfam, pspt, pspf are used for the Pfam, TIGRFAM, PROSITE
  pattern, and PROSITE profile database, respectively and for the first
  part of a motif id (e. g. pfam:aakinase). Another is the motif ids
  used to query the databases for genes that contain the motif, where
  only the first two letters of the abbreviations for databases form the
  first part of a motif id (e. g. pf:aakinase)  
}
\value{
  The function returns a named vector with the names of the vector being
  the textual definition of genes and values of the vector being the ids
  used by KEGG to represent genes
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Jianhua Zhang}

\seealso{\code{\link{get.motifs.by.gene}}}
\examples{
        genes <- get.genes.by.motifs(c("pf:DnaJ", "ps:DNAJ_2"), 1, 10)
}
\keyword{ datasets }


