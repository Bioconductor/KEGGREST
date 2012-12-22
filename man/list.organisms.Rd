\name{list.organisms}
\alias{list.organisms}
\alias{list.databases}
\alias{list.pathways}

\title{Client-side interface to obtain the names of organisms
  supported by KEGG databases}
\description{
  These functions provides an R interface to allow users to get the
  names/ids of organisms, databases, pathways that are available through
  KEGG SOAP services.
}

\arguments{
  \item{org}{\code{org} a character string for the id used by KEGG for
    organisms. The organism ids are normally three-letter codes with the
    first letter being the first letter of the genus name and the rest
    being the first two letters of the species name of the scientic name
    of the organism of concern}
}

\usage{
list.organisms()
list.pathways(org)
list.databases()
}

\details{
  Some queries agsinst the KEGG databases require abbreviations of
  organisms supported by KEGG. Although the abbreviations normally
  consist of three letters by truncating the first letter of the genus
  name and the first two letters of the species name (e. g. hsp for Homo
  sapiens), \code{\link{list.organisms}} obtains the abbreviations using
  the service provided by KEGG SOAP to make sure the abbreviations are
  correct and the organisms are indeed supported by KEGG databases. 
}
\value{
  \code{\link{list.organisms}} returns a named vector with names of the
  vector being the scientific names and the values of the vector being the
  abbreviations used by KEGG for the organisms supported by the
  databases.

  \code{\link{list.pathways}} returns a named vector with names of the
  vector being textual descriptions of KEGG pathways and the values of
  the vector being the ids used by KEGG to represent pathways.

  \code{\link{list.databases}} returns a named vector with names of the
  vector being textual descriptions of KEGG databases and the values of
  the vector being the ids used by KEGG to represent the databases.
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Jianhua Zhang}

\examples{
        list.organisms()
}
\keyword{datasets}

