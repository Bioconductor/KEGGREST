\name{keggCompounds}
\alias{keggCompounds}
\title{
Get list of compounds IDs for pathway
}
\description{
Get list of compounds IDs for pathway.
}
\usage{
keggCompounds(pathwayID)
}
\arguments{
  \item{pathwayID}{
  A KEGG pathway identifier with the prefix \code{map} and 5 digit number.
}

}
\value{
A list of KEGG compound identifiers
}
\references{
  \url{https://www.genome.jp/kegg/pathway.html}
}
\author{
Dan Tenenbaum, Kristina Riemer
}
\examples{
keggCompounds("map00361")
}
\keyword{ compounds }
