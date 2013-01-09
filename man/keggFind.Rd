\name{keggFind}
\alias{keggFind}
\title{
Finds entries with matching query keywords or other query data in a given 
database
}
\description{
Finds entries with matching query keywords or other query data in a given 
database.
}
\usage{
keggFind(database, query, option = c("formula", "exact_mass", 
    "mol_weight")) 
}
\arguments{
  \item{database}{
  Either the name of a single KEGG database (list available via
  \code{\link{listDatabases}()}, a "T number" genome identifier,
  or a KEGG organism code (lists of both available via
  \code{keggList("organism")}).
}
  \item{query}{
  One or more keywords, or a range of integers representing 
  molecular weights.
  If \code{query} includes identifiers not known to KEGG, 
  the results will not contain any information about those identifiers.
}
  \item{option}{
    \code{Optional.} If \code{database} is \code{compound} or \code{drug},
    \code{option} can be \code{formula}, \code{exact_mass}, or 
    \code{weight}. 
    Chemical formula search is a partial match irrespective of the
    order of atoms given. 
    The exact mass (or molecular weight) is checked by rounding off to the
    same decimal place as the query data.
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
keggFind("genes", c("shiga", "toxin")) ## for keywords "shiga" and "toxin"
keggFind("genes", "shiga toxin") ## for keywords "shiga toxin"
keggFind("compound", "C7H10O5", "formula") ## for chemical formula "C7H10O5"
keggFind("compound", "O5C7", "formula") ## for chemical formula containing
                                        ## "O5" and "C7"
keggFind("compound", 174.05, "exact_mass") ## for 174.045
                                           ## =< exact mass < 174.055
keggFind("compound", 300:310, "mol_weight") ## for 300 =<
                                            ## molecular weight =< 310
}
\keyword{ find }
