\name{KEGGserver}
\alias{KEGGserver}
\alias{KEGGaction}
\alias{KEGGxmlns}

\title{Definitions of objects used by other functions to access KEGG
  SOAP service}

\description{
  Definitions of KEGG SOAP server, KEGG SOAP action, and KEGG XML name
  space are made when the package is invoked so that they are available
  to other functions 
}

\details{
  All the functions that envoke KEGG SOAP services have a KEGG server,
  KEGG action, and KEGG XML name space as argument. These objects are
  defined in .First.lib and made available to the functions when the
  package is invoked
}

\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Jianhua Zhang}

\keyword{ misc }

