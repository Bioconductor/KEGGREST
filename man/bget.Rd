\name{bget}
\alias{bget}

\title{Client-side interface to obtain KEGG database entries by a list
  of entry identifiers}
\description{
  \code{bget} is used for retrieving KEGG database entries specified by
  a list of entry identifiers. It accepts all the KEGG bget commond line
  options as a character string. Number of entries retrieved at a time
  is restricted up to 100.
}

\arguments{
  \item{bget.command}{\code{bget.command} a character string of KEGG
    bget command}
}

\usage{
bget(bget.command)
}

\value{
  a character string of KEGG bget search result.
}
\references{\url{http://www.genome.jp/kegg/docs/keggapi_manual.html#label:40}}
\author{Nianhua Li}

\examples{
	# retrieve two KEGG/GENES entries
	bget("eco:b0002 hin:tRNA-Cys-1")
	# retrieve nucleic acid sequences in a FASTA format
	bget("-f -n n eco:b0002 hin:tRNA-Cys-1")
	# retrieve amino acid sequence in a FASTA format
	bget("-f -n a eco:b0002")
}
\keyword{datasets}

