\name{keggGet}
\alias{keggGet}
\title{
Retrieves given database entries
}
\description{
Retrieves given database entries.
}
\usage{
keggGet(dbentries, option = c("aaseq", "ntseq", "mol", "kcf", 
    "image"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbentries}{
  One or more (up to a maximum of 10) KEGG identifiers.
}
  \item{option}{
    \code{Optional.} Option governing the format of the output.
    \code{aaseq} is an amino acid sequence, \code{ntseq} is a nucleotide
    sequence. \code{image} returns an object which can be written
    to a PNG file.
}
}
\value{
A character vector containing a KEGG flat file.
If \code{option} is \code{aaseq} or \code{ntseq}, an
\code{AAStringSet} object.
If \code{option} is \code{image}, an object which can be written
to a PNG file.
}
\references{
  \url{http://www.kegg.jp/kegg/docs/keggapi.html}
}
\author{
Dan Tenenbaum
}
\examples{
keggGet(c("cpd:C01290", "gl:G00092")) ## retrieves a compound entry
                                  ## and a glycan entry
keggGet(c("C01290", "G00092")) ## same as above, without prefixes
keggGet(c("hsa:10458", "ece:Z5100")) ## retrieves a human gene entry
                                     ## and an E.coli O157 gene entry
keggGet(c("hsa:10458", "ece:Z5100"), "aaseq") ## retrieves amino acid sequences
                                              ## of a human gene and an 
                                              ## E.coli O157 gene
png <- keggGet("hsa05130", "image") ## retrieves the image file of a
                                    ## pathway map
t <- tempfile()
writePNG(png, t)

}
\keyword{ get }
