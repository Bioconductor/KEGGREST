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
    "image", "kgml"))
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
    to a PNG file, \code{kgml} returns a KGML document.
}
}
\details{
Retrieves all entries from the KEGG database for a set of KEGG identifers.

 \code{keggGet}() can only return 10 result sets at once (this limitation
is on the server side). If you supply more than 10 inputs to \code{keggGet()},
\code{KEGGREST} will warn that only the first 10 results will be returned.
}
\value{
A list wrapping a KEGG flat file.
If \code{option} is \code{aaseq}, an \code{AAStringSet} object.
If \code{option} is \code{ntseq}, a \code{DNAStringSet} object.
If \code{option} is \code{image}, an object which can be written
to a PNG file.
If \code{option} is \code{kgml}, a KGML document.
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
library(png)
writePNG(png, t)
keggGet("hsa05130", "kgml")

}
\keyword{ get }
