\name{get.paralogs.by.gene}
\alias{get.paralogs.by.gene}

\title{Client-side interface to obtain data for paralogous genes}
\description{
  Given a KEGG gene id, the function queries the KEGG Sequence Similarity
  Database (SSDB) for genes that are paralogous to the target
  gene. Paralogous genes result from duplication of existing genes and
  then function divergence  
}
\usage{
get.paralogs.by.gene(genes.id, start, max.results)
}

\arguments{
  \item{genes.id}{\code{genes.id} a character string for the id used by
    KEGG to represent the gene of interest. The id normally consists of
    three letters followed by a colon and then several numbers. The
    three letters are from the first letter of the genus name and the
    first two letters of the species name of the scientific name of the
    organism of concern (e. g. hsa:111 for Homo Sapiens)} 
  \item{start}{\code{start} an integer to indicate the location of the
    entry in the query results from which the results will be
    extracted and returned}
  \item{max.results}{\code{max.results} an integer to indicate the
    maximum number of entries that will be extracted from the query
    results and returned}
}
\details{
  A given gene may have several paralogous genes. A
  query to SSDB may have a list of genes that are paralogous to the
  target gene. \code{start} and \code{max.results} indicate where on the
  list to start and stop to extract data and return the results.
}
\value{
  The function returns a list of lists. Each sub-list contains data for
  a gene that is paralogous to the target gene with the following
  elements:

  \item{genes\_id1}{a character string for the id of the target gene used
    to query for hologous genes}
  \item{genes\_id2}{a character string for the id of the homologous gene
    found in another organism}
  \item{sw\_score}{an integer for Smith-Waterman score between genes\_id1
    and genes\_id2}
  \item{bit\_score}{a numeric value for the bit score between genes\_id1
    and genes\_id2} 
  \item{identity}{a numeric value between 0 and 1 for the degree of
    identity between genes\_id1 and genes\_id2}
  \item{overlap}{an integer for the overlapping length between genes\_id1
    and genes\_id2} 
  \item{start\_position1}{an integer for the start position of the
    alignment in genes\_id1} 
  \item{end\_position1}{an integer for the end position of the alignment
    in genes\_id1} 
  \item{start\_position2}{an integer for the start position of the
    alignment in genes\_id2} 
  \item{end\_position2}{an integer for the end position of the alignment
    in genes\_id2} 
  \item{best\_flag\_1to2}{a boolean that is TRUE if genes\_id2 is the best
    neighbor gene of genes\_id1}
  \item{best\_flag\_2to1}{a boolean that is TRUE if genes\_id1 is also
    the best neighbor gene of genes\_id2}
  \item{definition1}{a character string for the definition of genes\_id1}
  \item{definition2}{a character string for the definition of genes\_id2}
  \item{length1}{an integer for the amino acid length of the genes\_id1}
  \item{length2}{an integer for the amino acid length of the genes\_id2}
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Jianhua Zhang}

\seealso{\code{\link{get.best.neighbors.by.gene}}}
\examples{
        paraGenes <- get.paralogs.by.gene("eco:b0002", 1, 10)
}
\keyword{ datasets }

