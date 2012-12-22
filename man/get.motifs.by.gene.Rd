\name{get.motifs.by.gene}
\alias{get.motifs.by.gene}

\title{Client-side interface to obtain the name of genes that are
  homologous to a given gene}
\description{
  This function queries the Pfam, TIGRFAM, PROSITE pattern, and/or
  PROSITE profile databases for the motifs of a given gene. A motif is a
  locally conserved region of a sequence or a short sequence pattern
  shared by a set of sequences
}
\usage{
get.motifs.by.gene(genes.id, db)
}

\arguments{
  \item{genes.id}{\code{genes.id} a character string for the id used by
    KEGG to represent the gene of interest. The id normally consists of
    three letters followed by a colon and then several numbers. The
    three letters are from the first letter of the genus name and the
    first two letters of the species name of the scientific name of the
    organism of concern (e. g. hsa:111 for Homo Sapiens)} 
  \item{db}{\code{db} a character string for the name of the data to
    search for motifs. Valid database names include pfam, tfam, pspt,
    pspf for the Pfam, TIGRFAM, PROSITE pattern, and PROSITE profile,
    respectively, or all for all the four databases}
}
\details{
  The motif ids obtained can be used to search for the genes that
  contain the motif across organism using \code{\link{get.genes.by.motifs}}
}
\value{
  The function returns a list of lists with each of the sub-list having
  the following elements:
  \item{motif_id}{a character string for the id of the motif found} 
  \item{definition}{a character string for the definition of the motif}
  \item{genes_id}{a character string for the KEGG genes\_id of the gene
    that contains the motif and used to search the database(s)} 
  \item{start_position}{an integer for the start position of the motif match}
  \item{end.position}{an integer for the end position of the motif match}
  \item{score}{a numeric value for the score of the motif match for
    TIGRFAM and PROSITE databases}
  \item{evalue}{a numeric value for the E-value of the motif match for
    Pfam database}
}
\references{\url{http://www.genome.jp/kegg/soap/doc/keggapi_manual.html}}
\author{Jianhua Zhang}

\seealso{\code{\link{get.genes.by.motifs}}}
\examples{
        motifs <- get.motifs.by.gene("eco:b0002", "pfam")
}
\keyword{ datasets }

