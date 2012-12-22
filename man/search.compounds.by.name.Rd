\name{search.compounds.by.name}
\alias{search.compounds.by.name}
\alias{search.compounds.by.mass}
\alias{search.compounds.by.composition}
\alias{search.compounds.by.subcomp}

\title{Client-side interface to obtain a list of chemical compounds}
\description{
  The functions provide access to KEGG LIGAND database 
  \url{http://www.genome.jp/kegg/ligand.html}. Given a compound name, 
  a chemical formula, a molecular weight, or a common sub-structure, 
  one of the functions below can return a list of compounds 
  identifiers from KEGG LIGAND database.
}
\usage{
search.compounds.by.name(name)
search.compounds.by.composition(composition)
search.compounds.by.mass(mass, range)
search.compounds.by.subcomp(mol, offset, limit)
}

\arguments{
  \item{name}{\code{name} a character string to indicate a compound
    name}
  \item{composition}{\code{composition} a character string to indicate 
    a compound composition, usually expressed as chemical formula} 
  \item{mass}{\code{mass}a float to indicate a molecular weight
    around mass}
  \item{range}{\code{range} a float to indicate the range of molecular
    weight when searching compounds by mass}
  \item{mol}{\code{mol}a character string to indicate a MOL formatted
    structural data, more in details section}  
  \item{offset}{\code{offset} an integer}
  \item{limit}{\code{limit} an integer}  
}
\details{
  \code{search.compounds.by.name} returns a list of compounds having the
  specified name; 
  
  \code{search.compounds.by.composition} returns a list
  of compounds containing elements indicated by the composition. Order
  of the elements is insensitive; 

  \code{search.compounds.by.mass} returns
  a list of compounds having the molecular weight around "mass" with
  some ambiguity (range); 

  \code{search.compounds.by.subcomp} returns
  a list of compounds with the alignment having common sub-structure
  calculated by the subcomp program. You can obtain a MOL formatted
  structural data of matched compounds using \code{bget} with the
  "-f m" option to confirm the alignment. 
}
\value{
  All the functions return a character vector of chemical compound identifiers
  provided by KEGG LIGAND database
}
\references{\url{http://www.genome.jp/kegg/docs/keggapi_manual.html#label:105}}
\author{Nianhua Li}

\seealso{\code{\link{bget}}}
\examples{
	compounds_1 <- search.compounds.by.name("shikimic acid")
        compounds_2 <- search.compounds.by.composition("C7H10O5")
        compounds_3 <- search.compounds.by.mass(174.05, 0.1)
        mol <- bget("-f m cpd:C00111")
        compounds_4 <- search.compounds.by.subcomp(mol, 1, 5)
}
\keyword{ datasets }


