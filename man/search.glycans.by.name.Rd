\name{search.glycans.by.name}
\alias{search.glycans.by.name}
\alias{search.glycans.by.mass}
\alias{search.glycans.by.composition}
\alias{search.glycans.by.kcam}

\title{Client-side interface to obtain a list of chemical glycans}
\description{
  The functions provide access to KEGG LIGAND database 
  \url{http://www.genome.jp/kegg/ligand.html}. Given a glycan name, 
  a composition, a molecular weight, or a common sub-structure, 
  one of the functions below can return a list of glycans 
  identifiers from KEGG LIGAND database.
}
\usage{
search.glycans.by.name(name)
search.glycans.by.composition(composition)
search.glycans.by.mass(mass, range)
search.glycans.by.kcam(kcf, program, option, offset, limit)
}

\arguments{
  \item{name}{\code{name} a character string to indicate a glycan
    name}
  \item{composition}{\code{composition} a character string to indicate 
    the composition of monosaccharides} 
  \item{mass}{\code{mass} a float to indicate the mass computed from
    the composition, excluding those in parentheses}
  \item{range}{\code{range} a float to indicate the range of molecular
    weight when searching glycans by mass}
  \item{kcf}{\code{kcf} a character string to indicate the molecular
    structure (carbohydrate sequence) of a glycan in KCF format} 
  \item{program}{\code{program} a character string, either "gapped"
    or "ungaped"}
  \item{option}{\code{option} a character string, either "global" or
    "local"}
  \item{offset}{\code{offset} an integer}
  \item{limit}{\code{limit} an integer}  
}
\details{
  \code{search.glycans.by.name} returns a list of glycans having the
  specified name; 

  \code{search.glycans.by.composition} returns a list
  of glycans containing sugars indicated by the composition. Order
  of the sugars (in parenthesis with number) is insensitive; 

  \code{search.glycans.by.mass} returns a list of glycans having the
   molecular weight around "mass" with some ambiguity (range); 

   \code{search.glycans.by.subcomp} returns a list of glycans with 
  the alignment having common sub-structure calculated by the KCaM 
  program. You can obtain a KCF formatted structural data of matched 
  glycans using \code{bget} with the "-f m" option to confirm the alignment. 
}
\value{
  All the functions return a character vector of glycan identifiers
  provided by KEGG LIGAND database
}
\references{\url{http://www.genome.jp/kegg/docs/keggapi_manual.html#label:105}}
\author{Nianhua Li}

\seealso{\code{\link{bget}}}
\examples{
	glycans_1 <- search.glycans.by.name("Paragloboside")
        glycans_2 <- search.glycans.by.composition("(Man)4 (GalNAc)1")
        glycans_3 <- search.glycans.by.mass(689.6, 0.1)
        kcf <- bget("-f k gl:G12922")
        glycans_4 <- search.glycans.by.kcam(kcf, "gapped", "local", 1, 5)
}
\keyword{manip}


