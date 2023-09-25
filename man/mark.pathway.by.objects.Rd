\name{mark.pathway.by.objects}
\alias{mark.pathway.by.objects}
\alias{color.pathway.by.objects}

\title{Client-side interface to obtain an url for a KEGG pathway diagram
with a given set of genes marked}
\description{
  Given a KEGG pathway id and a set of KEGG gene ids, the functions
  return the URL of a KEGG pathway diagram with the elements
  corresponding to the genes marked by red or specified color
}
\usage{
mark.pathway.by.objects(pathway.id, object.id.list)
color.pathway.by.objects(pathway.id, object.id.list,
                                     fg.color.list, bg.color.list)
}

\arguments{
  \item{pathway.id}{\code{pathway.id} a character string for a KEGG
    pathway id. KEGG pathway ids consist of the string path followed by
    a colon, a three-letter code for the organism of concern, and then
    a number (e. g. "path:eco00020"). The three-letter organism code
    consists of the first letter of the genus name and the first two
    letters of the species name of the scientific name of the organism
    of concern}
  \item{object.id.list}{\code{object.id.list} a vector of character
    strings for KEGG gene ids. KEGG gene ids normally consist of
    three letters followed by a column and then several numeric
    numbers. The three letters are from the first letter of the genus
    name and the first two letters of the species name of the scientific
    name of the organism of concern (e. g. hsa:111 for Homo Sapiens)}
  \item{fg.color.list}{\code{fg.color.list} a vector of two character
    strings to indicate the color for the text and border, respectively,
    of the objects in a pathway diagram. The strings can either be a
    color code linke #ff0000 or letter link yellow}
  \item{bg.color.list}{\code{bg.color.list} a vector of character
    strings of the same length of \code{object.id.list} to indicate the
    background color of the objects in a pathway diagram. The strings
    can either be a color code like #ff0000 or letter like yellow}
}
\details{
  This function only returns the URL of the KEGG pathway diagram. Use
  the function \code{\link{browseURL}} to view the diagram.

  These functions are not part of the KEGG REST API; they are provided
  because they existed in \code{KEGGSOAP} and an alternative implementation
  was possible.
}
\value{
  This function returns a character string for the url
}
\references{\url{https://www.kegg.jp/kegg/docs/keggapi.html}}
\author{Jianhua Zhang}

\seealso{\code{\link{browseURL}}}
\examples{
 url <- mark.pathway.by.objects(
    "path:eco00260", c("eco:b0002", "eco:c00263")
)
if(interactive()){
    browseURL(url)
}
url <- color.pathway.by.objects(
    "path:eco00260", c("eco:b0002", "eco:c00263"),
    c("#ff0000", "#00ff00"),
    c("#ffff00", "yellow")
)
}
\keyword{ datasets }

