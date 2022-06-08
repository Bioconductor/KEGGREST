
.getRootUrl <- function()
{
    getOption("KEGG_REST_URL", "https://rest.kegg.jp")
}

.getGenomeUrl <- function()
{
    getOption("KEGG_GENOME_URL", "http://rest.genome.jp")
}

.printf <- function(...) message(noquote(sprintf(...)))

.cleanUrl <- function(url)
{
     url <- gsub(" ", "%20", url, fixed=TRUE)
     url <- gsub("#", "%23", url, fixed=TRUE)
     url <- gsub(":", "%3a", url, fixed=TRUE)
     sub("http(s)*%3a//", "http\\1://", url)
}

.getUrl <- function(url, parser, ...)
{
    url <- .cleanUrl(url)
    debug <- getOption("KEGGREST_DEBUG", FALSE)
    if (debug)
        .printf("url == %s", url)
    response <- GET(url)
    stop_for_status(response)
    content <- .strip(content(response, "text"))
    if (nchar(content) == 0)
        return(character(0))
    do.call(parser, list(content, ...))
}

.strip <- function(str)
{
    gsub("^\\s+|\\s+$", "", str)
}

.rstrip <- function(str)
{
    gsub("\\s+$", "", str)
}

.lstrip <- function(str)
{
    gsub("^\\s+", "", str)
}

.get.tmp.url <- function(url, use.httr=TRUE)
{
    if (use.httr)
    {
        content <- content(GET(url), type="text")
        lines <- strsplit(content, "\n", fixed=TRUE)[[1]]
    } else { ## https://github.com/hadley/httr/issues/27
        t <- tempfile()
        download.file(url, t, quiet=TRUE)
        lines <- readLines(t)
    }
    urlLine <- grep("<img src=\"/tmp", lines, value=TRUE)
    path <- strsplit(urlLine, '"', fixed=TRUE)[[1]][2]
    sprintf("https://www.kegg.jp%s", path)
}

.splitInGroups <- function(x, n)
{
    groups <- seq_len(ceiling(length(x) / n))
    members <- head(rep(groups, each = n), length(x))
    unname(split(x, members))
}
