
.getRootUrl <- function()
{
    getOption("KEGG_REST_URL", "http://rest.kegg.jp")
}



.printf <- function(...) message(noquote(sprintf(...)))

.cleanUrl <- function(url)
{
     url <- gsub(" ", "%20", url, fixed=TRUE)
     url <- gsub("#", "%23", url, fixed=TRUE)
     url <- gsub(":", "%3a", url, fixed=TRUE)
     sub("http%3a//", "http://", url, fixed=TRUE)
}

.getUrl <- function(url, parser, ...)
{
    url <- .cleanUrl(url)
    .printf("url == %s", url) ## for debugging, remove this later
    response <- GET(url)
    result <- http_status(response)
    if (result$category != "success")
        stop(sprintf("invalid request, server returned %s (%s)",
            result$message, url))
    content <- gsub("^\\s+|\\s+$", "", content(response, "text"))
    if (nchar(content) == 0)
        return(character(0))
    do.call(parser, list(content, ...))
}


keggInfo <- function(database)
{
    ## FIXME return an object instead of a character vector
    url <- sprintf("%s/info/%s", .getRootUrl(), database)
    .getUrl(url, .textParser)
}

keggList <- function(database, organism)
{
    database <- paste(database, collapse="+")
    if (missing(organism))
        url <- sprintf("%s/list/%s", .getRootUrl(), database)
    else
        url <- sprintf("%s/list/%s/%s", .getRootUrl(), database, organism)
    .getUrl(url, .listParser, nameColumn=1, valueColumn=2)
}

keggFind <- function(database, query,
    option=c("formula", "exact_mass", "mol_weight"))
{
    if (is.integer(query) && length(query) > 1)
        query <- sprintf("%s-%s", min(query), max(query))
    query <- paste(query, collapse="+")
    url <- sprintf("%s/find/%s/%s", .getRootUrl(), database, query)
    if (!missing(option))
        url <- sprintf("%s/%s", url, option)
    .getUrl(url, .listParser, nameColumn=1, valueColumn=2)
}

keggGet <- function(dbentries,
    option=c("aaseq", "ntseq", "mol", "kcf", "image"))
{
    dbentries <- paste(dbentries, collapse="+")
    url <- sprintf("%s/get/%s", .getRootUrl(), dbentries)
    if (!missing(option))
        url <- sprintf("%s/%s", url, option)
    if (!missing(option))
    {
        if (option == "image")
            return(content(GET(url), type="image/png"))
        if (option == "aaseq")
            ## FIXME returns amino acid sequence in FASTA format, 
            ## what's the best BioC class for this?
            return(.getUrl(url, .textParser)) 
    }
    ## FIXME, convert KEGG flat file to something useful
    .getUrl(url, .textParser)
}

keggConv <- function(target, source)
{
    url <- sprintf("%s/conv/%s/%s",
        .getRootUrl(), target, paste(source, collapse="+"))
    .getUrl(url, .listParser, nameColumn=1, valueColumn=2)
}


keggLink <- function(target, source)
{
    url <- sprintf("%s/link/%s/%s",
        .getRootUrl(), target, paste(source, collapse="+"))
    .getUrl(url, .listParser, nameColumn=1, valueColumn=2)
    ## FIXME?? keggLink("pathway",c("hsa:10458", "ece:Z5100"))
    ## returns a list with duplicate names
}


.listParser <- function(txt, valueColumn, nameColumn)
{
    lines <- strsplit(txt, "\n")[[1]]
    splits <- strsplit(lines, "\t")
    ret <- sapply(splits, "[[", valueColumn)
    if (!missing(nameColumn)) {
        nms <- sapply(splits, "[[", nameColumn)
        names(ret) <- nms
    }
    ret
}

.textParser <- function(txt)
{
    txt
}



listDatabases <- function()
{
    ## FIXME ADD unit test to ensure these all exist
    c("pathway", "brite", "module", "disease", "drug", "environ",
        "ko", "genome", "compound", "glycan", "reaction", "rpair",
        "rclass", "enzyme", "organism")
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
    urlLine <- grep("<img src=\"/tmp", lines, value=T)
    path <- strsplit(urlLine, '"', fixed=TRUE)[[1]][2]
    sprintf("http://www.kegg.jp%s", path)
}


## This is not strictly speaking an API supported by the KEGG REST
## server, but it seems useful, and does not use SOAP, so I'm leaving it in.
mark.pathway.by.objects <- function(pathway.id, object.id.list)
{
    ## example: http://www.kegg.jp/pathway/eco00260+b0002+c00263
    pathway.id <- sub("^path:", "", pathway.id)
    if (!missing(object.id.list)) {
        object.id.list <- paste(object.id.list, collapse="+")
        pathway.id <- sprintf("%s+%s", pathway.id, object.id.list)
    }
    url <- sprintf("http://www.kegg.jp/pathway/%s", pathway.id)
    .get.tmp.url(url)
}

## This is not strictly speaking an API supported by the KEGG REST
## server, but it seems useful, and does not use SOAP, so I'm leaving it in.
color.pathway.by.objects <- function(pathway.id, object.id.list,
                                     fg.color.list, bg.color.list)
{
    ## example: http://www.kegg.jp/kegg-bin/show_pathway?eco00260/b0002%09%23ff0000,%2300ff00/c00263%09%23ffff00,yellow
    ## also works to include organism code in gene IDs
    ## (but don't include path: in pathway id)
    ## documentation here: http://www.kegg.jp/kegg/rest/weblink.html
    ## and here: http://www.kegg.jp/kegg/tool/map_pathway2.html
    pathway.id <- sub("^path:", "", pathway.id)
    if (!(length(object.id.list)==length(fg.color.list) &&
          length(fg.color.list) == length(bg.color.list))) {
        stop(paste("object.id.list, fg.color.list, and bg.color.list must",
            "all be the same length."))
    }
    url <- sprintf("http://www.kegg.jp/kegg-bin/show_pathway?%s/", pathway.id)
    segs <- sprintf("%s%%09%s,%s", object.id.list,
                    fg.color.list, bg.color.list)
    url <- sprintf("%s%s", url, paste(segs, collapse="/"))
    url <- .cleanUrl(url)
    .get.tmp.url(url, use.httr=FALSE)
}

