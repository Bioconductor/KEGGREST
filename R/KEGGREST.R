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
    if (database == "organism")
        return(.organismListParser(url))
    .getUrl(url, .listParser, nameColumn=1, valueColumn=2)
}

keggFind <- function(database, query,
    option=c("formula", "exact_mass", "mol_weight"))
{
    if(missing(database))
        stop("'database' argument is required")
    if (!missing(option))
        option <- match.arg(option)
    if (is.integer(query) && length(query) > 1)
        query <- sprintf("%s-%s", min(query), max(query))
    query <- paste(query, collapse="+")
    url <- sprintf("%s/find/%s/%s", .getRootUrl(), database, query)
    if (!missing(option))
        url <- sprintf("%s/%s", url, option)
    .getUrl(url, .listParser, nameColumn=1, valueColumn=2)
}


keggGet <- function(dbentries,
    option=c("aaseq", "ntseq", "mol", "kcf", "image", "kgml"))
{
    if (length(dbentries) > 10)
        warning(paste("More than 10 inputs supplied, only the first",
            "10 results will be returned."))
    dbentries <- paste(dbentries, collapse="+")
    url <- sprintf("%s/get/%s", .getRootUrl(), dbentries)
    if (!missing(option))
    {
        url <- sprintf("%s/%s", url, option)

        if (option == "image")
            return(content(GET(url), type="image/png"))
        if (option %in% c("aaseq", "ntseq"))
        {
            t <- tempfile()
            cat(.getUrl(url, .textParser), file=t)
            if (option == "aaseq")
                return(readAAStringSet(t))
            else if (option == "ntseq")
                return(readDNAStringSet(t))
        }
        if (option %in% c("mol", "kcf", "kgml"))
            return(.getUrl(url, .textParser))
    }
    if (grepl("^br:", dbentries[1]))
        return(.getUrl(url, .textParser))
    .getUrl(url, .flatFileParser)
}

keggCompounds <- function(pathwayID)
{
    url <- sprintf("%s/link/cpd/%s", .getRootUrl(), pathwayID)
    .getUrl(url, .compoundParser)
}

.keggConv <- function(target, source)
{
    query <-paste(source, collapse = "+")
    url <- sprintf("%s/conv/%s/%s", .getRootUrl(), target, query)
    .getUrl(url, .listParser, nameColumn = 1, valueColumn = 2)
}

keggConv <- function (target, source, querySize = 100)
{
    groups <- .splitInGroups(source, querySize)
    answer <- lapply(groups, .keggConv, target = target)
    as(unlist(answer), "character")
}

keggLink <- function(target, source)
{
    if (missing(source))
    {
        url <- sprintf("%s/link/%s",
            .getGenomeUrl(), target)
        .getUrl(url, .matrixParser, ncol=3)
    } else {
        url <- sprintf("%s/link/%s/%s",
            .getRootUrl(), target, paste(source, collapse="+"))
    .getUrl(url, .listParser, nameColumn=1, valueColumn=2)

    }
    ## FIXME?? keggLink("pathway",c("hsa:10458", "ece:Z5100"))
    ## returns a list with duplicate names
}


listDatabases <- function()
{
    c("pathway", "brite", "module", "ko", "genome", "vg", "ag", "compound",
          "glycan", "reaction", "rclass", "enzyme", "disease", "drug",
          "dgroup", "environ", "genes", "ligand", "kegg")
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
    url <- sprintf("https://www.kegg.jp/pathway/%s", pathway.id)
    .get.kegg.url(url)
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

    ## Nov 2020: refactored to use form POST due to issues with long URLs when
    ## large identifier lists are passed.

    pathway.id <- sub("^path:", "", pathway.id)
    if (!(length(object.id.list)==length(fg.color.list) &&
          length(fg.color.list) == length(bg.color.list))) {
        stop(paste("object.id.list, fg.color.list, and bg.color.list must",
            "all be the same length."))
    }

    # format identifier/color list as expected by server
    payload <- paste(
        c("#ids", object.id.list),
        c("cols", paste(bg.color.list, fg.color.list, sep=',')),
        sep="\t",
        collapse="\n"
    )

    # fetch KEGG page from server, via a 302 redirect handled by httr
    # transparently
    res <- POST(
        url = "https://www.kegg.jp/kegg-bin/show_pathway",
        body = list(
            map = pathway.id,
            multi_query = payload,
            mode = 'color'
        ),
        encode="multipart"
    )
    res <- content(res, "text")

    # extract image URL from page
    img_matches <- regexpr(
        "(?<=<img src=\")[^\"]+",
        res,
        perl=T
    )
    img_url <- regmatches(res, img_matches)
    if (length(img_url) < 1) {
        stop(
            "'color.pathway.by.objects()' ",
            "failed to extract KEGG image path from response."
        )
    }
    if (length(img_url) > 1) {
        stop(
            "'color.pathway.by.objects()' ",
            "unexpectedly matched multiple KEGG image paths in response."
        )
    }
    sprintf("https://www.kegg.jp%s", img_url)

}

