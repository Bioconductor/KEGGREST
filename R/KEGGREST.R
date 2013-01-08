
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
    debug <- getOption("KEGGREST_DEBUG", FALSE)
    if (debug)
        .printf("url == %s", url)
    response <- GET(url)
    result <- http_status(response)
    if (result$category != "success")
        stop(sprintf("invalid request, server returned %s (%s)",
            result$message, url))
        content <- .strip(content(response, "text"))
    if (nchar(content) == 0)
        return(character(0))
    do.call(parser, list(content, ...))
}

.organismListParser <- function(url)
{
    lines <- readLines(url)
    split <- strsplit(lines, "\t")
    u <- unlist(split)
    m <- matrix(u, ncol=4, byrow=TRUE)
    colnames(m) <-  c("T.number", "organism", "species", "phylogeny")
    m
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
    if (database == "organism")
        return(.organismListParser(url))
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


.get_parser_NAME <- function(entry)
{
    ret <- list()
    for (value in names(entry))
    {
        ret[[value]] <- gsub("^;|;$", "", entry[[value]])
    }
    ret
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

.get_parser_REFERENCE <- function(refs)
{
    ret <- list()
    thisref <- list()
    sapply(refs, function(item) {
        if (item$refField == "REFERENCE")
        {
          if (length(thisref) > 0)
            ret <- c(ret, list(thisref))
          thisref <- list(id=item$value)
        } else {
          if (is.null(thisref[[item$refField]]))
            thisref[[item$refField]] <- list()
          thisref[[item$refField]] <- c(thisref[[item$refField]], 
                                        item$value)
        }
    })
    ret <- c(ret, list(thisref))
    ret
}


.get_parser_key_value <- function(entry)
{
    content <- c()
    lines <- strsplit(entry, "\n", fixed=TRUE)[[1]]
    for (line in lines)
    {
        tmp <- strsplit(line, "  ", fixed=TRUE)[[1]]
        key <- tmp[1]
        value <- paste(tmp[2:length(tmp)], collapse="  ")
        content <- c(content, c(.strip(key), .strip(value)))
    }
    content
}

.get_parser_list <- function(entry)
{
    tmp <- paste(entry, collapse=" ")
    strsplit(tmp, "\\s+")[[1]]
}

.flatFileParser <- function(txt)
{
    entry <- list() ## remove?
    refs <- list()
    allEntries <- list()
    last_field <- NULL
    lines <- strsplit(.strip(txt), "\n", fixed=TRUE)[[1]] ## ??
    for (line in lines)
    {
        if (line == "///")
        {
            if("ENTRY" %in% names(entry))
            {
                entry[["ENTRY"]] <- strsplit(entry[["ENTRY"]][1],
                    "\\s+", fixed=TRUE)[[1]]
            }
            if ("NAME" %in% names(entry))
            {
                #####entry[["NAME"]] <- .get_parser_NAME(entry[["NAME"]])
            }

            if (length(refs) > 0)
            {
                entry[["REFERENCE"]] <- .get_parser_REFERENCE(refs)
            }

            # if ("REFERENCE" %in% names(entry))
            # {
            #     entry[["REFERENCE"]] <-
            #         .get_parser_REFERENCE(entry[["REFERENCE"]])
            # }
            for (key in c("REACTION", "ENZYME"))
            {
                if (key %in% names(entry))
                {
                    entry[[key]] <- .get_parser_list(entry[[key]])
                }
            }
            for (key in c("PATHWAY", "ORTHOLOGY"))
            {
                if (key %in% names(entry))
                {
                    entry[[key]] <- .get_parser_key_value(entry[[key]])
                }
            }


            ## dreaded copy-and-append pattern
            allEntries <- c(allEntries, list(entry))
            entry <- list()
            last_field <- NULL
            refs <- list()
            next
        }

        tmp <- strsplit(line, "", fixed=TRUE)[[1]]
        field <- paste(tmp[1:12], collapse="")
        field <- .rstrip(field)
        refField <- .lstrip(field)
        value <- paste(tmp[13:length(tmp)], collapse="")
        value <- .strip(value)

        if (grepl("^ ", field) || field == "")
        ##if (field == "")
            field <- last_field
        else {
            last_field <- field
            entry[[field]] <- c()
        }
        #if (refField == "REFERENCE")
        #    thisref <- list(id=value)

        if (field == "REFERENCE")
        {
            refs <- c(refs,
                list(list(field=field, value=value, refField=refField)))
        }
        entry[[field]] <- c(entry[[field]], value)

    }
    allEntries
}

keggGet <- function(dbentries,
    option=c("aaseq", "ntseq", "mol", "kcf", "image"))
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
            return(readAAStringSet(t))
        }
        if (option %in% c("mol", "kcf"))
            return(.getUrl(url), .textParser)
    }
    .getUrl(url, .flatFileParser)
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
    lines <- strsplit(txt, "\n", fixed=TRUE)[[1]]
    splits <- strsplit(lines, "\t", fixed=TRUE)
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

