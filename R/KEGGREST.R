## Helper for cleaning up things that cannot be unlisted.

.getRootUrl <- function()
{
    getOption("KEGG_REST_URL", "http://rest.kegg.jp")
}

extractFromNamedListObj <- function(x) 
{
    res <- lapply(slotNames(x), slot, object=x)
    setNames(res, slotNames(x))
}

extractFromNamedListObjs <- function(x)
{
    lapply(x, extractFromNamedListObj)
}

## for best neighbors, investigate: 
## http://www.kegg.jp/pathway/eco00260+b0002+c00263
## referenced at http://www.kegg.jp/kegg/rest/weblink.html

get.best.neighbors.by.gene <-
    function(genes.id, start, max.results)
{
    getBestNeighbors(genes.id, start, max.results, "best")
}
## get.best.neighbors.by.gene("eco:b0002",1, 5)

get.best.best.neighbors.by.gene <- function(genes.id, start, max.results)
{
    getBestNeighbors(genes.id, start, max.results, "best_best")
}

getBestNeighbors <- function(genes.id, start, max.results,
                             what = c("best", "best_best"))
{
    stop("not yet implemented.")
    # what <- match.arg(what)

    # if(what == "best"){
    #     neighbors <-  .SOAP(KEGGserver, "get_best_neighbors_by_gene",
    #                         .soapArgs=list('genes_id' = genes.id,
    #                         start = start,  'max_results' = max.results),
    #                         action=KEGGaction, xmlns = KEGGxmlns,
    #                         nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    # }else{
    #     neighbors <-  .SOAP(KEGGserver, "get_best_best_neighbors_by_gene",
    #                         .soapArgs=list('genes_id' = genes.id,
    #                         start = start,
    #                              'max_results' = max.results),
    #                         action=KEGGaction, xmlns = KEGGxmlns,
    #                         nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    # }

    # extractFromNamedListObjs(neighbors)
}
## getBestNeighbors("eco:b0002",1, 5)

get.paralogs.by.gene <- function(genes.id, start, max.results)
{
    stop("not yet implemented.")
    # res <- .SOAP(KEGGserver, "get_paralogs_by_gene",
    #              .soapArgs=list('genes_id' = genes.id, start = start,
    #                'max_results' = max.results),
    #              action=KEGGaction, xmlns = KEGGxmlns,
    #              nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    # ## Only the 1st piece has data in it
    # extractFromNamedListObjs(res)
}
## library(KEGGSOAP); paraGenes <- get.paralogs.by.gene("eco:b0002", 1, 10)

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
    ##.printf("url == %s", url) ## for debugging, remove this later
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

.getX <- function(x, arg)
{
    arg <- paste(arg, collapse="+")
    url <- sprintf("%s/link/%s/%s", .getRootUrl(), x, arg)
    .getUrl(url, .linkParser)
}

.list <- function(arg, ...)
{
    url <- sprintf("%s/list/%s", .getRootUrl(), arg)
    .getUrl(url, .listParser, ...)
}

.conv <- function(arg, organism, ...)
{
    arg <- paste(arg, collapse="+")
    url <- sprintf("%s/conv/%s/%s", .getRootUrl(), organism, arg)
    .getUrl(url, .listParser, valueColumn=1, nameColumn=2)
}

.find <- function(database, arg)
{
    url <- sprintf("%s/find/%s/%s", .getRootUrl(), database, arg)
    .getUrl(url, .listParser, valueColumn=1)
}

.linkParser <- function(txt)
{
    lines <- strsplit(txt, "\n")[[1]]
    splits <- strsplit(lines, "\t")
    sapply(splits, "[[", 2)
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

get.genes.by.pathway <- function(pathway.id)
{
    .getX("genes", pathway.id)
}

get.enzymes.by.pathway <- function(pathway.id)
{
    ## the example returns nothing whereas in KEGGSOAP it returns 14 enzymes!
    .getX("enzyme", pathway.id)
}

get.compounds.by.pathway <- function(pathway.id)
{
    ## example: path:map00010
    .getX("compound", pathway.id)
}

get.reactions.by.pathway <- function(pathway.id)
{
    # example: path:map00010
    .getX("reaction", pathway.id)
}

get.motifs.by.gene <- function(genes.id, db)
{
    stop("not yet implemented.")
  # res <- .SOAP(KEGGserver, "get_motifs_by_gene",
  #              .soapArgs=list('gene_id' = genes.id,
  #                db = db), action=KEGGaction, xmlns = KEGGxmlns,
  #              nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  #  extractFromNamedListObjs(res)
}
## motifs <- get.motifs.by.gene("eco:b0002", "pfam")

## Helper for cleaning up things that cannot be unlisted.
extractFromDefinitions <- function(def)
{
    res <- sapply(def, slot, "entry_id")
    names(res) <- sapply(def, slot, "definition")
    res
}

get.genes.by.motifs <- function(motif.id.list, start, max.results)
{
    stop("not yet implemented.")
  # genes <- .SOAP(KEGGserver, "get_genes_by_motifs",
  #                .soapArgs=list('motif_id_list' = motif.id.list,
  #                  start = start, 'max_results' = max.results),
  #                action=KEGGaction, xmlns = KEGGxmlns,
  #                nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  # extractFromDefinitions(genes)
}

list.databases <- function()
{
    ## There does not seem to be an equivalent call in the REST API.
    ## There is this list of databases in the documentation:
    ##   <database> = pathway | brite | module | disease | drug | environ | ko | genome |
    ##              <org> | compound | glycan | reaction | rpair | rclass | enzyme |
    ##              organism
    ## <org> = KEGG organism code or T number
    ## So, either return that, or don't implement this function.
    stop("Not implemented.")
}

list.organisms <- function()
{
    .list("organism", nameColumn=3, valueColumn=2)
}

list.pathways <- function(org)
{
    .list("pathway", nameColumn=2, valueColumn=1)
}

get.genes.by.organism <- function(org)
{
    ## Not sure if it is worth implementing the start and max.results
    ## arguments as they were originally intended to limit the size of
    ## what comes back from the server, but it looks like we
    ## can't do that with the REST API (unless there is some
    ## HTTP trickery that would allow it), so there seems little
    ## point in trimming what's returned after the fact. The
    ## user can do that just as easily as we can. So I'm removing
    ## those arguments.
    .list(org, valueColumn=1)
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

## NOTE: get.pathways.by.genes() just gives the intersection of the pathways
## based on the genes passed in..
get.pathways.by.genes <- function(genes.id.list)
{
    .getX("pathway", genes.id.list)
}

get.pathways.by.enzymes <- function(enzyme.id.list)
{
    .getX("pathway", enzyme.id.list)
}

get.pathways.by.compounds <- function(compound.id.list)
{
    .getX("pathway", compound.id.list)
}

get.pathways.by.reactions <- function(reaction.id.list)
{
    .getX("pathway", reaction.id.list)
}

search.compounds.by.name <- function(name)
{
    .find("compound", name)
}

search.glycans.by.name <- function(name)
{
    .find("glycan", name)
}

search.compounds.by.composition <- function(composition)
{
    .find("compound", sprintf("%s/formula", composition))
}

search.glycans.by.composition <- function(composition)
{
    .find("glycan", composition)
}

search.compounds.by.mass <- function(mass)
{
    ## range argument goes away, because docs say:
    ## "The exact mass (or molecular weight) is checked by
    ## rounding off to the same decimal place as the query data."
    ## example: /find/compound/174.05/exact_mass  
    ## for 174.045 =< exact mass < 174.055
    .find("compound", sprintf("%s/exact_mass/", mass))
}

search.glycans.by.mass <- function(mass, range)
{
    stop("not yet implemented.")
    # unlist(.SOAP(KEGGserver, "search_glycans_by_mass",
    #              .soapArgs=list('mass' = mass, 'range' = range),
    #              action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

search.compounds.by.subcomp <- function(mol, offset, limit)
{
    stop("not yet implemented.")
    # res <- .SOAP(KEGGserver, "search_compounds_by_subcomp",
    #              .soapArgs=list('mol' = mol, 'offset' = offset, 
    #                'limit' = limit),
    #              action = KEGGaction, xmlns = KEGGxmlns,
    #              nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    # extractFromNamedListObjs(res)
}
## library(KEGGSOAP); mol <- bget("-f m cpd:C00111")
## c4 <- search.compounds.by.subcomp(mol, 1, 5)

search.glycans.by.kcam <- function(kcf, program, option, offset, limit)
{
    stop("not yet implemented.")
  # res <- .SOAP(KEGGserver, "search_glycans_by_kcam",
  #              .soapArgs=list('kcf'=kcf, 'program'=program, 'option'=option,
  #                'offset'=offset, 'limit'=limit),
  #              action = KEGGaction, xmlns = KEGGxmlns,
  #              nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  # extractFromNamedListObjs(res)
}
## library(KEGGSOAP); kcf <- bget("-f k gl:G12922")
## g4 <- search.glycans.by.kcam(kcf, "gapped", "local", 1, 5)

bget <- function(bget.command)
{
    stop("not yet implemented.")
    # unlist(.SOAP(KEGGserver, "bget",
    #              .soapArgs=list('str' = bget.command),
    #              action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}


bconv <- function(id.list, organism)
{
    .conv(id.list, organism)
}

## Support was requested for kegg orthology "ko" numbers.

get.ko.by.gene <- function(genes.id)
{
    .getX("ko", genes.id)
}


## next three need the return value split() into a named vector.

## TODO: compound object returned.  But it should be a simple vector...
get.ko.by.ko.class <- function(ko.class.id)
{
    stop("not yet implemented.")
    #   res<-.SOAP(KEGGserver, "get_ko_by_ko_class",
    #                   .soapArgs=list('ko_class_id' = ko.class.id),
    #                   action=KEGGaction, xmlns = KEGGxmlns,
    #                   nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    # extractFromDefinitions(res)
}
## library(KEGGSOAP); ko <- get.ko.by.ko.class("00524")

get.genes.by.ko.class <- function(ko.class.id, org, offset, limit)
{
    stop("not yet implemented.")
    # res<-unlist(.SOAP(KEGGserver, "get_genes_by_ko_class",
    #                   .soapArgs=list('ko_class_id'=ko.class.id, org=org,
    #                     offset=offset, limit=limit),
    #                   action = KEGGaction, xmlns = KEGGxmlns,
    #                   nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
    # ans <- lapply(res, function(x){x@entry_id})
    # names(ans) <- lapply(res, function(x){x@definition})
    # unlist(ans)
}

get.genes.by.ko  <- function(ko.id, org)
{
    ## how to filter by org?
    ## also, need long name (annotation)
    stop("not yet implemented.")
    # res<-unlist(.SOAP(KEGGserver, "get_genes_by_ko",
    #                   .soapArgs=list('ko_id'=ko.id, org=org),
    #                   action = KEGGaction, xmlns = KEGGxmlns,
    #                   nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
    # ans <- lapply(res, function(x){x@entry_id})
    # names(ans) <- lapply(res, function(x){x@definition})
    # unlist(ans)
}
                      
get.kos.by.pathway <- function(pathway.id)
{
    ## example returns nothing in rest and python,
    ## but returns 36 kos in KEGGSOAP...???
    .getX("ko", pathway.id)
}

get.pathways.by.kos <- function(ko.id.list, org)
{
    ## how to filter by org?
    ## http://rest.kegg.jp/link/pathway/ko%3AK00016+ko%3AK00382
    ## does not return the same thing as KEGGSOAP example
    stop("not yet implemented.")
    # unlist(.SOAP(KEGGserver, "get_pathways_by_kos",
    #              .soapArgs=list('ko_id_list'=ko.id.list, org=org),
    #              action = KEGGaction, xmlns = KEGGxmlns,
    #              nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

