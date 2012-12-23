#KEGGserver <- SOAPServer("http://soap.genome.ad.jp/keggapi/request.cgi")
#KEGGxmlns = "SOAP/KEGG"
#KEGGaction = "SOAP/KEGG"

## Helper for cleaning up things that cannot be unlisted.

.getRootUrl <- function()
{
    getOption("KEGG_REST_URL", "http://rest.kegg.jp")
}

extractFromNamedListObj <- function(x){
  names <- slotNames(x)
  res <- list()
  for(i in seq_len(length(names))){
    res[[i]] <- slot(x, names[i])
  }
  names(res)<- names
  res  
}

extractFromNamedListObjs <- function(x){
  lapply(x, extractFromNamedListObj)
}

## for best neighbors, investigate: 
## http://www.kegg.jp/pathway/eco00260+b0002+c00263
## referenced at http://www.kegg.jp/kegg/rest/weblink.html

get.best.neighbors.by.gene <- function(genes.id, start, max.results){
    getBestNeighbors(genes.id, start, max.results, "best")
}
## get.best.neighbors.by.gene("eco:b0002",1, 5)

get.best.best.neighbors.by.gene <- function(genes.id, start, max.results){
    getBestNeighbors(genes.id, start, max.results, "best_best")
}

getBestNeighbors <- function(genes.id, start, max.results,
                             what = c("best", "best_best")){
    what <- match.arg(what)

    if(what == "best"){
        neighbors <-  .SOAP(KEGGserver, "get_best_neighbors_by_gene",
                            .soapArgs=list('genes_id' = genes.id,
                            start = start,  'max_results' = max.results),
                            action=KEGGaction, xmlns = KEGGxmlns,
                            nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    }else{
        neighbors <-  .SOAP(KEGGserver, "get_best_best_neighbors_by_gene",
                            .soapArgs=list('genes_id' = genes.id,
                            start = start,
                                 'max_results' = max.results),
                            action=KEGGaction, xmlns = KEGGxmlns,
                            nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    }

    extractFromNamedListObjs(neighbors)
}
## getBestNeighbors("eco:b0002",1, 5)


get.paralogs.by.gene <- function(genes.id, start, max.results){
    res <- .SOAP(KEGGserver, "get_paralogs_by_gene",
                 .soapArgs=list('genes_id' = genes.id, start = start,
                   'max_results' = max.results),
                 action=KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    ## Only the 1st piece has data in it
    extractFromNamedListObjs(res)
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
    .printf("url == %s", url) ## for debugging, remove this later
    response <- GET(url)
    result <- http_status(response)
    if (!result$category == "success")
        stop(sprintf("Invalid request, server returned %s. (%s)",
            result$message, url))
    content <- gsub("^\\s+|\\s+$", "", content(response, "text"))
    if (nchar(content) == 0)
        return(character(0))
    do.call(parser, list(content, ...))
}

.get.x <- function(x, arg)
{
    if (length(arg) > 1)
        arg <- paste(arg, collapse="+")
    url <- sprintf("%s/link/%s/%s", .getRootUrl(), x, arg)
    .getUrl(url, .link.parser)
}

.list <- function(arg, ...)
{
    url <- sprintf("%s/list/%s", .getRootUrl(), arg)
    .getUrl(url, .list.parser, ...)
}


.find <- function(database, arg)
{
    url <- sprintf("%s/find/%s/%s", .getRootUrl(), database, arg)
    .getUrl(url, .list.parser, valueColumn=1)
}

.link.parser <- function(txt)
{
    lines <- strsplit(txt, "\n")[[1]]
    ret <- character()
    for (line in lines)
    {
      ret <- c(ret, strsplit(line, "\t")[[1]][2])
    }
    paste(ret, collapse="\n")
    ret
}

.list.parser <- function(txt, nameColumn, valueColumn)
{
    lines <- strsplit(txt, "\n")[[1]]
    ret <- character()
    nms <- character()
    for (line in lines)
    {
        segs <- strsplit(line, "\t")[[1]]
        ret <- c(ret, segs[valueColumn])
        if (!missing(nameColumn))
            nms <- c(nms, segs[nameColumn])
    }
    paste(ret, collapse="\n")
    if (!missing(nameColumn))
        names(ret) <- nms
    ret
}

get.genes.by.pathway <- function(pathway.id) {
    .get.x("genes", pathway.id)
}


get.enzymes.by.pathway <- function(pathway.id){
    ## the example returns nothing whereas in KEGGSOAP it returns 14 enzymes!
    .get.x("enzyme", pathway.id)
}



get.compounds.by.pathway <- function(pathway.id){
    ## example: path:map00010
    .get.x("compound", pathway.id)
}

get.reactions.by.pathway <- function(pathway.id){
    # example: path:map00010
    .get.x("reaction", pathway.id)
}


get.motifs.by.gene <- function(genes.id, db){
  res <- .SOAP(KEGGserver, "get_motifs_by_gene",
               .soapArgs=list('gene_id' = genes.id,
                 db = db), action=KEGGaction, xmlns = KEGGxmlns,
               nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
   extractFromNamedListObjs(res)
}
## motifs <- get.motifs.by.gene("eco:b0002", "pfam")


## Helper for cleaning up things that cannot be unlisted.
extractFromDefinitions <- function(def){
  res <- sapply(def, slot, "entry_id")
  names(res) <- sapply(def, slot, "definition")
  res
}


get.genes.by.motifs <- function(motif.id.list, start, max.results){
  genes <- .SOAP(KEGGserver, "get_genes_by_motifs",
                 .soapArgs=list('motif_id_list' = motif.id.list,
                   start = start, 'max_results' = max.results),
                 action=KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  extractFromDefinitions(genes)
}

list.databases <- function(){
    ## There does not seem to be an equivalent call in the REST API.
    ## There is this list of databases in the documentation:
    ##   <database> = pathway | brite | module | disease | drug | environ | ko | genome |
    ##              <org> | compound | glycan | reaction | rpair | rclass | enzyme |
    ##              organism
    ## <org> = KEGG organism code or T number
    ## So, either return that, or don't implement this function.
}

list.organisms <- function(){
    .list("organism", nameColumn=3, valueColumn=2)
}

list.pathways <- function(org){
  .list("pathway", nameColumn=2, valueColumn=1)
}



get.genes.by.organism <- function(org){
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


mark.pathway.by.objects <- function(pathway.id, object.id.list){
    ## example: http://www.kegg.jp/pathway/eco00260+b0002+c00263
    pathway.id <- sub("^path:", "", pathway.id)
    if (!missing(object.id.list))
    {
        object.id.list <- paste(object.id.list, collapse="+")
        pathway.id <- sprintf("%s+%s", pathway.id, object.id.list)
    }
    url <- sprintf("http://www.kegg.jp/pathway/%s", pathway.id)
    .get.tmp.url(url)
}

color.pathway.by.objects <- function(pathway.id, object.id.list,
                                     fg.color.list, bg.color.list){
    ## example: http://www.kegg.jp/kegg-bin/show_pathway?eco00260/b0002%09%23ff0000,%2300ff00/c00263%09%23ffff00,yellow
    ## also works to include organism code in gene IDs
    ## (but don't include path: in pathway id)
    ## documentation here: http://www.kegg.jp/kegg/rest/weblink.html
    ## and here: http://www.kegg.jp/kegg/tool/map_pathway2.html
    pathway.id <- sub("^path:", "", pathway.id)
    if (!(length(object.id.list)==length(fg.color.list) &&
        length(fg.color.list) == length(bg.color.list)))
    {
        stop(paste("object.id.list, fg.color.list, and bg.color.list must",
            "all be the same length."))
    }
    url <- sprintf("http://www.kegg.jp/kegg-bin/show_pathway?%s/", pathway.id)
    segs <- character(0)
    for (i in 1:length(object.id.list))
    {
        segs <- c(segs, sprintf("%s%%09%s,%s",
            object.id.list[i],
            fg.color.list[i],
            bg.color.list[i]))
    }
    url <- sprintf("%s%s", url, paste(segs, collapse="/"))
    url <- .cleanUrl(url)
    .get.tmp.url(url, use.httr=FALSE)
}

## NOTE: get.pathways.by.genes() just gives the intersection of the pathways
## based on the genes passed in..
get.pathways.by.genes <- function(genes.id.list){
    .get.x("pathway", genes.id.list)
}

get.pathways.by.enzymes <- function(enzyme.id.list){
    .get.x("pathway", enzyme.id.list)
}

get.pathways.by.compounds <- function(compound.id.list){
    .get.x("pathway", compound.id.list)
}

get.pathways.by.reactions <- function(reaction.id.list){
    .get.x("pathway", reaction.id.list)
}


search.compounds.by.name <- function(name){
    .find("compound", name)
}

search.glycans.by.name <- function(name){
    .find("glycan", name)
}

search.compounds.by.composition <- function(composition){
    .find("compound", sprintf("%s/formula", composition))
}

search.glycans.by.composition <- function(composition){
    .find("glycan", composition)
}

search.compounds.by.mass <- function(mass){
    ## range argument goes away, because docs say:
    ## "The exact mass (or molecular weight) is checked by
    ## rounding off to the same decimal place as the query data."
    ## example: /find/compound/174.05/exact_mass  
    ## for 174.045 =< exact mass < 174.055
    .find("compound", sprintf("%s/exact_mass/", mass))
}

search.glycans.by.mass <- function(mass, range){
    unlist(.SOAP(KEGGserver, "search_glycans_by_mass",
                 .soapArgs=list('mass' = mass, 'range' = range),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}


search.compounds.by.subcomp <- function(mol, offset, limit){
    res <- .SOAP(KEGGserver, "search_compounds_by_subcomp",
                 .soapArgs=list('mol' = mol, 'offset' = offset, 
                   'limit' = limit),
                 action = KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    extractFromNamedListObjs(res)
}
## library(KEGGSOAP); mol <- bget("-f m cpd:C00111")
## c4 <- search.compounds.by.subcomp(mol, 1, 5)


search.glycans.by.kcam <- function(kcf, program, option, offset, limit){
  res <- .SOAP(KEGGserver, "search_glycans_by_kcam",
               .soapArgs=list('kcf'=kcf, 'program'=program, 'option'=option,
                 'offset'=offset, 'limit'=limit),
               action = KEGGaction, xmlns = KEGGxmlns,
               nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  extractFromNamedListObjs(res)
}
## library(KEGGSOAP); kcf <- bget("-f k gl:G12922")
## g4 <- search.glycans.by.kcam(kcf, "gapped", "local", 1, 5)

bget <- function(bget.command){
    unlist(.SOAP(KEGGserver, "bget",
                 .soapArgs=list('str' = bget.command),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

## bconv() is handling the issue with only one element all by itself (older).
## I need to handle some formatting going and coming for this to work right.
## The input needs to be a vector (which I must parse into a space separated
## string).  and the output needs to be put into another vector of results (a
## named vector)
bconv <- function(id.list){
    if(length(id.list>1)){id.list=paste(id.list,collapse=" ")}
    res <- unlist(.SOAP(KEGGserver, "bconv",
                 .soapArgs=list('str' = id.list),
                 action = KEGGaction, xmlns = KEGGxmlns,
                        nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
    res <- strsplit(res, split = "\n")[[1]]
    res <- strsplit(res, split = "\t")
    unlist(lapply(res,function(e){
      v <- e[2]
      names(v) <-e[1]
      v
    }))
}



## Support was requested for kegg orthology "ko" numbers.

get.ko.by.gene <- function(genes.id) {
    .get.x("ko", genes.id)
}



## next three need the return value split() into a named vector.

## TODO: compound object returned.  But it should be a simple vector...
get.ko.by.ko.class <- function(ko.class.id) {
      res<-.SOAP(KEGGserver, "get_ko_by_ko_class",
                      .soapArgs=list('ko_class_id' = ko.class.id),
                      action=KEGGaction, xmlns = KEGGxmlns,
                      nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
    extractFromDefinitions(res)
}
## library(KEGGSOAP); ko <- get.ko.by.ko.class("00524")


get.genes.by.ko.class <- function(ko.class.id, org, offset, limit){
    res<-unlist(.SOAP(KEGGserver, "get_genes_by_ko_class",
                      .soapArgs=list('ko_class_id'=ko.class.id, org=org,
                        offset=offset, limit=limit),
                      action = KEGGaction, xmlns = KEGGxmlns,
                      nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
    ans <- lapply(res, function(x){x@entry_id})
    names(ans) <- lapply(res, function(x){x@definition})
    unlist(ans)
}

get.genes.by.ko  <- function(ko.id, org){
    ## how to filter by org?
    ## also, need long name (annotation)
    res<-unlist(.SOAP(KEGGserver, "get_genes_by_ko",
                      .soapArgs=list('ko_id'=ko.id, org=org),
                      action = KEGGaction, xmlns = KEGGxmlns,
                      nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
    ans <- lapply(res, function(x){x@entry_id})
    names(ans) <- lapply(res, function(x){x@definition})
    unlist(ans)
}
                      
get.kos.by.pathway <- function(pathway.id) {
    ## example returns nothing in rest and python,
    ## but returns 36 kos in KEGGSOAP...???
    .get.x("ko", pathway.id)
}

get.pathways.by.kos <- function(ko.id.list, org){
    ## how to filter by org?
    ## http://rest.kegg.jp/link/pathway/ko%3AK00016+ko%3AK00382
    ## does not return the same thing as KEGGSOAP example
    unlist(.SOAP(KEGGserver, "get_pathways_by_kos",
                 .soapArgs=list('ko_id_list'=ko.id.list, org=org),
                 action = KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

