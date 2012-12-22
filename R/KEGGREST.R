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

.get.x.by.y <- function(x, y, arg)
{
#  url <- sprintf("%s/link/%s/%s/"
}

.printf <- function(...) message(noquote(sprintf(...)))

.get.x <- function(x, arg)
{
    url <- sprintf("%s/link/%s/%s", .getRootUrl(), x, arg)
    response <- GET(url)
    result <- http_status(response)
    if (!result$category == "success")
        stop(sprintf("Invalid request, server returned %s. (%s)",
            result$message, url))
    .parse(content(response, "text"))
}


.parse <- function(txt)
{
    lines <- strsplit(txt, "\n")
    ret <- character()
    for (line in lines[[1]])
    {
      ret <- c(ret, strsplit(line, "\t")[[1]][2])
    }
    paste(ret, collapse="\n")
    ret
}

get.genes.by.pathway <- function(pathway.id) {
    .get.x("genes", pathway.id)
}


OLDget.genes.by.pathway <- function(pathway.id) {
    unlist(.SOAP(KEGGserver, "get_genes_by_pathway",
                 .soapArgs=list('pathway_id' = pathway.id), action=KEGGaction,
                 xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}


OLDget.genes.by.pathway <- function(pathway.id) {
    unlist(.SOAP(KEGGserver, "get_genes_by_pathway",
                 .soapArgs=list('pathway_id' = pathway.id), action=KEGGaction,
                 xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

get.enzymes.by.pathway <- function(pathway.id){
    unlist(.SOAP(KEGGserver, "get_enzymes_by_pathway",
                 .soapArgs=list('pathway_id' = pathway.id), action=KEGGaction,
                 xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

get.compounds.by.pathway <- function(pathway.id){
    unlist(.SOAP(KEGGserver, "get_compounds_by_pathway",
                 .soapArgs=list('pathway_id' = pathway.id), action=KEGGaction,
                 xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

get.reactions.by.pathway <- function(pathway.id){
    unlist(.SOAP(KEGGserver, "get_reactions_by_pathway",
                        .soapArgs=list('pathway_id' = pathway.id), action=KEGGaction,
                        xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
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
  dbs <- .SOAP(KEGGserver, "list_databases", "",
               action = KEGGaction,
               xmlns = KEGGxmlns,
               nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  extractFromDefinitions(dbs)
}

list.organisms <- function(){
  orgs <- .SOAP(KEGGserver, "list_organisms", "",
                action = KEGGaction, xmlns = KEGGxmlns,
                nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  extractFromDefinitions(orgs)
}

list.pathways <- function(org){
  paths <- .SOAP(KEGGserver, "list_pathways",
                 .soapArgs=list(org = org), action = KEGGaction,
                 xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns))
  extractFromDefinitions(paths)
}



get.genes.by.organism <- function(org, start, max.results){
    unlist(.SOAP(KEGGserver, "get_genes_by_organism",
                        .soapArgs=list(org = org,
          start = start, 'max_results' = max.results),
          action = KEGGaction, xmlns = KEGGxmlns,
                        nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}


## single args with ""
mark.pathway.by.objects <- function(pathway.id, object.id.list){
  f = KEGGIFace@functions[["mark_pathway_by_objects"]]  
    .SOAP(KEGGserver, "mark_pathway_by_objects",
                 .soapArgs = list('pathway_id' = pathway.id,
                                  'object_id_list' = object.id.list),
                 action = KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns),
          .types = environment(f)$.operation@parameters )
}

color.pathway.by.objects <- function(pathway.id, object.id.list,
                                     fg.color.list, bg.color.list){
  f = KEGGIFace@functions[["color_pathway_by_objects"]]
  .SOAP(KEGGserver, "color_pathway_by_objects",
                        .soapArgs=list('pathway_id' = pathway.id,
                             'object_id_list' = object.id.list,
                             'fg_color_list' = fg.color.list,
                             'bg_color_list' = bg.color.list),
                  action = KEGGaction, xmlns = KEGGxmlns,
                  nameSpaces = SOAPNameSpaces(version=KEGGsoapns),
          .types = environment(f)$.operation@parameters )
}

## NOTE: get.pathways.by.genes() just gives the intersection of the pathways
## based on the genes passed in..
get.pathways.by.genes <- function(genes.id.list){
  f = KEGGIFace@functions[["get_pathways_by_genes"]]
  unlist(.SOAP(KEGGserver, "get_pathways_by_genes",
                 .soapArgs=list('genes_id_list' = genes.id.list),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns),
                 .types = environment(f)$.operation@parameters ))
}

get.pathways.by.enzymes <- function(enzyme.id.list){
  f = KEGGIFace@functions[["get_pathways_by_enzymes"]]
  unlist(.SOAP(KEGGserver, "get_pathways_by_enzymes",
                 .soapArgs=list('enzymes_id_list' = enzyme.id.list),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns),
                 .types = environment(f)$.operation@parameters ))
}

get.pathways.by.compounds <- function(compound.id.list){
  f = KEGGIFace@functions[["get_pathways_by_compounds"]]
  unlist(.SOAP(KEGGserver, "get_pathways_by_compounds",
                 .soapArgs=list('compounds_id_list' = compound.id.list),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns),
                 .types = environment(f)$.operation@parameters ))
}

get.pathways.by.reactions <- function(reaction.id.list){
  f = KEGGIFace@functions[["get_pathways_by_reactions"]]
  unlist(.SOAP(KEGGserver, "get_pathways_by_reactions",
                 .soapArgs=list('reactions_id_list' = reaction.id.list),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns),
                 .types = environment(f)$.operation@parameters ))
}


search.compounds.by.name <- function(name){
    unlist(.SOAP(KEGGserver, "search_compounds_by_name",
                 .soapArgs=list('name' = name),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

search.glycans.by.name <- function(name){
    unlist(.SOAP(KEGGserver, "search_glycans_by_name",
                 .soapArgs=list('name' = name),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

search.compounds.by.composition <- function(composition){
    unlist(.SOAP(KEGGserver, "search_compounds_by_composition",
                 .soapArgs=list('composition' = composition),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

search.glycans.by.composition <- function(composition){
    unlist(.SOAP(KEGGserver, "search_glycans_by_composition",
                 .soapArgs=list('composition' = composition),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

search.compounds.by.mass <- function(mass, range){
    unlist(.SOAP(KEGGserver, "search_compounds_by_mass",
                 .soapArgs=list('mass' = mass, 'range' = range),
                 action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
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
    unlist(.SOAP(KEGGserver, "get_ko_by_gene",
                 .soapArgs=list('genes_id' = genes.id),
                 action=KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
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
    res<-unlist(.SOAP(KEGGserver, "get_genes_by_ko",
                      .soapArgs=list('ko_id'=ko.id, org=org),
                      action = KEGGaction, xmlns = KEGGxmlns,
                      nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
    ans <- lapply(res, function(x){x@entry_id})
    names(ans) <- lapply(res, function(x){x@definition})
    unlist(ans)
}
                      
get.kos.by.pathway <- function(pathway.id) {
    unlist(.SOAP(KEGGserver, "get_kos_by_pathway",
                 .soapArgs=list('pathway_id' = pathway.id),
                 action=KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

get.pathways.by.kos <- function(ko.id.list, org){
    unlist(.SOAP(KEGGserver, "get_pathways_by_kos",
                 .soapArgs=list('ko_id_list'=ko.id.list, org=org),
                 action = KEGGaction, xmlns = KEGGxmlns,
                 nameSpaces = SOAPNameSpaces(version=KEGGsoapns)))
}

