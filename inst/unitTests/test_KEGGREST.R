library(KEGGREST)
library(RUnit)

## checker helper
.checkLOL <- function(res){
  all(checkTrue(class(res)=="list"),
      checkTrue(class(res[[1]])=="list"),
      checkTrue(length(res) > 0))
}

.checkCharVec <- function(res){
  all(checkTrue(class(res)=="character"),
      checkTrue(length(res) > 0))
}

.checkPlainText <- function(res) {
    all(checkTrue(class(res)=="character"),
        checkTrue(length(res) == 1))
}

test_keggInfo <- function()
{
    res <- keggInfo("kegg")
    .checkPlainText(res)
    res <- keggInfo("pathway")
    .checkPlainText(res)
    res <- keggInfo("hsa")
    .checkPlainText(res)

}

test_keggList <- function()
{
    res <- keggList("pathway")
    .checkCharVec(res)
    res <- keggList("pathway", "hsa")
    .checkCharVec(res)
    res <- keggList("organism")
    checkTrue("matrix" %in% class(res))
    .checkCharVec(res)
    res <- keggList("hsa")
    .checkCharVec(res)
    res <- keggList("T01001")
    .checkCharVec(res)
    res <- keggList(c("hsa:10458", "ece:Z5100"))
    .checkCharVec(res)
    res <- keggList(c("cpd:C01290","gl:G00092"))
    .checkCharVec(res)
    res <- keggList(c("C01290+G00092"))
    .checkCharVec(res)
}

## The thorough thing to do would be to hit /list/x for each
## x in listDatabases, but that might slam KEGG too hard and
## make them mad. Instead we hit /info. KEGG does not like
## /info/organism for some reason so we will test /list/organism.
test_listDatabases <- function()
{
    dbs <- listDatabases()
    for (db in dbs)
    {
        if (db != "organism")
        {
            res <- keggInfo(db)
            .checkPlainText(res)
        }
    }
    res <- keggList("organism")
    .checkCharVec(res)
}


test_keggFind <- function()
{
    res <- keggFind("genes", c("shiga", "toxin"))
    .checkCharVec(res)
    res <- keggFind("genes", "shiga toxin")
    .checkCharVec(res)
    res <- keggFind("compound", "C7H10O5", "formula")
    .checkCharVec(res)
    res <- keggFind("compound", "O5C7", "formula")
    .checkCharVec(res)
    res <- keggFind("compound", 174.05, "exact_mass")
    .checkCharVec(res)
    res <- keggFind("compound", 300:310, "mol_weight")
    .checkCharVec(res)
}

test_keggGet <- function()
{
    res <- keggGet(c("cpd:C01290", "gl:G00092"))
    .checkLOL(res)
    res <- keggGet(c("C01290", "G00092"))
    .checkLOL(res)
    res <- keggGet(c("hsa:10458", "ece:Z5100"))
    .checkLOL(res)
    res <- keggGet("ec:1.1.1.1")
    .checkLOL(res)
    .checkLOL(res[[1]]$REFERENCE)
    res <- keggGet(c("hsa:10458", "ece:Z5100"), "aaseq")
    checkTrue("AAStringSet" %in% class(res))
    png <- keggGet("hsa05130", "image")
    checkTrue("array" %in% class(png))
}

test_keggConv <- function()
{
    res <- keggConv("eco", "ncbi-geneid")
    .checkCharVec(res)
    res <- keggConv("ncbi-geneid", "eco")
    .checkCharVec(res)
    res <- keggConv("ncbi-gi", c("hsa:10458", "ece:Z5100"))
    .checkCharVec(res)
}

test_keggLink <- function()
{
    res <- keggLink("pathway", "hsa")
    .checkCharVec(res)
    res <- keggLink("hsa", "pathway")
    .checkCharVec(res)
    res <- keggLink("pathway", c("hsa:10458", "ece:Z5100"))
    .checkCharVec(res)
}

test_mark_and_color_pathways_by_objects  <- function(){
  url <- mark.pathway.by.objects("path:eco00260",
                                 c("eco:b0002", "eco:c00263"))
  .checkCharVec(url)
  checkTrue(grep("http://", url)==1)
  url <- color.pathway.by.objects("path:eco00260",
                                  c("eco:b0002", "eco:c00263"),
                                  c("#ff0000", "#00ff00"),
                                  c("#ffff00", "yellow"))
  .checkCharVec(url)
  checkTrue(grep("http://", url)==1)
}
