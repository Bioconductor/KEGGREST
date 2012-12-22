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

test_bconv <- function(){
  res <- bconv("ncbi-geneid:10")
  .checkCharVec(res)
  res <- bconv(c("ncbi-geneid:100008586", "ncbi-geneid:10")) 
  .checkCharVec(res)
}

test_bget <- function(){
  res <- bget("eco:b0002 hin:tRNA-Cys-1")
  .checkCharVec(res)
  res <- bget("-f -n n eco:b0002 hin:tRNA-Cys-1")
  .checkCharVec(res)
  res <- bget("-f -n a eco:b0002")
  .checkCharVec(res)
}

test_getBestNeighbors <- function(){
  res <-bestGenes <- get.best.neighbors.by.gene("eco:b0002",1, 5)
  .checkLOL(res)
}

test_get.genes.by.motifs <- function(){
  res <- get.genes.by.motifs(c("pf:DnaJ", "ps:DNAJ_2"), 1, 10)
  .checkCharVec(res)
}

test_get.genes.by.organism  <- function(){
  res <- get.genes.by.organism("hsa", 1, 10)
  .checkCharVec(res)
}


test_get.genes.by.pathway  <- function(){
  g <- get.genes.by.pathway("path:eco00020")
  .checkCharVec(g)
  e <- get.enzymes.by.pathway("path:eco00020")
  .checkCharVec(e)
  c <- get.compounds.by.pathway("path:eco00020")
  .checkCharVec(c)
  r <- get.reactions.by.pathway("path:eco00020")
  .checkCharVec(r)
}

test_get.ko.by.gene  <- function(){
  k1 <- get.ko.by.gene("eco:b0002")
  .checkCharVec(k1)
  k2 <- get.ko.by.ko.class("00524")
  .checkCharVec(k2)
  g1 <- get.genes.by.ko.class("00903", "hsa" , 1, 100)
  .checkCharVec(g1)
  g2 <- get.genes.by.ko("ko:K12524", "eco")
  .checkCharVec(g2)
  ks <- get.kos.by.pathway("path:hsa00010")
  .checkCharVec(ks)
  p <- get.pathways.by.kos(c("ko:K00016","ko:K00382"), "hsa")
  .checkCharVec(p)
}

test_get.motifs.by.gene  <- function(){
  m <- get.motifs.by.gene("eco:b0002", "pfam")
  .checkLOL(m)
}

test_get.paralogs.by.gene  <- function(){
  p <- get.paralogs.by.gene("eco:b0002", 1, 10)
  .checkLOL(p)
}


test_get.pathways.by.genes  <- function(){
  p1 <- try(get.pathways.by.genes(c("eco:b0077", "eco:b0078")))
  .checkCharVec(p1)
  p2 <- try(get.pathways.by.enzymes("ec:1.3.99.1"))
  .checkCharVec(p2)
  p3 <- try(get.pathways.by.compounds(c("cpd:C00033", "cpd:C00158"))) 
  .checkCharVec(p3)
  p4 <- try(get.pathways.by.reactions(c("rn:R00959",
                                 "rn:R02740", "rn:R00960", "rn:R01786")))
  .checkCharVec(p4)
}


test_list.organisms  <- function(){
  res <- list.organisms()
  .checkCharVec(res)
}

test_  <- function(){
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

test_search.compounds.by.name  <- function(){
  c1 <- search.compounds.by.name("shikimic acid")
  .checkCharVec(c1)
  c2 <- search.compounds.by.composition("C7H10O5")
  .checkCharVec(c2)
  c3 <- search.compounds.by.mass(174.05, 0.1)
  .checkCharVec(c3)
  mol <- bget("-f m cpd:C00111")
  c4 <- search.compounds.by.subcomp(mol, 1, 5)
  .checkLOL(c4)
}

test_search.glycans.by.name  <- function(){
  g1 <- search.glycans.by.name("Paragloboside")
  .checkCharVec(g1)
  g2 <- search.glycans.by.composition("(Man)4 (GalNAc)1")
  .checkCharVec(g2)
  g3 <- search.glycans.by.mass(689.6, 0.1)
  .checkCharVec(g3)
  kcf <- bget("-f k gl:G12922")
  g4 <- search.glycans.by.kcam(kcf, "gapped", "local", 1, 5)
  .checkLOL(g4)
}


