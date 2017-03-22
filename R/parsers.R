
.matrixParser <- function(txt, ncol)
{
    lines <- strsplit(txt, "\n")[[1]]
    split <- strsplit(lines, "\t")
    u <- unlist(split)
    matrix(u, ncol=ncol, byrow=TRUE)
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

.get_parser_NAME <- function(entry)
{
    ret <- list()
    for (value in names(entry))
    {
        ret[[value]] <- gsub("^;|;$", "", entry[[value]])
    }
    ret
}

.get_parser_ENTRY <- function(entry)
{
    segs <- strsplit(unlist(entry[[1]]), "   +")[[1]]
    ret <- c(segs[1])
    names(ret) <- segs[2]
    ret
}


.get_parser_REFERENCE <- function(refs)
{
    ret <- list()
    thisref <- list()
    for (i in 1:length(refs)) {
    #sapply(refs, function(item) {
        item <- refs[[i]]
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
    #})
    }
    ret <- c(ret, list(thisref))
    ret
}


.get_parser_key_value <- function(entry)
{
    content <- c()
    names <- c()
    lines <- unlist(strsplit(unname(unlist(entry)), "\n", fixed=TRUE))
    for (line in lines)
    {
        tmp <- strsplit(line, "  ", fixed=TRUE)[[1]]
        key <- tmp[1]
        value <- paste(tmp[2:length(tmp)], collapse="  ")
        if (is.na(value))
            value <- ""
        content <- c(content, .strip(value))
        names <- c(names, .strip(key))
    }
    names(content) <- names
    content
}

.get_parser_list <- function(entry)
{
    unname(unlist(strsplit(unlist(entry), " {2,}")))
}

.get_parser_list_or_key_value <- function(entry)
{
    x <- unlist(entry)
    if (any(grepl(" {2,}", x)))
        .get_parser_key_value(entry)
    else
        .get_parser_list(entry)
##        unlist(unname(sapply(entry, strsplit, " ")))
}


.get_parser_biostring <- function(entry, type)
{
    ntseq <- unname(unlist(entry))
    tmp <- ntseq[2:length(ntseq)]
    seq <- paste(tmp, collapse="")
    if (type=="AAStringSet")
        AAStringSet(seq)
    else if (type == "DNAStringSet")
        DNAStringSet(seq)
}


.flatFileParser <- function(txt)
{
    entry <- list()
    refs <- list()
    allEntries <- c()
    last_field <- NULL
    lines <- strsplit(.strip(txt), "\n", fixed=TRUE)[[1]]
    ffrec <- flatFileRecordGen()
    for (line in lines)
    {
        if (line == "///")
        {
            ffrec$flush()
            for (name in ffrec$names())
            {
                item <- ffrec$get(name)
                if (name == "ENTRY")
                    ffrec$set("ENTRY", .get_parser_ENTRY(item))
                if (name %in% c("ENZYME", "MARKER", "ALL_REAC",
                    "RELATEDPAIR", "DBLINKS", "DRUG", "GENE"))
                    ffrec$set(name, .get_parser_list(item))
                if (name %in% c("PATHWAY", "ORTHOLOGY", "PATHWAY_MAP", "MODULE",
                    "DISEASE", "REL_PATHWAY", "COMPOUND",
                    "REACTION", "ORGANISM"))
                    {
                        ffrec$set(name, .get_parser_key_value(item))
                    }
                if (name %in% c("REACTION"))
                {
                    ffrec$set(name, .get_parser_list_or_key_value(item))
                }
                item <- ffrec$get(name)
                if(length(item) == 1 && "list" %in% class(item))
                {
                    item <- unlist(item)
                    item <- unname(item)
                    ffrec$set(name, item)
                }
            }
            if ("NTSEQ" %in% ffrec$names())
            {
                ffrec$set("NTSEQ",
                    .get_parser_biostring(ffrec$get("NTSEQ"), "DNAStringSet"))
            }
            if ("AASEQ" %in% ffrec$names())
            {
                ffrec$set("AASEQ",
                    .get_parser_biostring(ffrec$get("AASEQ"), "AAStringSet"))
            }

            ## dreaded copy-and-append pattern
            allEntries <- c(allEntries, list(ffrec$getFields()))
            ffrec <- flatFileRecordGen()
        } else {
            subfield <- NULL
            tmp <- strsplit(line, "", fixed=TRUE)[[1]]
            fs <- tmp[1:12]
            fs <- fs[!is.na(fs)]
            first12 <- .strip(paste(fs, collapse=""))
            if(is.na(tmp[13]))
                value <- ""
            else
                value <- .rstrip(paste(tmp[13:length(tmp)], collapse=""))
            if (!grepl("^ ", line))
            {
                field <- strsplit(line, " ", fixed=TRUE)[[1]][1]
                ffrec$setField(field)
            } else {
                if (first12 != "")
                {
                    subfield <- first12
                    ffrec$setSubfield(first12)
                }
            }
            ffrec$setBody(value)
        }
    }
    allEntries
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


flatFileRecordGen <- setRefClass("KEGGFlatFileRecord", 
    fields=list("fields"="list",
        lastField="character",
        lastSubfield="character",
        lastReference="list",
        references="list"),
    methods=list(
        initialize=function()
        {
            .self$fields <- list()
            .self$references <- list()
            .self$lastField <- character(0)
            .self$lastSubfield <- character(0)
            .self$lastReference <- list()
        },
        setField=function(field)
        {
            .self$flush()
            .self$lastField <- field
            .self$lastSubfield <- character(0)
            .self
        },
        setSubfield=function(subfield)
        {
            .self$lastSubfield <- subfield
            .self
        },
        setBody=function(body)
        {
            if (.self$lastField == "REFERENCE")
            {
                if(length(.self$lastSubfield))
                {
                    if(is.null(.self$lastReference[[.self$lastSubfield]]))
                        .self$lastReference[[.self$lastSubfield]] <- c()
                    .self$lastReference[[.self$lastSubfield]] <- c(
                        .self$lastReference[[.self$lastSubfield]],
                        body)
                } else {
                    if(is.null(.self$lastReference[[.self$lastField]]))
                        .self$lastReference[[.self$lastField]] <- c()
                    .self$lastReference[[.self$lastField]] <- c(
                        .self$lastReference[[.self$lastField]],
                        body)
                }
            } else{
                if (is.null(.self$fields[[.self$lastField]]))
                    .self$fields[[.self$lastField]] <- list()

                if(length(.self$lastSubfield))
                {
                    if(is.null(.self$fields[[.self$lastField]][[.self$lastSubfield]]))
                        .self$fields[[.self$lastField]][[.self$lastSubfield]] <- c()
                    .self$fields[[.self$lastField]][[.self$lastSubfield]] <- c(
                        .self$fields[[.self$lastField]][[.self$lastSubfield]],
                        body
                    )
                } else {
                    if (is.null(.self$fields[[.self$lastField]][[.self$lastField]]))
                        .self$fields[[.self$lastField]][[.self$lastField]] <- c()
                    .self$fields[[.self$lastField]][[.self$lastField]] <- c(
                        .self$fields[[.self$lastField]][[.self$lastField]], body)
                }
            }
            .self
        },
        flush = function()
        {
            .self$fields[["///"]] <- NULL
            if (length(.self$lastReference))
            {
                .self$references[[length(.self$references)+1]] <- .self$lastReference
                .self$lastReference <- list()
            }
            .self
        },
        names = function()
        {
            nms <- base::names(.self$fields)
            if (length(.self$references))
                nms <-c(nms, "REFERENCE")
            nms
        },
        get = function(name)
        {
            if (name == "REFERENCE")
                return(.self$references)
            return(.self$fields[[name]])
        },
        set = function(name, value)
        {
            .self$fields[[name]] <- value
            .self
        }, getFields = function()
        {
            f <- .self$fields
            if (length(.self$references))
                f[["REFERENCE"]] <- .self$references
            f
        }
    )
)
