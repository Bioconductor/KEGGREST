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

        if (field == "REFERENCE")
        {
            refs <- c(refs,
                list(list(field=field, value=value, refField=refField)))
        }
        entry[[field]] <- c(entry[[field]], value)

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


