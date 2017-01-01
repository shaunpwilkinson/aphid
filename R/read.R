#' Read dendrogram from parenthetic text
#'
#' Reads a file or character string in Newick (New Hampshire) format into
#' an object of class \code{dendrogram}.
#'
#' @param file the name of the file to read the data from.
#' @param text character string: if a text argument is provided instead of a
#' file path then data are read from the value of text via a text connection.
#' @param edges a logical value indicating whether edge weights
#' provided in the Newick string should be retained (defaults to TRUE).
#' @param ... further arguments to be passed to \code{scan}.
#' @details discards comments enclosed in square brackets
#'
#' @return an object of class \code{"dendrogram"}.
#'
#' @seealso \code{\link{write.dendrogram}} to write an object of
#' class \code{"dendrogram"} to a text string.
#'
#' @examples mynewick <- "(A:1,(B:0.5,C:0.2):0.6);"
#' mydendrogram <- read.dendrogram(text = mynewick)
#' plot(mydendrogram)
#'
#'

read.dendrogram <- function(file = "", text = NULL, edges = TRUE, ...){
  if(!is.null(text)){
    if(!is.character(text))
      stop("Argument 'text' must be of mode character.")
    x <- text
  }else{
    x <- scan(file = file, what = "", sep = "\n", quiet = TRUE, ...)
  }
  if (identical(x, character(0))) {
    warning("Empty character string.")
    return(NULL)
  }
  # collapse vector to a single string if necessary
  x <- paste0(x, collapse = "")
  # enclose entire string in brackets (excluding ;)
  xsplit <- strsplit(x, split = "")[[1]]
  if(xsplit[length(xsplit) - 1] != ")"){
    xsplit <- c("(", xsplit[-length(xsplit)], ");")
    x <- paste0(xsplit, collapse = "")
  }
  has.comments <- grepl("\\[", x) | grepl("\\]", x)
  if(has.comments){
    opens <- which(xsplit == "[")
    closes <- which(xsplit == "]")
    if(length(opens) != length(closes)) stop("Invalid metacharacters in Newick string")
    comments <- unlist(mapply(":", opens, closes))
    xsplit <- xsplit[-comments]
    x <- paste0(xsplit, collapse = "")
  }
  if(!edges){
    x <- gsub(":([-0-9Ee.]+)", "", x)
    has.edges <- FALSE
  }else{
    has.edges <- grepl(":", x)
  }
  has.unmatched.singlequotes <- grepl("''", x)
  if(has.unmatched.singlequotes) x <- gsub("''", "singlequote", x)
  # rectified later with fixnames function

  fun1 <- function(s, has.edges){ # a string, applied to odds (not enclosed in single quotes)
    # Underscore characters outside unquoted labels are converted to blanks.
    #s <- gsub("_", "", s)
    s <- gsub(" ", "", s)
    # blank leaves are renamed, rectified later with fixnames function
    s <- gsub("\\( *,", "\\(unnamedleaf,", s)
    while(grepl(", *,", s)) s <- gsub(", *,", ",unnamedleaf,", s)
    s <- gsub(", *\\)", ",unnamedleaf\\)", s)
    # remove inner nodes (for now)
    s <- gsub("\\)[^,:;\\(\\)]+([,:;\\(\\)])", "\\)\\1", s)

    if(has.edges){
      s <- gsub("([\\(,])([^\\(\\),]+):", "\\1'\\2':", s)
      s <- gsub(";", ":1", s)
    }else{
      s <- gsub("([\\(,])([^\\(\\),]+)", "\\1\'\\2':1", s)
      s <- gsub("\\)","\\):1", s)
      s <- gsub(";", "", s)
    }
    return(s)
  }
  has.singlequotes <- grepl("'", x)
  if(has.singlequotes){
    tmp <- strsplit(x, split = "'")[[1]]
    evens <- seq(from = 2, to = length(tmp), by = 2)
    odds <- seq(from = 1, to = length(tmp), by = 2) # not names in single quotes
    tmp[odds] <- unname(sapply(tmp[odds], fun1, has.edges = has.edges))
    tmp[evens] <- unname(sapply(tmp[evens], function(s) paste0("'", s, "'", if(has.edges) "" else ":1")))
    tmp <- paste0(tmp, collapse = "")
  }else{
    tmp <- fun1(x, has.edges = has.edges)
  }
  tmp2 <- strsplit(tmp, split = "'")[[1]]
  evens <- seq(from = 2, to = length(tmp2), by = 2)
  leafnames <- tmp2[evens]
  has.unnamed.leaves <- any(grepl("unnamedleaf", leafnames))
  newleafnames <- paste0("L", seq(100001, 100000 + length(leafnames)))
  tmp2[evens] <- newleafnames
  tmp2 <- paste0(tmp2, collapse = "")

  tmp3 <- tmp2
  tree <- lapply(leafnames, function(e) e)
  names(tree) <- newleafnames
  innernodecount <- 100001
  while(grepl("\\([-LI0123456789Ee,:.]+\\)", tmp3)){
    # I is inner node, L is leaf, max 900000 leaves
    tojoin <- gsub(".*\\(([-LI0123456789Ee,:.]+)\\).*", "\\1", tmp3)
    tojoin2 <- strsplit(tojoin, split = ",")[[1]]
    ntojoin2 <- length(tojoin2)
    whichntj <- integer(ntojoin2)
    newnode <- vector(mode = "list", length = ntojoin2)
    for(s in seq_along(tojoin2)){
      nameedge <- strsplit(tojoin2[s], split = ":")[[1]]
      whichntj[s] <- match(nameedge[1], names(tree))
      attr(tree[[whichntj[s]]], "edge") <- as.numeric(nameedge[2])
      newnode[[s]] <- tree[[whichntj[s]]]
    }
    newnodename <- paste0("I", innernodecount)
    tree[[newnodename]] <- newnode
    tree <- tree[-whichntj]
    tmp3 <- gsub(paste0("\\(", tojoin, "\\)"), newnodename, tmp3)
    innernodecount <- innernodecount + 1
  }
  if(!grepl("^[LI][0-9]{6}:([-0-9Ee.]+)$", tmp3)) warning("Incomplete tree parse")
  #if(length(tree) == 1 & is.list(tree[[1]])) tree <- tree[[1]]
  tree <- tree[[1]]
  attr(tree, "edge") <- 0
  # convert nested list to dendrogram object by setting attributes recursvely

  setnodeattr <- function(x, leafnames){ # x is a nested list with 'edge' attributes, leafnmes is a character vector
    if(is.list(x)){
      cladesizes <- sapply(x, function(y) length(unlist(y)))
      attr(x, "members") <- sum(cladesizes)
      attr(x, "midpoint") <- if(length(cladesizes) > 1){
        ((cladesizes[1] - 1)/2 + (cladesizes[1] + (cladesizes[2] - 1)/2))/2
      }else 0
      if(is.null(attr(x, "height"))) attr(x, "height") <- 0
      setleafattr <- function(y, leafnames){
        attr(y, "height") <- attr(x, "height") - attr(y, "edge")
        attr(y, "edge") <- NULL
        if(!(is.list(y))){
          attr(y, "label") <- y[1]
          attr(y, "leaf") <- TRUE
          attr(y, "members") <- 1
          y[] <- match(y, leafnames)
          mode(y) <- "integer"
        }
        y
      }
      x[] <- lapply(x, setleafattr, leafnames = leafnames)
      x[] <- lapply(x, setnodeattr, leafnames = leafnames)
    }
    x
  }
  res <- setnodeattr(tree, leafnames = leafnames)
  if(!is.list(res)){
    attr(res, "height") <- attr(res, "edge")
    attr(res, "label") <- res[1]
    attr(res, "members") <- 1
    attr(res, "leaf") <- TRUE
  }
  attr(res, "edge") <- NULL
  # now convert to dendrogram
  attr(res, "class") <- "dendrogram"
  min.height <- min(unlist(dendrapply(res, attr, "height")))
  reposition <- function(y, min.height){ # y is a dendrogram
    attr(y, "height") <- attr(y, "height") - min.height
    y
  }
  res <- dendrapply(res, reposition, min.height = min.height)
  fixnames <- function(y){
    if(!(is.list(y))){
      attr(y, "label") <- gsub("unnamedleaf", "", attr(y, "label"))
      attr(y, "label") <- gsub("singlequote", "'", attr(y, "label"))
    }
    y
  }
  if(has.unnamed.leaves | has.unmatched.singlequotes) res <- dendrapply(res, fixnames)
  ultrametricize <- function(y){
    if(is.leaf(y)) attr(y, "height") <- 0
    y
  }
  if(!has.edges) res <- dendrapply(res, ultrametricize)
  return(res)
}





read.PHMM <- function(file = "", ...){
  x <- scan(file = file, what = "", sep = "\n", quiet = TRUE, ... = ...)
  # make this part a new fun and lapply in list case
  if(identical(x, character(0))){
    warning("Empty character string.")
    return(NULL)
  }
  modstart <- grep("^HMMER3/f", x)[1]
  if(!(length(modstart) > 0)) stop("only HMMER3 save file format
                                   version f supported at this stage")
  modend <- modstart + grep("^//", x[modstart:length(x)])[1] - 1

  x <- x[modstart:modend]
  out <- list()

  out$name <- gsub("^NAME +", "", x[2]) #Mandatory

  i <- grep("^ACC ", x[1:30]) # Optional
  if(length(i) > 0) out$accession <- gsub("ACC +", "", x[i[1]])

  i <- grep("^DESC ", x[1:30]) # Optional
  if(length(i) > 0) out$description <- gsub("DESC +", "", x[i[1]])

  out$size <- as.integer(gsub("^LENG +", "", x[grep("LENG ", x[1:30])[1]])) #Mandatory

  i <- grep("^MAXL ", x[1:30]) # Optional
  if(length(i) > 0) out$maxlength <- as.integer(gsub("MAXL +", "", x[i[1]]))

  out$alphabet <- trimws(gsub("ALPH", "", x[grep("^ALPH ", x[1:30])[1]]))  #Mandatory

  i <- grep("^RF ", x[1:30]) #Optional
  hasref <- if(length(i) > 0) trimws(toupper(gsub("RF", "", x[i[1]]))) == "YES" else F

  i <- grep("^MM ", x[1:30]) #Optional
  hasmm <- if(length(i) > 0) trimws(toupper(gsub("MM", "", x[i[1]]))) == "YES" else F

  i <- grep("^CONS ", x[1:30]) #Optional
  hascons <- if(length(i) > 0) trimws(toupper(gsub("CONS", "", x[i[1]]))) == "YES" else F

  i <- grep("^CS ", x[1:30]) #Optional
  hascs <- if(length(i) > 0) trimws(toupper(gsub("CS", "", x[i[1]]))) == "YES" else F

  i <- grep("^MAP ", x[1:30]) #Optional
  hasmap <- if(length(i) > 0) trimws(toupper(gsub("MAP", "", x[i[1]]))) == "YES" else F

  i <- grep("^DATE ", x[1:30]) #Optional
  if(length(i) > 0) out$date <- gsub("DATE +", "", x[i[1]]) # leave as string

  i <- grep("^NSEQ ", x[1:30])#Optional
  if(length(i) > 0) out$nseq <- as.integer(gsub("NSEQ +", "", x[i]))

  i <- grep("^EFFN ", x[1:30])#Optional
  if(length(i) > 0) out$effn <- as.numeric(gsub("EFFN +", "", x[i]))

  i <- grep("^CKSUM ", x[1:30]) #Optional
  if(length(i) > 0) out$checksum <- as.integer(gsub("CKSUM +", "", x[i]))

  i <- grep("^HMM ", x[1:30])[1] #Mandatory, note use ^ for start, $ for end
  out$residues <- unlist(strsplit(gsub("HMM +", "", x[i]), split = " +"))
  transitions <- unlist(strsplit(x[i + 1], split = " +"))
  transitions <- transitions[transitions != ""]
  trans2 <- toupper(gsub("->", "", transitions))
  out$DI <- "d->i" %in% transitions
  out$ID <- "i->d" %in% transitions
  whichqe <- i + 2
  if(grepl("COMPO ", x[whichqe])){
    compo <- as.numeric(unlist(strsplit(gsub("COMPO +", "", x[whichqe]), split = " +")))
    compo <- compo[!is.na(compo)]
    names(compo) <- out$residues
    out$compo <- -compo
    whichqe <- whichqe + 1
  }
  qe <- as.numeric(unlist(strsplit(x[whichqe], split = " +")))
  qe <- qe[!is.na(qe)]
  names(qe) <- out$residues
  out$qe <- -qe
  out$A <- matrix(-Inf, nrow = length(trans2), ncol = out$size + 1)
  rownames(out$A) <- trans2 #c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  colnames(out$A) <- 0:out$size
  out$E <- matrix(-Inf, nrow = length(out$residues), ncol = out$size)
  rownames(out$E) <- out$residues
  colnames(out$E) <- 1:out$size
  begintrans <- unlist(strsplit(gsub("^ +", "", x[whichqe + 1]), split = " +"))
  tmp <- if(out$ID) 6 else 5
  out$A[1:tmp, 1] <- -as.numeric(begintrans[1:tmp])
  if(hasmap) out$alignment <- integer(out$size)
  if(hascons) out$consensus <- character(out$size)
  if(hasref) out$reference <- character(out$size)
  if(hasmm) out$mask <- logical(out$size)
  if(hascs) out$construct <- character(out$size)
  wms <- whichqe + 2 # which module start
  for(i in 1:out$size){
    matchprobs <- unlist(strsplit(gsub("^ +", "", x[wms]), split = " +"))
    stopifnot(as.integer(matchprobs[1]) == i)
    counter <- length(out$residues) + 1
    out$E[, i] <- -as.numeric(matchprobs[2:counter])
    counter <- counter + 1
    if(hasmap) out$alignment[i] <- as.integer(matchprobs[counter])
    if(hascons) out$consensus[i] <- matchprobs[counter + 1]
    if(hasref) out$reference[i] <- matchprobs[counter + 2]
    if(hasmm) out$mask[i] <- switch(matchprobs[counter + 3], m = TRUE, FALSE)
    if(hascs) out$construct[i] <- matchprobs[counter + 4]
    if(i < out$size){
      out$A[, i + 1] <- -as.numeric(unlist(strsplit(gsub("^ +", "", x[wms + 2]), split = " +")))
      wms <- wms + 3
    }else{
      stopifnot(grepl("^//", x[wms + 3]))
      endprobs <- unlist(strsplit(gsub("^ +", "", x[wms + 2]), split = " +"))
      out$A[endprobs != "*", i + 1] <- -as.numeric(endprobs[endprobs != "*"])
    }
  }
  if(!out$DI | !out$ID) {
    appdg <- matrix(nrow = 0, ncol = out$size + 1)
    if(!out$DI) appdg <- rbind(appdg, -Inf)
  }
  if(!out$DI){
    tmp <- matrix(rep(-Inf, out$size + 1), nrow = 1)
    rownames(tmp) <- "DI"
    out$A <- rbind(out$A, tmp)
  }
  if(!out$ID){
    tmp <- matrix(rep(-Inf, out$size + 1), nrow = 1)
    rownames(tmp) <- "ID"
    out$A <- rbind(out$A, tmp)
  }
  tmp <- match(c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II"), rownames(out$A))
  out$A <- out$A[tmp, ]
  out$A[is.na(out$A)] <- -Inf
  class(out) <- "PHMM"
  return(out)
}






