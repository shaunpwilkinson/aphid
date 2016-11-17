#' Read dendrogram from parenthetic text
#'
#' Reads a file or character string in Newick (New Hampshire) format into
#' an object of class \code{dendrogram}.
#'
#' @param file the name of the file to read the data from.
#' @param text character string: if a text argument is provided instead of a
#' file path then data are read from the value of text via a text connection.
#' @param strip.edges a logical value indicating whether edge weights
#' provided in the Newick string should be ignored.
#' @param ... further arguments to be passed to \code{scan}.
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
read.dendrogram <- function(file = "", text = NULL, strip.edges = FALSE, ...){
  if(!is.null(text)){
    if(!is.character(text))
      stop("argument 'text' must be of mode character.")
    x <- text
  }else{
    x <- scan(file = file, what = "", sep = "\n", quiet = TRUE, ...)
  }
  if (identical(x, character(0))) {
    warning("empty character string.")
    return(NULL)
  }
  x <- paste0(x, collapse = "")
  if(strip.edges) x <- gsub(":([0-9.]+)", "", x)
  newick.has.edgeweights <- grepl(":", x)
  if(newick.has.edgeweights){
    tmp <- gsub("([[:alnum:]]+):", "\\(\"\\1\"):", x)
    tmp <- gsub(";", ":1", tmp)
  }else{
    tmp <- gsub("([[:alnum:]]+)", "\\(\"\\1\")", x)
    tmp <- gsub("(\\))","\\1:1", tmp)
    tmp <- gsub(";", "", tmp)
  }
  tmp <- gsub("\\(", "structure\\(\\(", tmp)
  tmp <- gsub("ture\\(\\(struc", "ture\\(list\\(struc", tmp)
  tmp <- gsub(":([0-9.]+)", ",edge = \\1)", tmp)
  tmp <- eval(parse(text = tmp))
  attr(tmp, 'edge') <- 0
  # convert nested list to dendrogram object by setting attributes recursvely
  leaflabels <- unlist(tmp)
  setAttributes <- function(x){ # x is a nested list with 'edge' attributes
    if(is.list(x)){
      clade_sizes <- sapply(x, function(y) length(unlist(y)))
      attr(x, 'members') <- sum(clade_sizes)
      attr(x, 'midpoint') <- ((clade_sizes[1] - 1)/2 +
                                (clade_sizes[1] + (clade_sizes[2] - 1)/2))/2
      if(is.null(attr(x, 'height'))) attr(x, 'height') <- 0
      # if(!(exists("leaf.values"))) leaf.values <- cbind(unlist(x), 1:length(unlist(x)))
      x[] <- lapply(x, function(y){
        attr(y, 'height') <- attr(x, 'height') - attr(y, 'edge')
        attr(y, 'edge') <- NULL
        if(!(is.list(y))){
          attr(y, "label") <- y[1]
          attr(y, "leaf") <- TRUE
          attr(y, 'members') <- 1
          # y[1] <- leaf.values[which(leaf.values[,1] == y), 2]
          if(!is.null(leaflabels)){
            y[] <- match(y,  leaflabels)
            mode(y) <- 'integer'
          }
        }
        y
      })
      x[] <- lapply(x, setAttributes)
    }
    x
  }
  res <- setAttributes(tmp)
  attr(res, 'edge') <- NULL
  attr(res, 'class') <- 'dendrogram'
  min.height <- min(unlist(dendrapply(res, attr, 'height')))
  res <- dendrapply(res, function(y){
    attr(y, 'height') <- attr(y, 'height') - min.height
    y
  })
  if(!(newick.has.edgeweights)){
    res <- dendrapply(res, function(y){
      if(is.leaf(y)){
        attr(y, 'height') <- 0
      }
      y
    })
  }
  return(res)
}



read.PHMM <- function(file = "", ...){
  x <- scan(file = file, what = "", sep = "\n", quiet = TRUE, ...)
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
  out$name <- gsub("NAME +", "", x[2]) #Mandatory
  whichacc <- grep("^ACC ", x[1:30]) # Optional
  if(length(whichacc) > 0) out$accession <- gsub("ACC +", "", x[whichacc[1]])
  whichdesc <- grep("DESC ", x[1:30]) # Optional
  if(length(whichdesc) > 0) out$description <- gsub("DESC +", "", x[whichdesc[1]])
  whichsize <- grep("LENG ", x[1:30])[1] #Mandatory
  out$size <- as.integer(gsub("LENG +", "", x[whichsize]))
  whichmaxl <- grep("^MAXL ", x[1:30]) # Optional
  if(length(whichmaxl) > 0) out$maxlength <- as.integer(gsub("MAXL +", "", x[whichmaxl[1]]))
  whichalpha <- grep("ALPH ", x[1:30])[1] #Mandatory
  out$alphabet <- gsub("ALPH +", "", x[whichalpha])
  whichref <- grep("^RF ", x[1:30]) #Optional
  if(length(whichref) > 0){
    ref <- gsub("RF +", "", x[whichref[1]])
    out$has.reference <- switch(toupper(ref), YES = TRUE, NO = FALSE, NULL)
  }
  whichmm <- grep("^MM ", x[1:30]) #Optional
  if(length(whichmm) > 0){
    mm <- gsub("MM +", "", x[whichmm[1]])
    out$has.mask <- switch(toupper(mm), YES = TRUE, NO = FALSE, NULL)
  }
  whichcons <- grep("^CONS ", x[1:30]) #Optional
  if(length(whichcons) > 0){
    cons <- gsub("CONS +", "", x[whichcons[1]])
    out$has.consensus <- switch(toupper(cons), YES = TRUE, NO = FALSE, NULL)
  }
  whichcs <- grep("^CS ", x[1:30]) #Optional
  if(length(whichcs) > 0){
    cs <- gsub("CS +", "", x[whichcs[1]])
    out$has.construct <- switch(toupper(cs), YES = TRUE, NO = FALSE, NULL)
  }
  whichmap <- grep("^MAP ", x[1:30]) #Optional
  if(length(whichmap) > 0){
    map. <- gsub("MAP +", "", x[whichmap[1]])
    out$map <- switch(toupper(map.), YES = TRUE, NO = FALSE, NULL)
  }
  whichnseq <- grep("^NSEQ ", x[1:30])#Optional
  if(length(whichnseq) > 0){
    nseq <- gsub("NSEQ +", "", x[whichnseq[1]])
    out$nseq <- as.integer(gsub("NSEQ +", "", x[whichnseq]))
  }
  whicheffn <- grep("^EFFN ", x[1:30])#Optional
  if(length(whicheffn) > 0){
    effn <- gsub("EFFN +", "", x[whicheffn[1]])
    out$effn <- as.integer(gsub("EFFN +", "", x[whicheffn]))
  }
  whichchecksum <- grep("^CKSUM ", x[1:30]) #Optional
  if(length(whichchecksum) > 0){
    effn <- gsub("CKSUM +", "", x[whichchecksum[1]])
    out$checksum <- as.integer(gsub("CKSUM +", "", x[whichchecksum]))
  }
  whichmodstart <- grep("^HMM ", x[1:30])[1] #Mandatory, note use ^ for start, $ for end
  out$residues <- unlist(strsplit(gsub("HMM +", "", x[whichmodstart]), split = " +"))
  transitions <- unlist(strsplit(x[whichmodstart + 1], split = " +"))
  transitions <- transitions[transitions != ""]
  trans2 <- toupper(gsub("->", "", transitions))
  out$DI <- "d->i" %in% transitions
  out$ID <- "i->d" %in% transitions
  whichqe <- whichmodstart + 2
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
  if(out$map) out$alignment <- integer(out$size)
  if(out$has.consensus) out$consensus <- character(out$size)
  if(out$has.reference) out$reference <- character(out$size)
  if(out$has.mask) out$masked <- logical(out$size)
  if(out$has.construct) out$construct <- character(out$size)
  wms <- whichqe + 2 # which module start
  for(i in 1:out$size){
    matchprobs <- unlist(strsplit(gsub("^ +", "", x[wms]), split = " +"))
    stopifnot(as.integer(matchprobs[1]) == i)
    counter <- length(out$residues) + 1
    out$E[, i] <- -as.numeric(matchprobs[2:counter])
    counter <- counter + 1
    if(out$map) out$alignment[i] <- as.integer(matchprobs[counter])
    if(out$has.consensus) out$consensus[i] <- matchprobs[counter + 1]
    if(out$has.reference) out$reference[i] <- matchprobs[counter + 2]
    if(out$has.mask) out$masked[i] <- switch(matchprobs[counter + 3], m = TRUE, FALSE)
    if(out$has.construct) out$construct[i] <- matchprobs[counter + 4]
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






