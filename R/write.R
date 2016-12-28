#' Write dendrogram as a Newick string.
#'
#' Writes a dendrogram object to a file or connection in Newick
#' (New Hampshire) format.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param file a character string naming a file or connection to write the output to.
#' If no file path is specified or \code{file = ""} the result will be printed to the console.
#' @param append a logical value indicating whether the output should be appended to the file.
#' If \code{append = FALSE} the contents of the file will be overwritten (the default setting).
#' @param edges a logical value indicating whether edge weights should be
#' included in the output string.
#' @param ... further arguments to be passed to \code{format} to specify the numbering style
#' of the edge weights (assuming edges = TRUE).
#' @seealso \code{\link{read.dendrogram}} to create a \code{"dendrogram"} object from a
#' text file.
#' @examples arrests.hc <- hclust(dist(USArrests[1:6,]), "ave")
#' arrests.den <- as.dendrogram(arrests.hc)
#' write.dendrogram(arrests.den)
#'
write.dendrogram <- function(x, file = "", append = FALSE, edges = TRUE, ...){
  if(!(inherits(x, "dendrogram"))) stop("Input object must be of class 'dendrogram'")
  renameLeaves <- function(y){
    if(is.leaf(y)){
      tmp <- attr(y, "label")
      if(grepl("[^abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]", tmp)){
        tmp <- paste0(c("'", tmp, "'"), collapse = "")
      }
      y[1] <- tmp
    }
    y
  }
  x <- dendrapply(x, renameLeaves)
  xnames <- unlist(x)
  ynames <- paste0("S", seq_along(xnames) + 10000)
  renameLeaves2 <- function(y){
    if(is.leaf(y)){
      y[1] <- ynames[match(y[1], xnames)]
    }
    y
  }
  x <- dendrapply(x, renameLeaves2)
  x <- dendrapply(x, unclass)
  addEdges <- function(y){
    if(is.list(y)){
      y[] <- lapply(y, function(z){
        attr(z, "edge") <- format(attr(y, "height") - attr(z, "height"), ... = ...)#scientific = FALSE)
        z
      })
      attributes(y)[names(attributes(y)) != "edge"] <- NULL
      y[] <- lapply(y, addEdges)
    }else{
      attributes(y)[names(attributes(y)) != "edge"] <- NULL
    }
    y
  }
  x <- addEdges(x)
  attr(x, "edge") <- 0
  tmp <- deparse(x)
  tmp <- paste0(tmp, collapse = "")
  tmp <- gsub(" ", "", tmp)
  tmp <- gsub("edge=\"([-0-9Ee.]+)\"", "edge=\\1", tmp)
  #for(i in seq_along(xnames)) tmp <- gsub(xnames[i], ynames[i], tmp)
  tmp2 <- gsub("list", ";", tmp) # needs to be a special single char - alternative?
  while(grepl("structure", tmp2)){
    tmp2 <- gsub("structure\\(\"([^\"]*)\",edge=([-0-9Ee.]+)\\)", "\\1:\\2", tmp2)
    tmp2 <- gsub(";\\(([^;\\)]*)\\)", "\"openbracket\\1closebracket\"", tmp2)
  }
  tmp3 <- gsub("openbracket", "\\(", tmp2)
  tmp3 <- gsub("closebracket", "\\)", tmp3)
  res <- gsub("(.*):0$", "\\1;", tmp3)
  if(!edges){
    res <- gsub(":[-0-9Ee.]+", "", res)
  }
  for(i in seq_along(xnames)) res <- gsub(ynames[i], xnames[i], res)
  # fun <- function(taxon, edge=1) paste0("(", taxon, ":", edge, ")")
  # tmp <- gsub("structure", "fun", tmp)
  # tmp <- gsub("list", "paste", tmp)
  # res <- eval(parse(text = tmp))
  # res <- gsub(") \\(", ",", res)
  # res <- substr(res, start = 2, stop = nchar(res)-3)
  # res <- paste0(res, ";")
  # if(strip.edges){
  #   res <- gsub(":[-0-9Ee.]+", "", res)
  # }else{
  #   target <- paste0(c("(\\.", rep("[0-9]", dec.places), ")[0-9]+"), collapse = "")
  #   res <- gsub(target, "\\1", res) ### maybe signif places better? format?
  #   if(dec.places == 0) res <- gsub("\\.", "", res)
  # }
  if(file == ""){
    return(res)
  }else{
    cat(res, file = file, append = append, sep = "\n")
  }
}


write.PHMM <- function(x, file = "", append = FALSE, form = "HMMER3", vers = "f"){
  options(digits = 7, scipen = 10)
  stopifnot(form == "HMMER3" & vers == "f")
  if(!(inherits(x, "PHMM"))) stop("Input object must be of class 'PHMM'")
  cat(paste0(form, "/", vers), file = file, sep = "\n", append = append) ##[produced by R::profile?]
  cat(paste0("NAME  ", x$name), file = file, sep = "\n", append = T)
  if(!is.null(x$accession)) cat(paste0("ACC   ", x$accession), file = file, sep = "\n", append = T)
  if(!is.null(x$description)) cat(paste0("DESC  ", x$description), file = file, sep = "\n", append = T)
  cat(paste0("LENG  ", x$size), file = file, sep = "\n", append = T)
  if(!is.null(x$maxlength)) cat(paste0("MAXL  ", x$maxlength), file = file, sep = "\n", append = T)
  cat(paste0("ALPH  ", x$alphabet), file = file, sep = "\n", append = T)

  rfl <- if(is.null(x$reference)) "RF    no" else "RF    yes"
  mml <- if(is.null(x$mask))      "MM    no" else "MM    yes"
  cl <- if(is.null(x$consensus))  "CONS  no" else "CONS  yes"
  csl <- if(is.null(x$construct)) "CS    no" else "CS    yes"
  mpl <- if(is.null(x$alignment)) "MAP   no" else "MAP   yes"
  cat(rfl, mml, cl, csl, mpl, file = file, sep = "\n", append = TRUE)

  if(!is.null(x$date)) cat(paste0("DATE  ", x$date), file = file, append = T, sep = "\n")
  if(!is.null(x$nseq)) cat(paste0("NSEQ  ", x$nseq), file = file, append = T, sep = "\n")
  if(!is.null(x$effn)) cat(paste0("EFFN  ", round(x$effn, 6)), file = file, append = T, sep = "\n")
  if(!is.null(x$checksum)) cat(paste0("CKSUM ", x$checksum), file = file, append = T, sep = "\n")
  ### placeholder for stats, cutoffs, etc
  residues <- rownames(x$E)
  cat(c("HMM          ", paste0(residues, "        "), "\n"), file = file, append = T, sep = "")
  cat("            m->m     m->i     m->d     i->m     i->i     d->m     d->d",
      file = file, append = T, sep = "\n")
  if(!is.null(x$compo)){
    compoline <- formatC(-x$compo, digits = 5, format = "f")
    compoline[trimws(compoline) == "Inf"] <- "      *"
    compoline <- gsub("(.......).+", "\\1", compoline) # in case any logprobs are > 10
    compoline = paste0("  COMPO   ", paste(compoline, collapse = "  "))
    cat(compoline, file = file, sep = "\n", append = TRUE)
  }

  qe <- formatC(-x$qe, digits = 5, format = "f")
  qe[trimws(qe) == "Inf"] <- "      *"
  qe <- gsub("(.......).+", "\\1", qe)
  qeline = paste0("          ", paste(qe, collapse = "  "))

  # qe <- as.character(format(-x$qe, digits = 6, scientific = FALSE))
  # qeline <- paste0("          ", paste(qe, collapse = "  "))
  cat(qeline, file = file, sep = "\n", append = TRUE)

  A <- formatC(-x$A[c(5, 6, 4, 8, 9, 2, 1), ], digits = 6, format = "f")
  A[trimws(A) == "Inf"] <- "      *"
  A <- gsub("(.......).+", "\\1", A)
  dim(A) <- c(7, x$size + 1)
  A[6, 1] <- A[6, x$size + 1] <- "0.00000"

  #E <- as.character(format(-x$E, digits = 6, scientific = FALSE))
  E <- formatC(-x$E, digits = 6, format = "f")
  E[trimws(E) == "Inf"] <- "      *"
  E <- gsub("(.......).+", "\\1", E)
  dim(E) <- c(length(residues), x$size)

  trans0line = paste("         ", paste(A[, 1], collapse = "  "))
  cat(trans0line, file = file, sep = "\n", append = TRUE)
  for(i in 1:x$size){
    emline <- paste0(if(i < 10) "      " else if(i < 100) "     " else "    ", i, "   ")
    emline <- paste0(emline, paste(E[, i], collapse = "  "))
    if(!is.null(x$alignment)){
      tmp <- x$alignment[i]
      emline <- paste0(emline, if(tmp < 10) "      " else if(tmp < 100) "     " else  "    ")
      emline <- paste0(emline, tmp)
    }else{
      emline <- paste0(emline, "      -")
    }
    emline <- paste0(emline, " ", if(!is.null(x$consensus)) x$consensus[i] else "-")
    emline <- paste0(emline, " ", if(!is.null(x$reference)) x$reference[i] else "-")
    emline <- paste0(emline, " ", if(!is.null(x$mask)) x$mask[i] else "-")
    emline <- paste0(emline, " ", if(!is.null(x$construct)) x$construct[i] else "-")
    cat(emline, file = file, sep = "\n", append = TRUE)
    cat(qeline, file = file, sep = "\n", append = TRUE)
    transline <- paste0("          ", paste(A[, i + 1], collapse = "  "))
    cat(transline, file = file, sep = "\n", append = TRUE)
  }
  cat("//", file = file, sep = "", append = TRUE)
}
