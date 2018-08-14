#' Import profile hidden Markov models into R.
#'
#' The \code{readPHMM} function parses a HMMER3 text file into R and creates
#'   an object of class \code{"PHMM"}.
#'
#' @param file the name of the file from which to read the model.
#' @param ... further arguments to be passed to \code{"scan"}.
#' @return an object of class \code{"PHMM"}.
#' @details
#'   This function scans a HMMER3/f text file and creates an object of
#'   class \code{"PHMM"} in R. Note that unlike HMMER, the \pkg{aphid}
#'   package does not currently support position-specific background
#'   emission probabilities, and so only a single vector the same length
#'   as the reside alphabet is included as an element of the returned
#'   object. Also the function currently only parses the first profile
#'   HMM encountered in the text file, with subsequent models ignored.
#' @author Shaun Wilkinson
#' @references
#'   Finn RD, Clements J & Eddy SR (2011) HMMER web server: interactive sequence
#'   similarity searching. \emph{Nucleic Acids Research}. \strong{39}, W29-W37.
#'   \url{http://hmmer.org/}.
#'
#'   HMMER: biosequence analysis using profile hidden Markov models.
#'   \url{http://www.hmmer.org}.
#'
#' @seealso \code{\link{writePHMM}} for writing PHMM objects in HMMER3 text format.
#' @examples
#'   ## Derive a profile hidden Markov model from the small globin alignment
#'   data(globins)
#'   x <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
#'   fl <- tempfile()
#'   writePHMM(x, file = fl)
#'   readPHMM(fl)
################################################################################
readPHMM <- function(file = "", ...){
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
  if(length(i) > 0) out$checksum <- gsub("CKSUM +", "", x[i])
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
  if(hasmap) out$map <- integer(out$size)
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
    if(hasmap) out$map[i] <- as.integer(matchprobs[counter])
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
  tmp <- match(c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II"),
               rownames(out$A))
  out$A <- out$A[tmp, ]
  out$A[is.na(out$A)] <- -Inf
  class(out) <- "PHMM"
  return(out)
}
################################################################################
