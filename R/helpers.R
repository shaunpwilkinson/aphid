#' Detect residue alphabet.
#'
#' \code{"alphabet"} performs checks on the format of the "residues" argument
#' to be passed to a variety of other functions. Single-element string arguments
#' such as "bases" and "aminos" are
#' converted to their respective alphabet in a character vector format.
#' @param x a character matrix or vector, or a list of character matrices and/or vectors
#'
alphabet <- function(x, residues = "autodetect", gapchar = "-"){
  if(identical(residues, "autodetect")) residues <- sort(unique(as.vector(unlist(x))))
  if(identical(residues, "aminos")) residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  if(identical(residues, "bases")) residues <- c("A", "C", "G", "T")
  if(gapchar %in% residues) residues <- residues[residues != gapchar]
  if(!(length(residues) > 1)) stop("invalid residues argument")
  return(residues)
}


#' Perform diagnostic model checks.
#'
valid <- function(x){
  if(inherits(x, "HMM")){
    states <- rownames(x$A)
    if(is.null(colnames(x$A)) | is.null(colnames(x$E))){
      message('both A and E must have dimnames attributes')
    } else if(!(identical(states, colnames(x$A)) & identical(states[-1], rownames(x$E)))){
      message("rownames and colnames of the transitions matrix (A), and
         rownames of the emissions matrix (E) must all be non-NULL
         and identical")
    } else if(rownames(x$A)[1] != "BeginEnd") {
      message("transition probability matrix (A) should include a 'BeginEnd' state at row 1 and col 1")
    } else {
      return(TRUE)
    }
    return(FALSE)
  } else if(inherits(x, "PHMM")){
    states <- c("M", "I", "D")
    if(!(all(dimnames(x$A)[[1]] == states) & all(dimnames(x$A)[[3]] == states))){
      message("names for dimensions 1 and 3 of transitions
         array (A) must be c('D', 'M', 'I')")
    } else if(!(ncol(x$A) == ncol(x$E) + 1)){
      message("transitions array (A) should include a begin state")
    } else if(!all(names(x$qe) == rownames(x$E))){
      message("residue names for background emissions vector (qe) should be identical
         to rownames of emissions matrix (E)")
    } else if(!(all(colnames(x$qa) == states) & all(rownames(x$qa) == states))){
      message("rownames and colnames of background transitions matrix
         must be c('D', 'M', 'I')")
    } else {
      return(TRUE)
    }
    return(FALSE)
  } else {
    stop("x must be an object of class 'HMM' or 'PHMM'")
  }
}

#' Detect if model parameters are in log space.
#'
logdetect <- function(x){
  if(inherits(x, "HMM")){
    if(all(x$A <= 0) & all(x$E <= 0)){
      return(TRUE)
    } else if(all(x$A >= 0) & all(x$A <= 1) & all(x$E >= 0) & all(x$E <= 1)){
      return(FALSE)
    } else stop("unable to detect if model probabilities are in log space")
  } else if(inherits(x, "PHMM")){
    if(all(x$A <= 0) & all(x$E <= 0) & all(x$qa <= 0) & all(x$qe <= 0)){
      return(TRUE)
    } else if(all(x$A >= 0) & all(x$A <= 1) & all(x$E >= 0) & all(x$E <= 1) &
             all(x$qa >= 0) & all(x$qa <= 1) & all(x$qe >= 0) & all(x$qe <= 1)){
      return(FALSE)
    } else stop("unable to detect if model probabilities are in log space")
  } else stop("x must be an object of class 'HMM' or 'PHMM'")
}


