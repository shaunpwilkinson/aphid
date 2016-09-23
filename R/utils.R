#' Utilities.
#'

decimal <- function(x, from) sum(x * from^rev(seq_along(x) - 1))


whichismax <- function(v){
  ind <- which(v == max(v, na.rm = TRUE))
  if(length(ind) > 1) ind <- sample(ind, 1)
  ind
}


#' Detect residue alphabet.
#'
#' \code{"alphabet"} performs checks on the format of the "residues" argument
#' to be passed to a variety of other functions. Single-element string arguments
#' such as "DNA" and "AA" are
#' converted to their respective alphabet in a character vector format.
#' @param sequences a character matrix or vector, or a list of character matrices
#' and/or vectors
#'
alphabet <- function(sequences, residues = "autodetect", gapchar = "-"){
  if(inherits(sequences, "DNAbin") | identical(residues, "DNA")){
    residues <- c("A", "C", "G", "T")
  } else if(inherits(sequences, "AAbin") | identical(residues, "AA")){
    residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  }
  else if(identical(residues, "autodetect")) {
    residues <- sort(unique(as.vector(unlist(sequences))))
    residues <- residues[residues != gapchar]
  }else{
    residues <- residues[residues != gapchar]
  }
  if(!(length(residues) > 1 & mode(residues) == "character")) {
    stop("invalid residues argument")
  }
  return(residues)
}


tab <- function(v, residues, seqweights = 1){
  if(identical(seqweights, 1)) seqweights <- rep(1, length(v))
  stopifnot(length(seqweights) == length(v) & sum(seqweights) == length(v))
  res <- structure(integer(length(residues)), names = residues)
  for(i in residues) res[i] <- sum(seqweights[v == i], na.rm = TRUE)
  return(res)
}
#-------------------------------------------------------------------------------------------
# x is a DNAbin vector
tabDNA <- function(x, ambiguities = FALSE, seqweights = 1){
  if(identical(seqweights, 1)) seqweights <- rep(1, length(x))
  stopifnot(length(seqweights) == length(x) & sum(seqweights) == length(x))
  res <- structure(numeric(4), names = c("A", "C", "G", "T"))
  res["A"] <- sum(seqweights[x == 136])
  res["C"] <- sum(seqweights[x == 40])
  res["G"] <- sum(seqweights[x == 72])
  res["T"] <- sum(seqweights[x == 24])
  if(ambiguities){
    xambigs <- x != 4 & (x & as.raw(8)) != 8
    if(any(xambigs)){
      truncx <- x[xambigs]
      truncweights <- seqweights[xambigs]
      R <- sum(truncweights[truncx == 192])/2 # A or G
      M <- sum(truncweights[truncx == 160])/2 # A or C
      W <- sum(truncweights[truncx == 144])/2 # A or T
      S <- sum(truncweights[truncx == 96])/2 # G or C
      K <- sum(truncweights[truncx == 80])/2 # G or T
      Y <- sum(truncweights[truncx == 48])/2 # C or T
      V <- sum(truncweights[truncx == 224])/3 # A, G or C
      H <- sum(truncweights[truncx == 176])/3 # A, C or T
      d <- sum(truncweights[truncx == 208])/3 # A, G or T
      B <- sum(truncweights[truncx == 112])/3 # G, C or T
      N <- sum(truncweights[truncx == 240])/4 # A, G, C or T
      res["A"] <- res["A"] + R + M + W + V + H + d + N
      res["C"] <- res["C"] + M + S + Y + V + H + B + N
      res["G"] <- res["G"] + R + S + K + V + d + B + N
      res["T"] <- res["T"] + W + K + Y + H + d + B + N
    }
  }
  return(res)
}


#' Diagnostic model checks.
#'
validate <- function(x){
  if(inherits(x, "HMM")){
    states <- rownames(x$A)
    if(is.null(colnames(x$A)) | is.null(colnames(x$E))){
      message('both A and E must have dimnames attributes')
    } else if(!(identical(states, colnames(x$A)) & identical(states[-1], rownames(x$E)))){
      message("rownames and colnames of the transitions matrix (A), and
              rownames of the emissions matrix (E) must all be non-NULL
              and identical")
    } else if(rownames(x$A)[1] != "BeginEnd") {
      message("transition probability matrix (A) should include a 'BeginEnd'
              state at row 1 and col 1")
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

DNAprob <- function(a, probs){
  # a is a raw byte in Paradis (2007) format
  # probs is a 4-element numeric vector of probabilities for the set {a,c,g,t}
  # returns the weighted average probability
  if(a <= 4){
    stop("Input sequence contains gaps or unknown characters")
  }else{
    if((a & as.raw(55)) == as.raw(0)){ # is purine?
      if(a == 136){
        probs[1] # A
      } else if(a == 72){
        probs[3] # G
      } else{
        mean(probs[c(1, 3)]) # A or G  ### what if probs are logged??
      }
    } else if((a & as.raw(199)) == as.raw(0)){
      if(a == 40){
        probs[2] # C
      } else if(a == 24){
        probs[4] # T
      } else{
        mean(probs[c(2, 4)]) # A or G  ### what if probs are logged??
      }
    } else if(a == 224){
      mean(probs[-4]) # V (A or C or G)
    } else if(a == 176){
      mean(probs[-3]) # H (A or C or T)
    } else if(a == 208){
      mean(probs[-2]) # D (A or G or T)
    } else if(a == 112){
      mean(probs[-1]) # B (C or G or T)
    } else if(a == 240){
      mean(probs) #N
    } else stop("invalid byte")
  }
}

is.ambiguous <- function(a) a != 4 & (a & as.raw(8)) != 8

disambiguate <- function(a, probs = rep(0.25, 4)){
  # a is a raw byte in Paradis (2007) format
  # probs is a 4-element numeric vector of background probabilities for the set {a,c,g,t}
  # returns a sampled base
  if(a <= 4){
    stop("Input sequence contains gaps or unknown characters")
  }else{
    if((a & as.raw(55)) == as.raw(0)){ # is purine?
      if(a != 136 & a != 72){
        sample(as.raw(c(136, 72)), size = 1, prob = probs[c(1, 3)]) # unknown A or G
      }else{
        return(a) #known base A or G
      }
    }else if((a & as.raw(199)) == as.raw(0)){ # is pyrimidine
      if(a != 40 & a != 24){
        sample(as.raw(c(40, 24)), size = 1, prob = probs[c(2, 4)]) # unknown base C or T
      }else{
        return(a) # known base C or T
      }
      # a,c,g,t = 136, 40, 72, 24
    } else if(a == 224){ # V (A or C or G)
      sample(as.raw(c(136, 40, 72)), size = 1, prob = probs[-4])
    } else if(a == 176){
      sample(as.raw(c(136, 40, 24)), size = 1, prob = probs[-3]) # H (A or C or T)
    } else if(a == 208){
      sample(as.raw(c(136, 72, 24)), size = 1, prob = probs[-2]) # D (A or G or T)
    } else if(a == 112){
      sample(as.raw(c(40, 72, 24)), size = 1, prob = probs[-1]) # B (C or G or T)
    } else if(a == 240){
      sample(as.raw(c(136, 40, 72, 24)), size = 1, prob = probs) #N
    } else stop("invalid byte")
  }
}

decimal.DNAbin <- function(x){
  # x is a short DNABIN vector
  # returns a numeric matrix with decimal pointers in col 1 and their probabilities in col 2
  xlen <- length(x)
  out <- matrix(nrow = 1, ncol = 2)
  if(all((x & as.raw(8)) == 8)){

  }
}



