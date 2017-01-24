#' Derive a HMM from training sequences.
#'
#' \code{derive.HMM} calculates the maximum likelihood hidden Markov model from
#' a list of training sequences, each a vector of residues named according
#' the state from which they were emitted.
#'
#' @param x a list of named character vectors representing residue emissions
#' from the model. The 'names' attribute should represent the hidden state
#' that each residue was emitted from.
#' @param residues either NULL (default; emitted residues are automatically
#' detected from the input sequences), or a case sensitive character vector specifying the
#' residue alphabet (e.g. A, C, G, T for DNA).
#' Note that the former option can be slow for large sequence lists;
#' therefore specifying the residue alphabet can increase speed in these cases.
#' Also note that the default setting \code{residues = NULL} will not
#' detect rare residues that are not present in the input sequences, and thus will
#' not assign them emission probabilities.
#' @param states either NULL (default; states are automatically detected from the 'names'
#' attributes of the input sequences), or a case sensitive character
#' vector matching the hidden states (these should match or be a superset of
#' the unique states provided in the 'names' attributes in the input sequences).
#' @param pseudocounts used to account for the possible absence of certain transition
#' and/or emission types in the input sequences.
#' either \code{'Laplace'} (adds one of each possible transition and emission type to the
#' training dataset; default), \code{'none'}, or a two-element list containing a matrix of
#' transition pseudocounts as its first element and a matrix of emission pseudocounts
#' as its second. If a list is supplied both matrices must have row and column names
#' according to the residues (column names of emission matrix) and states
#' (row and column names of the transition matrix and row names of the emission matrix).
#' The first row and column of the transition matrix must be 'Begin'.
#' @param modelend logical indicating whether transitions to the 'end' state should be
#' modeled. Defaults to FALSE.
#' @return an object of class \code{"HMM"}
#' @export
#'
derive.HMM <- function(x, seqweights = NULL, residues = NULL, states = NULL, modelend = FALSE,
                      pseudocounts = "background", logspace = FALSE, k = 1){
  if(!(is.list(x))) stop("x must be a list of named vectors")
  # x is a list of named character vectors
  # includes start and or end states?
  DNA <- is.DNA(x)
  AA <- is.AA(x)
  namesok <- all(sapply(x, function(y) !is.null(names(y)) | length(y) == 0))
  if(!(namesok)) stop("x must be a list of named vectors")
  residues <- alphadetect(x, residues = residues)
  if(is.null(states)) states <- unique(unlist(lapply(x, names)))
  if(states[1] != "Begin") states <- c("Begin", states)
  nres <- length(residues)
  nstates <- length(states)
  n <- length(x)
  if(is.null(seqweights)){
    seqweights <- rep(1, n)
  }else{
    if(round(sum(seqweights), 2) != n){
      if(round(sum(seqweights), 2) == 1){
        seqweights <- seqweights * n
      }else stop("invalid seqweights argument")
    }
  }
  # code states as integers
  indices <- setNames(seq_along(states) - 1, states)
  #indices <- 0:(nstates - 1)
  #names(indices) <- states
  pathscoded <- lapply(x, function(e) indices[c("Begin", names(e), "Begin")])
  Acounts <- transitioncount(pathscoded[[1]], arity = nstates) * seqweights[1]
  if(n > 1){
    for(i in 2:n){
      Acounts <- Acounts + transitioncount(pathscoded[[i]], arity = nstates) * seqweights[i]
    }
  }
  indices <- setNames(0:(nstates - 2), states[-1])
  statescoded <- lapply(x, function(e) indices[names(e)])
  if(DNA){
    #rescoded <- lapply(x, DNA2quaternary)
    rescoded <- lapply(x, encode.DNA, arity = 4)
  }else if(AA){
    #rescoded <- lapply(x, AA2vigesimal)
    rescoded <- lapply(x, encode.AA, arity = 20)
  }else{
    indices <- setNames(0:(nres - 1), residues)
    rescoded <- lapply(x, function(e) indices[e])
  }
  Ecounts <- emissioncount(statescoded[[1]], nstates - 1, rescoded[[1]], nres) * seqweights[1]
  if(n > 1) for(i in 2:n){
    Ecountsi <- emissioncount(statescoded[[i]], nstates - 1, rescoded[[i]], nres) * seqweights[i]
    Ecounts <- Ecounts + Ecountsi
  }
  # add pseudocounts
  if(identical(pseudocounts, "background")){
    Apseudocounts <- Acounts + 1
    Apseudocounts <- Apseudocounts/sum(Apseudocounts) * (nstates^2 - if(modelend) 0 else nstates)
    Acounts <- Acounts + Apseudocounts
    Epseudocounts <- Ecounts + 1
    Epseudocounts <- Epseudocounts/sum(Epseudocounts) * (nstates - 1) * nres
    Ecounts <- Ecounts + Epseudocounts
  } else if(identical(pseudocounts, "Laplace")){
    Acounts <- Acounts + 1
    Ecounts <- Ecounts + 1
  }else if(is.list(pseudocounts)){
    stopifnot(length(pseudocounts == 2))
    Acounts <- Acounts + pseudocounts[[1]]
    Ecounts <- Ecounts + pseudocounts[[2]]
  }else if(!identical(pseudocounts, "none")) stop("invalid pseudocounts argument")
  if(!modelend) Acounts[, 1] <- 0
  # calculate transition and emission matrices
  A <- Acounts/apply(Acounts, 1, sum) #rows must sum to 1
  dimnames(A) <- list(from = states, to = states)
  E <- Ecounts/apply(Ecounts, 1, sum)
  dimnames(E) <- list(state = states[-1], residue = residues)
  if(logspace){
    A <- log(A)
    E <- log(E)
  }
  res <- structure(list(A = A, E = E), class = "HMM")
  return(res)
}
