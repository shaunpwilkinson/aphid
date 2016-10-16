#' Derive a HMM from training sequences.
#'
#' \code{deriveHMM} calculates the maximum likelihood hidden Markov model from
#' a list of training sequences, each a vector of emitted residues named according
#' the state from which they were emitted.
#'
#' @param x a list of named character vectors representing residue emissions
#' from the model. The 'names' rattribute should represent the hidden state
#' that each residue was emitted from.
#' @param residues either \code{'auto'} (default), or a case
#' sensitive character vector matching the alphabet of residues
#' emitted by the model. Specifying the residue alphabet can increase speed
#' for larger training datasets. Note that setting \code{residues = 'auto'} will not
#' detect rare residues that are not present in the training data and thus will
#' not assign them emission probabilities. These should match or be a superset of
#' the unique residues present in the list of training sequences.
#' @param states either \code{'auto'} (default), or a case sensitive character
#' vector matching the hidden states. These should match or be a superset of
#' the unique states provided in the 'names' attributes in
#' the list of training
#' sequences.
#' @param pseudocounts used to account for the possible absence of certain transition
#' and/or emission types in the training dataset.
#' either \code{'Laplace'} (adds one of each possible transition and emission type to the
#' training dataset; default), \code{'none'}, or a two-element list containing a matrix of
#' transition pseudocounts as its first element and a matrix of emission pseudocounts
#' as its second. If a list is supplied both matrices must have row and column names
#' according to the residues (column names of emission matrix) and states
#' (row and column names of the transition matrix and row names of the emission matrix).
#' The first row and column of the transition matrix must be 'BeginEnd'.
#' @param modelend logical indicating whether transitions to the 'end' state should be
#' modeled. Defaults to FALSE.
#'

deriveHMM <- function(x, residues = 'auto', states = 'auto', modelend = FALSE,
                      pseudocounts = "Laplace", logspace = FALSE){
  if(!(is.list(x))) stop("x must be a list of named vectors")
  # x is a list of named character vectors
  # includes start and or end states?
  namesok <- all(sapply(x, function(y) !is.null(names(y))  | length(y) == 0))
  if(!(namesok)) stop("all elements of x must be named vectors")
  #if(!modelend) Apseudocounts[, 1] <- 0
  unlistx <- unlist(x)
  # states are coerced to character:
  if(identical(residues, "auto")) residues <- as.character(sort(unique(unlistx)))
  if(!(length(residues) > 1)) stop("invalid 'residues' argument")
  if(identical(states, "auto")) states <- c("BeginEnd", sort(unique(names(unlistx))))
  if(!(length(states) > 1)) stop("invalid 'states' argument")
  nres <- length(residues)
  nstates <- length(states)
  # check pseudocounts ok
  Apseudocounts <- matrix(nrow = nstates, ncol = nstates)
  Epseudocounts <- matrix(nrow = nstates - 1, ncol = nres)
  dimnames(Apseudocounts) <- list(from = states, to = states)
  dimnames(Epseudocounts) <- list(state = states[-1], residue = residues)
  if(identical(pseudocounts, "Laplace")){
    Apseudocounts[] <- Epseudocounts[] <- 1
  } else if(identical(pseudocounts, "none")){
    AApseudocounts[] <- Epseudocounts[] <- 0
  } else if(is.list(pseudocounts)){
    stopifnot(length(pseudocounts == 2))
    Apseudocounts[] <- pseudocounts[[1]]
    Epseudocounts[] <- pseudocounts[[2]]
  } else stop("invalid 'pseudocounts' argument")
  # code states as integers
  indices <- 0:(nstates - 1)
  names(indices) <- states
  pathscoded <- lapply(x, function(v) indices[c("BeginEnd", names(v), "BeginEnd")])
  # Transition probabilities
  Afun <- function(v, states){
    tuples <- rbind(v[-length(v)], v[-1])
    decs <- apply(tuples, 2, decimal, length(states))
    # note there are nstates^2 possible 2-tuples
    counts <- rep(0, length(states)^2)
    for(i in 1:length(counts)) counts[i] <- sum(decs == i - 1)
    return(counts)
  }
  Acounts <- apply(sapply(pathscoded, Afun, states), 1, sum)
  Acounts <- matrix(Acounts, nrow = nstates, byrow = TRUE)
  Acounts <- Acounts + Apseudocounts
  if(!(modelend)) Acounts[, 1] <- 0
  A <- Acounts/apply(Acounts, 1, sum) #rows must sum to 1
  # Emission probabilities
  Efun <- function(v, states, residues){
    counts <- matrix(0, length(states) - 1, length(residues))
    dimnames(counts) <- list(states[-1], residues)
    for(i in states[-1]){
      for(j in residues){
        counts[i, j] <- sum(names(v) == i & v == j)
      }
    }
    return(counts)
  }
  Ecounts <- apply(sapply(x, Efun, states, residues), 1, sum)
  Ecounts <- matrix(Ecounts, nrow = nstates - 1, byrow = FALSE)
  Ecounts <- Ecounts + Epseudocounts
  E <- Ecounts/apply(Ecounts, 1, sum)
  # compile hidden Markov model object
  if(logspace){
    A <- log(A)
    E <- log(E)
  }
  res <- structure(list(A = A, E = E), class = "HMM")
  return(res)
}
