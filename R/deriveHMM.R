#' Derive a standard hidden Markov model from a set of sequences.
#'
#' \code{deriveHMM} calculates the maximum likelihood hidden Markov model from
#'   a list of training sequences, each a vector of residues named according
#'   the state from which they were emitted.
#'
#' @param x a list of named character vectors representing emissions
#'   from the model. The 'names' attribute should represent the hidden state
#'   from which each residue was emitted. "DNAbin" and "AAbin" list
#'   objects are also supported for modeling DNA or amino acid sequences.
#' @param seqweights either NULL (all sequences are given
#'   weights of 1) or a numeric vector the same length as \code{x} representing
#'   the sequence weights used to derive the model.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Furthermore, the default setting
#'   \code{residues = NULL} will not detect rare residues that are not present
#'   in the sequences, and thus will not assign them emission probabilities
#'   in the model. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param states either NULL (default; the unique Markov states are
#'   automatically detected from the 'names' attributes of the input
#'   sequences), or a case sensitive character vector specifying the unique
#'   Markov states (or a superset of the unique states) to appear in the
#'   model. The latter option is recommended since it saves computation time
#'   and ensures that all valid Markov states appear in the model,
#'   regardless of their possible absence from the training dataset.
#' @param modelend logical indicating whether transition probabilites
#'   to the end state of the standard hidden Markov model should be
#'   modeled (if applicable). Defaults to FALSE.
#' @param pseudocounts character string, either "background", Laplace"
#'   or "none". Used to account for the possible absence of certain
#'   transition and/or emission types in the input sequences.
#'   If \code{pseudocounts = "background"} (default), pseudocounts
#'   are calculated from the background transition and emission
#'   frequencies in the training dataset.
#'   If \code{pseudocounts = "Laplace"} one of each possible transition
#'   and emission type is added to the training dataset (default).
#'   If \code{pseudocounts = "none"} no pseudocounts are added (not
#'   usually recommended, since low frequency transition/emission types
#'   may be excluded from the model).
#'   Alternatively this argument can be a two-element list containing
#'   a matrix of transition pseudocounts
#'   as its first element and a matrix of emission pseudocounts as its
#'   second. If this option is selected, both matrices must have row and column
#'   names corresponding with the residues (column names of emission matrix)
#'   and states (row and column names of the transition matrix and
#'   row names of the emission matrix). For downstream applications
#'   the first row and column of the transition matrix should be named
#'   "Begin".
#' @param logspace logical indicating whether the emission and transition
#'   probabilities in the returned model should be logged. Defaults to TRUE.
#' @return an object of class \code{"HMM"}.
#' @details
#'   This function creates a standard hidden Markov model (object class:
#'   \code{"HMM"}) using the method described in Durbin et al (1998) chapter
#'   3.3. It assumes the state sequence is known
#'   (as opposed to the \code{\link{train.HMM}} function, which is used
#'   when the state sequence is unknown) and provided as the names attribute(s)
#'   of the input sequences. The output object is a simple list with elements
#'   "A" (transition probability matrix) and "E" (emission probability matrix),
#'   and the "class" attribute "HMM". The emission matrix has the same number
#'   of rows as the number of states, and the same number of columns as the
#'   number of unique symbols that can be emitted (i.e. the residue alphabet).
#'   The number of rows and columns in the transition probability matrix
#'   should be one more the number of states, to include the silent "Begin"
#'   state in the first row and column. Despite its name, this state is
#'   also used when modeling transitions to the (silent)
#'   end state, which are entered in the first column.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @seealso \code{\link{derivePHMM}}
#' @examples
#'  data(casino)
#'  deriveHMM(list(casino))
################################################################################
deriveHMM <- function(x, seqweights = NULL, residues = NULL, states = NULL,
                      modelend = FALSE, pseudocounts = "background",
                      logspace = TRUE){
  if(!(is.list(x))) stop("x must be a list of named vectors")
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  namesok <- all(sapply(x, function(y) !is.null(names(y)) | length(y) == 0))
  if(!(namesok)) stop("x must be a list of named vectors")
  residues <- .alphadetect(x, residues = residues)
  if(is.null(states)) states <- unique(unlist(lapply(x, names), use.names = FALSE))
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
  indices <- structure(seq_along(states) - 1, names = states)
  pathscoded <- lapply(x, function(e) indices[c("Begin", names(e), "Begin")])
  Acounts <- .acount(pathscoded[[1]], arity = nstates) * seqweights[1]
  if(n > 1){
    for(i in 2:n){
      Acounts <- Acounts + .acount(pathscoded[[i]], arity = nstates) * seqweights[i]
    }
  }
  indices <- structure(0:(nstates - 2), names = states[-1])
  statescoded <- lapply(x, function(e) indices[names(e)])
  if(DNA){
    rescoded <- lapply(x, .encodeDNA, arity = 4)
  }else if(AA){
    rescoded <- lapply(x, .encodeAA, arity = 20)
  }else{
    indices <- structure(0:(nres - 1), names = residues)
    rescoded <- lapply(x, function(e) indices[e])
  }
  Ecounts <- .ecount(statescoded[[1]], nstates - 1, rescoded[[1]], nres) * seqweights[1]
  if(n > 1) for(i in 2:n){
    Ecountsi <- .ecount(statescoded[[i]], nstates - 1, rescoded[[i]], nres) * seqweights[i]
    Ecounts <- Ecounts + Ecountsi
  }
  # add pseudocounts
  if(identical(pseudocounts, "background")){
    Apseudocounts <- Acounts + 1
    Apseudocounts <- Apseudocounts/sum(Apseudocounts) *
      (nstates^2 - if(modelend) 0 else nstates)
    Acounts <- Acounts + Apseudocounts
    Epseudocounts <- Ecounts + 1
    Epseudocounts <- Epseudocounts/sum(Epseudocounts) * (nstates - 1) * nres
    Ecounts <- Ecounts + Epseudocounts
  }else if(identical(pseudocounts, "Laplace")){
    Acounts <- Acounts + 1
    Ecounts <- Ecounts + 1
  }else if(is.list(pseudocounts)){
    stopifnot(length(pseudocounts == 2))
    Acounts <- Acounts + pseudocounts[[1]]
    Ecounts <- Ecounts + pseudocounts[[2]]
  }else if(!identical(pseudocounts, "none")){
    stop("invalid pseudocounts argument")
  }
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
