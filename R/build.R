#' Manually build HMMs or profile HMMs
#'
#' \code{buildPHMM} and \code{buildHMM} are functions for the manual creation of
#' objects of class \code{PHMM} and \code{HMM}, respectively.
#'
#' @param A matrix of transition probabilities with \code{"dimnames"} attributes
#' following \code{dimnames(A) = List(from = states, to = states)}
#' @param E matrix of emission probabilities with \code{"dimnames"} attributes
#' following \code{dimnames(E) = List(state = states, residue = residues)}
#' @param qe a named vector of profile HMM background emission frequencies
#' @param qa a 3 x 3 matrix of profile HMM background transition frequencies
#'
buildPHMM <- function(E, A, qe, qa, n){
  states <- c("M", "I", "D")
  if(!(all(dimnames(A)[[1]] == states) & all(dimnames(A)[[3]] == states))){
    stop("names for dimensions 1 and 3 of transitions
         array (A) must be c('M', 'I', 'D')")
  }
  if(!(ncol(A) == ncol(E) + 1)){
    stop("transitions array (A) should include a begin state")
  }
  if(!all(names(qe) == rownames(E))){
    stop("residue names for background emissions vector (qe) should be identical
         to rownames of emissions matrix (E)")
  }
  if(!(all(colnames(qa) == states) & all(rownames(qa) == states))){
    stop("rownames and colnames of background transitions matrix
         must be c('M', 'I', 'D')")
  }
  res <- structure(list(A = A, E = E, qe = qe, qa = qa, n), class = 'PHMM')
  return(res)
}
buildHMM <- function(A, E){
  if(is.null(colnames(A)) | is.null(colnames(E))) stop(
    'both A and E must have dimnames attributes')
  states <- rownames(A)
  if(!(identical(states, colnames(A)) &
       identical(states[-1], rownames(E)))){
    stop("rownames and colnames of the transitions matrix (A), and
         rownames of the emissions matrix (E) must all be non-NULL
         and identical")
  }
  if(rownames(A)[1] != "BeginEnd") stop(
    "Transitions matrix (A) should include a 'BeginEnd' state at row 1
    and col 1")
  res <- structure(list(A = A, E = E), class = 'HMM')
  return(res)
}


