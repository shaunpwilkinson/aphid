#' Build HMMs or profile HMMs
#' 
#' \code{buildPHMM} and \code{buildHMM} are functions for the manual creation of
#' objects of class \code{PHMM} and \code{HMM}, respectively. 
#' 
#' @param E a matrix of emission probabilities
#' @param A a matrix of transition probabilities
#' @param qe a named vector of profile HMM background emission frequencies
#' @param qa a 3 x 3 matrix of profile HMM background transition frequencies
#' @param s a named vector of HMM start probabilities
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
buildHMM <- function(s, A, E){
  if(is.null(names(s))) stop('all arguments must have names attributes')
  states <- names(s)
  if(!(identical(states, colnames(A)) & 
       identical(states, rownames(A)))){
    stop("rownames and colnames of the transitions matrix (A) must be 
         the same as the names of the startprobs argument (s)")
  }
  if(!(identical(states, rownames(E)))){
    stop("rownames of the emissions matrix (E) must be the 
         same as the names of the startprobs argument (s)")
  }
  res <- structure(list(s = s, A = A, 
                        E = E), class = 'HMM')
  return(res)
}


