#' Iterative estimation of model parameters.
#'
#' \code{BaumWelch} uses a special case of the expectation-maximization (EM) algorithm
#' to find the locally-optimal parameters of a HMM or PHMM.
#' @param x an object of class \code{'HMM'} or \code{'PHMM'} specifying the
#' starting parameter values.
#' @param y a list of training sequences whose hidden states are unknown.
#' @param maxiter the maximum number of EM iterations before the process is terminated
#' @param deltaLL the change in log likelihood specified as the convergence threshold
#' @param logspace logical argument indicating whether the emission
#' and transmission probabilities povided for the model(s) are logged.
#' @param quiet logical argument indicating whether the iteration progress should
#' be suppressed (TRUE) or printed to the console (FALSE; default).
#' @param pseudocounts used to account for the possible absence of certain transition
#' and/or emission types in the training dataset.
#' either \code{'Laplace'} (adds one of each possible transition and emission type to the
#' training dataset; default), \code{'none'}, or a two-element list containing a matrix of
#' transition pseudocounts as its first element and a matrix of emission pseudocounts
#' as its second. If a list is supplied both matrices must have row and column names
#' according to the residues (column names of emission matrix) and states
#' (row and column names of the transition matrix and row names of the emission matrix).
#' The first row and column of the transition matrix must be 'BeginEnd'.
#' @seealso \code{\link{deriveHMM}} and \code{\link{derivePHMM}} for
#' maximum-likelihood parameter estimation when training sequence states are
#' known.
#' @name BaumWelch
#'
BaumWelch <- function(x, y, maxiter = 100, deltaLL = 1E-07,
                          logspace = FALSE, quiet = FALSE,
                          pseudocounts = "Laplace"){
  UseMethod("BaumWelch")
}

#' @rdname BaumWelch
BaumWelch.HMM <- function(x, y, maxiter = 100, deltaLL = 1E-07,
                      logspace = FALSE, quiet = FALSE,
                      modelend = FALSE, pseudocounts = "Laplace"){
  if(is.list(y)){
  } else if(is.vector(y, mode = "character")){
    y <- list(y)
  } else stop("invalid y argument")
  n <- length(y)
  seqlengths <- lapply(y, length)
  states <- rownames(x$A)
  nstates <- length(states)
  residues <- colnames(x$E)
  nres <- length(residues)
  Apseudocounts <- matrix(nrow = nstates, ncol = nstates)
  Epseudocounts <- matrix(nrow = nstates - 1, ncol = nres)
  dimnames(Apseudocounts) <- list(from = states, to =  states)
  dimnames(Epseudocounts) <- list(state = states[-1], residue = residues)
  if(identical(pseudocounts, "Laplace")){
    Apseudocounts[] <- Epseudocounts[] <- 1
    if(!modelend) Apseudocounts[, 1] <- 0
  } else if(identical(pseudocounts, "none")){
    AApseudocounts[] <- Epseudocounts[] <- 0
  } else if(is.list(pseudocounts)){
    stopifnot(length(pseudocounts == 2))
    Apseudocounts[] <- pseudocounts[[1]]
    Epseudocounts[] <- pseudocounts[[2]]
  } else stop("invalid 'pseudocounts' argument")
  out <- x
  if(!logspace){
    out$E <- log(out$E)
    out$A <- log(out$A)
  }
  E <- out$E
  A <- out$A
  LL <- -1E06
  for(i in 1:maxiter){
    tmpA <- Apseudocounts
    tmpE <- Epseudocounts
    tmplogPx <- rep(NA, n)
    for(j in 1:n){
      yj <- y[[j]]
      nj <- seqlengths[[j]]
      if(nj == 0){
        tmpA[1, 1] <- tmpA[1, 1] + if(modelend) 1 else 0
      }else{
        forwj <- forward(out, yj, logspace = TRUE)
        Rj <- forwj$array
        logPxj <- forwj$score
        tmplogPx[j] <- logPxj
        backj <- backward(out, yj, logspace = TRUE)
        Bj <- backj$array
        for(k in states[-1]){
          tmpA[1, -1] <- tmpA[1, -1] + exp(A[1, -1] + E[, yj[1]] + Bj[, 1] - logPxj)
          tmpA[-1, 1] <- tmpA[-1, 1] + exp(Rj[, nj] + A[-1, 1] - logPxj)
          for(l in states[-1]){
            tmpA[k, l] <- tmpA[k, l] + exp(logsum(Rj[k, -nj] + A[k, l] +
                                                    E[l, yj[-1]] + Bj[l, -1])
                                           - logPxj)
          }
          for(b in residues){
            cond <- yj == b
            tmpE[k, b] <- tmpE[k, b] + exp(logsum(Rj[k, cond] + Bj[k, cond]) - logPxj)
          }
        }
      }
    }
    A[] <- log(tmpA/apply(tmpA, 1, sum))
    E[] <- log(tmpE/apply(tmpE, 1, sum))
    out$A <- A
    out$E <- E
    logPx <- sum(tmplogPx) # page 62 eq 3.17
    if(!quiet) cat("Iteration", i, "log likelihood = ", logPx, "\n")
    if(abs(LL - logPx) < deltaLL){
      if(!logspace){
        out$A <- exp(out$A)
        out$E <- exp(out$E)
      }
      if(!quiet) cat("Convergence threshold reached after", i, "EM iterations\n")
      return(out)
    }
    LL <- logPx
  }
  stop("Failed to converge on a local maximum. Try increasing 'maxiter',
        decreasing 'deltaLL' or modifying start parameters")
}

#' @rdname BaumWelch
BaumWelch.PHMM <- function(x, y, maxiter = 100, deltaLL = 1E-07,
                          logspace = FALSE, quiet = FALSE,
                          pseudocounts = "Laplace"){
  stop("BaumWelch method for PHMMs is work in progress")
}

