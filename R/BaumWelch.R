BaumWelch <- function(x, obs, iterations = 100, 
                      delta.llk = 0.0000001, logspace = FALSE,
                      quiet = FALSE){
  if(!(inherits(x, 'HMM'))) stop("x must be an object of class 'HMM'")
  n <- length(obs)
  states <- names(x$s)
  H <- length(states)
  residues <- colnames(x$E)
  symblen <- length(residues)
  res <- x
  if(!logspace){
    res$E <- log(res$E)
    res$A <- log(res$A)
    res$s <- log(res$s)
  }
  E <- res$E
  A <- res$A
  s <- res$s
  llk <- -1000000
  Akl <- Vectorize(function(k, l) 
    exp(logsum(R[k, -n] + A[k, l] + E[l, obs[-1]] + B[l, -1]) - logPx))
  Ekb <- Vectorize(function(k, b){
    cond <- obs == residues[b]
    return(exp(logsum(R[k, cond] + B[k, cond]) - logPx))
  })
  for(i in 1:iterations){
    forw <- forward(res, obs, logspace = TRUE)
    R <- forw$forwardArray
    logPx <- forw$logFullProb
    if(!quiet) cat("Iteration", i, "log likelihood = ", logPx, "\n")
    if(llk - logPx < -delta.llk){
      back <- backward(res, obs, logspace = TRUE)
      B <- back$backwardArray
      tmp <- outer(1:H, 1:H, Akl)
      A[] <- log(tmp/apply(tmp, 1, sum))
      tmp <- outer(1:H, 1:symblen, Ekb)
      E[] <- log(tmp/apply(tmp, 1, sum))
      res$A <- A
      res$E <- E
      llk <- logPx
    }else{
      if(!logspace){
        res$A <- exp(res$A)
        res$E <- exp(res$E)
        res$s <- exp(res$s)
      }
      if(!quiet) cat("Model converged after", i, "iterations")
      return(res)
    }
  }
  stop("Model did not converge")
}

