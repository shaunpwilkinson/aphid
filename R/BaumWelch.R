BaumWelch <- function(x, y, maxiter = 100, deltaLL = 1E-07,
                      logspace = FALSE, quiet = FALSE){
  if(!(inherits(x, 'HMM'))) stop("x must be an object of class 'HMM'")
  n <- length(y)
  states <- rownames(x$A)[-1]
  H <- length(states)
  residues <- colnames(x$E)
  reslen <- length(residues)
  out <- x
  if(!logspace){
    out$E <- log(out$E)
    out$A <- log(out$A)
    #out$s <- log(out$s)
  }
  E <- out$E
  A <- out$A
  #s <- out$s
  LL <- -1E06
  Akl <- Vectorize(function(k, l)
    exp(logsum(R[k, -n] + A[k + 1, l + 1] + E[l, y[-1]] + B[l, -1]) - logPx))
  Ekb <- Vectorize(function(k, b){
    cond <- y == residues[b]
    return(exp(logsum(R[k, cond] + B[k, cond]) - logPx))
  })
  for(i in 1:maxiter){
    forw <- forward(out, y, logspace = TRUE)
    R <- forw$array
    logPx <- forw$score
    if(!quiet) cat("Iteration", i, "log likelihood = ", logPx, "\n")
    if(LL - logPx < -deltaLL){
      back <- backward(out, y, logspace = TRUE)
      B <- back$array
      tmp <- outer(1:H, 1:H, Akl)
      A[-1, -1] <- log(tmp/apply(tmp, 1, sum))
      tmp <- outer(1:H, 1:reslen, Ekb)
      E[] <- log(tmp/apply(tmp, 1, sum))
      out$A <- A
      out$E <- E
      LL <- logPx
    }else{
      if(!logspace){
        out$A <- exp(out$A)
        out$E <- exp(out$E)
        #out$s <- exp(out$s)
      }
      if(!quiet) cat("Model converged after", i, "EM iterations\n")
      return(out)
    }
  }
  stop("Model did not converge")
}

