
backward.HMM <- function (x, obs, logspace = FALSE){
  n <- length(obs)
  states <- names(x$s)
  H <- length(states)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  s <- if(logspace) x$s else log(x$s)
  B <- array(NA, dim = c(H, n), dimnames = list(states = states, index = 1:n))
  B[, n] <- 0
  fun <- Vectorize(function(k, l) A[k, l] + E[l, obs[i]] + B[l, i]) 
  for (i in n:2) B[, i - 1] <- apply(outer(1:H, 1:H, fun), 1, logsum)
  logFullProbs <- rep(NA, H)
  names(logFullProbs) <- states
  for(i in states) logFullProbs[i] <- s[i] + E[i, obs[1]] + B[i, 1] #termination
  res <- structure(list(logFullProb = logsum(logFullProbs),
                        backwardArray = B), class = 'backward')
  return(res)
}


