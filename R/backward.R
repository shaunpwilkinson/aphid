
backward.HMM <- function (x, y, logspace = FALSE){
  n <- length(y)
  states <- rownames(x$A)[-1]
  H <- length(states)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  #s <- if(logspace) x$s else log(x$s)
  s <-
  B <- array(NA, dim = c(H, n), dimnames = list(state = states, rolls = 1:n))
  B[, n] <- 0
  fun <- Vectorize(function(k, l) A[k + 1, l + 1] + E[l, y[i]] + B[l, i])
  for (i in n:2) B[, i - 1] <- apply(outer(1:H, 1:H, fun), 1, logsum)
  logprobs <- rep(NA, H)
  names(logprobs) <- states
  for(i in states) logprobs[i] <- A[1, i] + E[i, y[1]] + B[i, 1] #termination
  res <- structure(list(score = logsum(logprobs), array = B), class = 'backward')
  return(res)
}


