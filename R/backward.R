
backward.HMM <- function (x, y, logspace = FALSE){
  n <- length(y)
  states <- rownames(x$E)
  H <- length(states)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  #s <- if(logspace) x$s else log(x$s)
  if(length(y) == 0) return(structure(list(score = A[1, 1], array = NULL),
                                      class = 'backward'))
  B <- array(NA, dim = c(H, n), dimnames = list(state = states, rolls = 1:n))
  B[, n] <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H) #ak0
  fun <- Vectorize(function(k, l) A[k + 1, l + 1] + E[l, y[i]] + B[l, i])
  for (i in n:2) B[, i - 1] <- apply(outer(1:H, 1:H, fun), 1, logsum)
  logprobs <- rep(NA, H)
  names(logprobs) <- states
  for(l in states) logprobs[l] <- A[1, l] + E[l, y[1]] + B[l, 1] #termination
  score <- logsum(logprobs)
  res <- structure(list(score = score, array = B), class = 'backward')
  return(res)
}


