posterior.HMM <- function(x, obs){
  n <- length(obs)
  states <- names(x$s)
  H <- length(states)
  back <- backward(x, obs)
  B <- back$backwardArray
  forw <- forward(x, obs)
  R <- forw$forwardArray
  logPx <- forw$logFullProb
  postProbs <- exp(R + B - logPx)
  return(postProbs)
}
