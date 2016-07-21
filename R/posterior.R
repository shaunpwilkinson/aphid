posterior.HMM <- function(x, y){
  n <- length(y)
  states <- rownames(x$A)[-1]
  H <- length(states)
  back <- backward(x, y)
  B <- back$array
  forw <- forward(x, y)
  R <- forw$array
  logPx <- forw$score
  postprobs <- exp(R + B - logPx)
  return(postprobs)
}
