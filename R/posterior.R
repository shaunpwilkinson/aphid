#' Posterior decoding.
#'
#' Calculates the posterior probability of a sequence given a model.
#'
#' @param x an object of class \code{'HMM'} or \code{'PHMM'}.
#' @param y a character vector...
#' @param logspace logical argument indicating whether the emission and transition
#' probabilities of x are logged (base e; TRUE) or raw (FALSE). Alternatively, if
#' \code{logspace = "autodetect"} (default), the function will automatically detect
#' if the probabilities are in log space, returning an error if inconsistencies are found.
#' Note that choosing the latter option increases the computational
#' overhead; therefore specifying \code{TRUE} or \code{FALSE} can reduce the running time.

posterior.HMM <- function(x, y, logspace = "autodetect"){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  n <- length(y)
  states <- rownames(x$E)
  H <- length(states)
  back <- backward(x, y, logspace = logspace)
  B <- back$array
  forw <- forward(x, y, logspace = logspace)
  R <- forw$array
  logPx <- forw$score
  postprobs <- exp(R + B - logPx)
  return(postprobs)
}
