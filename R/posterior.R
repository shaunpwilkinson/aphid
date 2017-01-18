#' Posterior decoding.
#'
#' Calculate the posterior probability of a sequence given a model.
#'
#' @param x an object of class \code{'HMM'} or \code{'PHMM'}.
#' @param y a character vector.
#' @param logspace logical argument indicating whether the emission and transition
#' probabilities of x are logged (base e; TRUE) or raw (FALSE). Alternatively, if
#' \code{logspace = "autodetect"} (default), the function will automatically detect
#' if the probabilities are in log space, returning an error if inconsistencies are found.
#' Note that choosing the latter option increases the computational
#' overhead; therefore specifying \code{TRUE} or \code{FALSE} can reduce the running time.
#' @return a vector or matrix of posterior probabilities.
#' @name posterior
#' @export
#'
posterior <- function(x, y, logspace = "autodetect", cpp = TRUE){
  UseMethod("posterior")
}


#' @rdname posterior
#' @export
#'
posterior.HMM <- function(x, y, logspace = "autodetect", cpp = TRUE){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  back <- backward(x, y, logspace = logspace, cpp = cpp)
  B <- back$array
  forw <- forward(x, y, logspace = logspace, cpp = cpp)
  R <- forw$array
  logPx <- forw$score
  postprobs <- exp(R + B - logPx)
  return(postprobs)
}


#' @rdname posterior
#' @export
#'
posterior.PHMM <- function(x, y, logspace = "autodetect", cpp = TRUE){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  back <- backward(x, y, logspace = logspace, cpp = cpp)
  B <- back$array
  forw <- forward(x, y, logspace = logspace, cpp = cpp)
  R <- forw$array
  logPx <- forw$score
  postprobs <- exp(R[,,"M"] + B[,,"M"] - logPx)
  maxprobs <- apply(postprobs, 2, max)
  return(maxprobs)
}

