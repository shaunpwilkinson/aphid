#' Posterior decoding.
#'
#' Calculate the posterior probability of a sequence given a model.
#'
#' @param x an object of class \code{'HMM'} or \code{'PHMM'}.
#' @param y a vector of mode "character" or "raw" (a "DNAbin" or "AAbin"
#'   object) representing a single sequence hypothetically emitted by
#'   the model in \code{x}.
#' @param logspace logical indicating whether the emission and transition
#'   probabilities of x are logged. If \code{logspace = "autodetect"}
#'   (default setting), the function will automatically detect
#'   if the probabilities are logged, returning an error if
#'   inconsistencies are found. Note that choosing the latter option
#'   increases the computational overhead; therefore specifying
#'   \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param cpp logical, indicates whether the dynamic programming matrix
#'   should be filled using compiled C++ functions (default; many times faster).
#'   The FALSE option is primarily retained for bug fixing and experimentation.
#' @return a vector or matrix of posterior probabilities.
#' @details TBA
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @examples
#'   ## Posterior decoding for standard hidden Markov models
#'   ## The dishonest casino example from Durbin et al. (1998) chapter 3.2
#'   A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9),
#'               nrow = 3) # transition probability matrix
#'   dimnames(A) <- list(from = c("Begin", "Fair", "Loaded"),
#'                       to = c("Begin", "Fair", "Loaded"))
#'   E <- matrix(c((1/6), (1/6), (1/6), (1/6), (1/6), (1/6),
#'                 (1/10), (1/10), (1/10), (1/10), (1/10), (1/2)),
#'               nrow = 2, byrow = TRUE) # emission probability matrix
#'   dimnames(E) <- list(states = c('Fair', 'Loaded'), residues = paste(1:6))
#'   x <- structure(list(A = A, E = E), class = "HMM") # create hidden Markov model
#'   plot(x, main = "Dishonest casino HMM")
#'   data(casino)
#'   casino.post <- posterior(x, casino)
#'   plot(1:300, casino.post[1, ], type = "l", xlab = "Roll number",
#'        ylab = "Posterior probability of dice being fair",
#'        main = "The dishonest casino")
#' @seealso \code{\link{forward}}, \code{\link{backward}}, \code{\link{Viterbi}}
#' @name posterior
################################################################################
posterior <- function(x, y, logspace = "autodetect", cpp = TRUE){
  UseMethod("posterior")
}
################################################################################
#' @rdname posterior
################################################################################
posterior.HMM <- function(x, y, logspace = "autodetect", cpp = TRUE){
  if(identical(logspace, 'autodetect')) logspace <- .logdetect(x)
  back <- backward(x, y, logspace = logspace, cpp = cpp)
  B <- back$array
  forw <- forward(x, y, logspace = logspace, cpp = cpp)
  R <- forw$array
  logPx <- forw$score
  postprobs <- exp(R + B - logPx)
  return(postprobs)
}
################################################################################
#' @rdname posterior
################################################################################
posterior.PHMM <- function(x, y, logspace = "autodetect", cpp = TRUE){
  if(identical(logspace, 'autodetect')) logspace <- .logdetect(x)
  back <- backward(x, y, logspace = logspace, cpp = cpp)
  B <- back$array
  forw <- forward(x, y, logspace = logspace, cpp = cpp)
  R <- forw$array
  logPx <- forw$score
  postprobs <- exp(R[,,"M"] + B[,,"M"] - logPx)
  maxprobs <- apply(postprobs, 2, max)
  return(maxprobs)
}
################################################################################
