#' Print summary methods.
#' @param x object of various classes.
#' @param ... additional arguments to be passed between methods.
#' @return NULL (invisibly)
#' @author Shaun Wilkinson
#' @name print
print.PHMM <- function(x, ...){
  cat("Profile hidden Markov model (object class: 'PHMM')\n",
      "with ",
      x$size,
      " internal positions emitting ",
      nrow(x$E),
      " unique residues\n",
      "(",
      paste(rownames(x$E), collapse = ", "),
      ").\n",
      sep = "")
}
#' @rdname print
print.HMM <- function(x, ...){
  cat("Hidden Markov model (object class: 'HMM') with ",
      nrow(x$E),
      " hidden states (",
      paste(rownames(x$E), collapse = ", "),
      ") emitting ",
      ncol(x$E),
      " unique residues (",
      paste(colnames(x$E), collapse = ", "),
      ").\n",
      sep = "")
}
#' @rdname print
print.fullprob <- function(x, ...){
  if(x$odds){
    cat("Log odds score: ", x$score)
  } else cat("Full (log) probability of sequence given model =", x$score)
}
#' @rdname print
print.Viterbi <- function(x, ...){
  cat("Optimal path with length",
      length(x$path),
      "and score",
      x$score)
}


# digits integer, for compatibility with other print methods.
