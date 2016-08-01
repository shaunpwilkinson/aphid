#' Viterbi training.
#'
#' Update model parameters using a list of training sequences. Also known as the
#' segmental K-means algorithm (Juang & Rabiner 1990)
#'
#' @param x the initial model with starting parameters to be updated.
#' @param y a list of character vectors...
#' @references Juang B-H & Rabiner L R (1990) The segmental K-means
#' algorithm for estimating parameters of hidden Markov models.
#' IEEE transactions on Acoustics,
#' Speech...
#'
train <- function(x, y, maxiter = 100, logspace = FALSE, quiet = FALSE,
                  modelend = FALSE, pseudocounts = "Laplace"){
  UseMethod("train")
}

#' @rdname train
train.HMM <- function(x, y, maxiter = 100, logspace = FALSE, quiet = FALSE,
                  modelend = FALSE, pseudocounts = "Laplace"){
  if(is.list(y)){
  } else if(is.vector(y, mode = "character")){
    y <- list(y)
  } else stop("invalid y argument")
  n <- length(y)
  states <- rownames(x$A)
  residues <- colnames(x$E)
  out <- x
  if(!logspace){
    out$E <- log(out$E)
    out$A <- log(out$A)
  }
  for(i in 1:maxiter){
    samename <- rep(FALSE, n)
    for(j in 1:n){
      vitj <- Viterbi(out, y[[j]], logspace = TRUE)
      if(identical(vitj$path, names(y[[j]]))) samename[j] <- TRUE
      names(y[[j]]) <- vitj$path
    }
    if(all(samename)){
      if(!logspace){
        out$A <- exp(out$A)
        out$E <- exp(out$E)
      }
      if(!quiet) cat("Iteration", i, "\nConverged after", i, "iterations\n")
      return(out)
    }else{
      if(!quiet) cat("Iteration", i, "\n")
      out <- deriveHMM(y, residues = residues, states = states, modelend = modelend,
                          pseudocounts = pseudocounts, logspace = TRUE)
    }
  }
  stop("Failed to converge. Try increasing 'maxiter' or modifying start parameters")
}
