#' Viterbi training.
#'
#' Update model parameters using a list of training sequences. Also known as the
#' segmental K-means algorithm (Juang & Rabiner 1990)
#'
#' @param x the initial model with starting parameters to be updated.
#' @param y a list of training sequences (character vectors) whose hidden states are unknown.
#' @references Juang B-H & Rabiner L R (1990) The segmental K-means
#' algorithm for estimating parameters of hidden Markov models.
#' IEEE transactions on Acoustics,
#' Speech...
#' @references Durbin...
#'
train <- function(x, y, maxiter = 100, logspace = "autodetect", quiet = FALSE,
                  modelend = FALSE, pseudocounts = "Laplace", fixqa = FALSE, fixqe = FALSE,
                  inserts = "threshold",
                  threshold = 0.5, lambda = 0, DI = TRUE){
  UseMethod("train")
}

#' @rdname train
train.PHMM <- function(x, y, maxiter = 100,
                       logspace = "autodetect", quiet = FALSE,
                       pseudocounts = "background",
                       fixqa = FALSE, fixqe = FALSE,
                       inserts = "threshold",
                       threshold = 0.5, lambda = 0, DI = TRUE){
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  if(is.list(y)){
  } else if(is.vector(y, mode = "character")){
    y <- list(y)
  } else stop("invalid y argument")
  n <- length(y)
  states <- c("D", "M", "I")
  residues <- rownames(x$E)
  out <- x
### what if input model has no qe/qa?
  if(!logspace){
    out$E <- log(out$E)
    out$A <- log(out$A)
    out$qa <- log(out$qa)
    out$qe <- log(out$qe)
  }
  alignment <- align2phmm(y, out, logspace = TRUE)
  for(i in 1:maxiter){
    if(!quiet) cat("Iteration", i, "\n")
    out <- derivePHMM(alignment, residues = residues, inserts = inserts,
               pseudocounts = pseudocounts, logspace = TRUE,
               qa = if(fixqa) out$qa else NULL,
               qe = if(fixqe) out$qe else NULL) ### add others too
    newalig <- align2phmm(y, out, logspace = TRUE)
    if(!identical(alignment, newalig)){
      alignment <- newalig
    } else{
      if(!logspace){
        out$A <- exp(out$A)
        out$E <- exp(out$E)
        out$qa <- exp(out$qa)
        out$qe <- exp(out$qe)
      }
      if(!quiet) cat("Converged after", i, "iterations\n")
      return(out)
    }
  }
  stop("Failed to converge. Try increasing 'maxiter' or modifying start parameters")
}


#' @rdname train
train.HMM <- function(x, y, maxiter = 100, logspace = "autodetect", quiet = FALSE,
                  modelend = FALSE, pseudocounts = "Laplace"){
  if(is.list(y)){
  } else if(is.vector(y, mode = "character")){
    y <- list(y)
  } else stop("invalid y argument")
  n <- length(y)
  states <- rownames(x$A)
  residues <- colnames(x$E)
  out <- x
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  if(!logspace){
    out$E <- log(out$E)
    out$A <- log(out$A)
  }
  for(i in 1:maxiter){
    samename <- logical(n)
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
