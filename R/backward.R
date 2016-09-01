#' Backward algorithm.
#'
#' Backward algorithm for calculating the full (log) probability
#' of a sequence given a hidden Markov model.
#'
#' @param x an object of class \code{PHMM} or \code{HMM}.
#' @param y a character vector representing a single instance of a sequence
#' hypothetically emitted by the model.
#' @param logspace logical argument indicating whether the emission and transition
#' probabilities of x are logged (base e; TRUE) or raw (FALSE). Alternatively, if
#' \code{logspace = "autodetect"} (default), the function will automatically detect
#' if the probabilities are in log space, returning an error if inconsistencies are found.
#' Note that choosing the latter option increases the computational
#' overhead; therefore specifying \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param odds logical indicating whether the full (log) probability is required
#' (\code{FALSE}; default) or if the log-odds score should be returned
#' (\code{TRUE}).
#' @param type character string indicating whether insert and delete states
#' at the beginning and end of the path should count towards the final score
#' ('global'; default). Note that semiglobal and local models
#' are not currently supported in this version.
#' @name backward
#'
backward <- function(x, y, logspace = "autodetect",  odds = FALSE, type = "global"){
  UseMethod("backward")
}

#' @rdname backward
backward.PHMM <- function (x, y, logspace = "autodetect",  odds = FALSE, type = "global"){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
  if(type %in% c('semiglobal','local'))
    stop("semiglobal and local models are not supported in this version")
  pp <- inherits(y, 'PHMM')
  n <- ncol(x$E)
  m <- if(pp) ncol(y$E) else length(y)
  states <- if(pp) c("MI", "DG", "MM", "GD", "IM") else c("D", "M", "I")
  B <- array(-Inf, dim = c(n + 1, m + 1, length(states)))
  dimnames(B) <- list(x = 0:n, y = 0:m, state = states)
  if(pp){
    stop("backward algorithm for PHMM vs PHMM not supported yet")
    #
    #
  }else{
    if(!is.null(x$qe)){
      qe <- if(logspace) x$qe else log(x$qe)
    }else{
      qe <- log(rep(1/nrow(x$E), nrow(x$E)))
      names(qe) <- rownames(x$E)
    }
    A <- if(logspace) x$A else log(x$A)
    E <- if(logspace) x$E else log(x$E)
    if(odds) E <- E - qe
    B[n + 1, m + 1, ] <- A[ , n + 1, "M"]
    if(type == "global"){ #key: D = 1, M = 2, I = 3
     for(i in n:1) {
       B[i, m + 1, "D"] <- B[i + 1, m + 1, "D"] + A["D", i, "D"]
       B[i, m + 1, "M"] <- B[i + 1, m + 1, "D"] + A["M", i, "D"]
       B[i, m + 1, "I"] <- B[i + 1, m + 1, "D"] + A["I", i, "D"]
      }
     for(j in m:1) {
       B[n + 1, j, "D"] <- B[n + 1, j + 1, "I"] + A["D", n + 1, "I"] +
         if(odds) 0 else qe[y[j]]
       B[n + 1, j, "M"] <- B[n + 1, j + 1, "I"] + A["M", n + 1, "I"] +
         if(odds) 0 else qe[y[j]]
       B[n + 1, j, "I"] <- B[n + 1, j + 1, "I"] + A["I", n + 1, "I"] +
         if(odds) 0 else qe[y[j]]
      }
    }else{
      # B[-1, 1, 1] <- B[1, -1, 3] <- 0 ### check this
    }
    for(i in n:1){
      for(j in m:1){
        sij <- E[y[j], i]
        Dcdt <- c(B[i + 1, j, "D"] + A["D", i, "D"],
                  B[i + 1, j + 1, "M"] + A["D", i, "M"] + sij,
                  B[i, j + 1, "I"] + A["D", i, "I"] + if(odds) 0 else qe[y[j]])
        Mcdt <- c(B[i + 1, j, "D"] + A["M", i, "D"],
                  B[i + 1, j + 1, "M"] + A["M", i, "M"] + sij,
                  B[i, j + 1, "I"] + A["M", i, "I"] + if(odds) 0 else qe[y[j]])
        Icdt <- c(B[i + 1, j, "D"] + A["I", i, "D"],
                  B[i + 1, j + 1, "M"] + A["I", i, "M"] + sij,
                  B[i, j + 1, "I"] + A["I", i, "I"] + if(odds) 0 else qe[y[j]])
        B[i, j, "D"] <- logsum(Dcdt)
        B[i, j, "M"] <- logsum(Mcdt)
        B[i, j, "I"] <- logsum(Icdt)
      }
    }
    score <- B[1, 1, "M"]
    res <- structure(list(score = score, array = B, odds = odds),
                     class = 'fullprob')
    return(res)
  }
}

#' @rdname backward
backward.HMM <- function (x, y, logspace = "autodetect"){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
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
  res <- structure(list(score = score, array = B, odds = FALSE),
                   class = 'fullprob')
  return(res)
}


