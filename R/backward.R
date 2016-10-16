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
#' @param odds logical indicating whether the scores in th edynamic programming matrix
#' should be odds ratios (TRUE; default) or full (log) probabilities.
#' (\code{FALSE}; default) or if the log-odds score should be returned
#' (\code{TRUE}).
#' @param type character string indicating whether insert and delete states
#' at the beginning and end of the path should count towards the final score
#' ('global'; default). Note that semiglobal and local models
#' are not currently supported in this version.
#' @param DI logical. should delete-insert transitions be allowed? Only applicable for
#' objects of class \code{"PHMM"}.
#' @param ID logical. should insert-delete transitions be allowed? Only applicable for
#' objects of class \code{"PHMM"}.
#' @name backward
#'
backward <- function(x, y, qe = NULL, logspace = "autodetect",  odds = TRUE,
                     type = "global", DI = TRUE, ID = TRUE){
  UseMethod("backward")
}

#' @rdname backward
backward.PHMM <- function (x, y, qe = NULL, logspace = "autodetect",  odds = TRUE,
                           type = "global", DI = TRUE, ID = TRUE){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
  if(type %in% c('semiglobal','local'))
    stop("semiglobal and local models are not supported in this version")
  pp <- inherits(y, 'PHMM')
  if(pp) stop("PHMM/PHMM backward algorithm not implemented in this version")
  pd <- inherits(y, "DNAbin")
  if(pd){
    x$E <- x$E[order(rownames(x$E)),]
    if(!(identical(rownames(x$E), c("a", "c", "g", "t")) |
         identical(rownames(x$E), c("A", "C", "G", "T")))){
      stop("invalid model, residue alphabet does not correspond to
           nucleotide alphabet")
    }
    if(is.list(y)) y <- y[[1]]
    if(is.matrix(y)){
      if(nrow(y) == 1) {
        y <- as.vector(y)
      } else {
        stop("backward algorithm can only process one sequence at a time")
      }
    }
  }
  n <- ncol(x$E)
  m <- if(pp) ncol(y$E) else length(y)
  states <- if(pp) c("MI", "DG", "MM", "GD", "IM") else c("D", "M", "I")
  B <- array(-Inf, dim = c(n + 1, m + 1, length(states)))
  dimnames(B) <- list(x = 0:n, y = 0:m, state = states)
  # background emission probabilities
  if(!(is.null(qe))){
    if(all(qe >= 0) & all(qe <= 1)) qe <- log(qe)
  }else if(pp){
    if(!is.null(x$qe) | !is.null(y$qe)){
      if(!is.null(x$qe) & !is.null(y$qe)){
        qe <- if(logspace) (exp(x$qe) + exp(y$qe))/2 else (x$qe + y$qe)/2
        qe <- log(qe)
      }else if(!is.null(x$qe)){
        qe <- if(logspace) x$qe else log(x$qe)
      }else{
        qe <- if(logspace) y$qe else log(y$qe)
      }
    }else{
      qe <- log(rep(1/nrow(x$E), nrow(x$E)))
      names(qe) <- rownames(x$E)
    }
  }else{
    if(!is.null(x$qe)){
      qe <- if(logspace) x$qe else log(x$qe)
    }else{
      qe <- log(rep(1/nrow(x$E), nrow(x$E)))
      names(qe) <- rownames(x$E)
    }
  }
  if(pp){
    ### placeholder
  }else{
    qey <- if(odds) rep(0, m) else if(pd) sapply(y, DNAprobC, qe) else qe[y]
    A <- if(logspace) x$A else log(x$A)
    E <- if(logspace) x$E else log(x$E)
    if(odds) E <- E - qe
    B[n + 1, m + 1, ] <- A[c("DM", "MM", "IM"), n + 1]
    if(type == "global"){ #key: D = 1, M = 2, I = 3
     for(i in n:1) {
       B[i, m + 1, "D"] <- B[i + 1, m + 1, "D"] + A["DD", i]
       B[i, m + 1, "M"] <- B[i + 1, m + 1, "D"] + A["MD", i]
       B[i, m + 1, "I"] <- if(ID) B[i + 1, m + 1, "D"] + A["ID", i] else -Inf
      }
     for(j in m:1) {
       B[n + 1, j, "D"] <- if(DI) B[n + 1, j + 1, "I"] + A["DI", n + 1] + qey[j] else -Inf
       B[n + 1, j, "M"] <- B[n + 1, j + 1, "I"] + A["MI", n + 1] + qey[j]
       B[n + 1, j, "I"] <- B[n + 1, j + 1, "I"] + A["II", n + 1] + qey[j]
      }
    }else{
      # B[-1, 1, 1] <- B[1, -1, 3] <- 0 ### check this
    }
    for(i in n:1){
      for(j in m:1){
        sij <- if(pd) DNAprobC(y[j], E[, i]) else E[y[j], i]
        Dcdt <- c(B[i + 1, j, "D"] + A["DD", i],
                  B[i + 1, j + 1, "M"] + A["DM", i] + sij,
                  if(DI) B[i, j + 1, "I"] + A["DI", i] + qey[j] else -Inf)
        Mcdt <- c(B[i + 1, j, "D"] + A["MD", i],
                  B[i + 1, j + 1, "M"] + A["MM", i] + sij,
                  B[i, j + 1, "I"] + A["MI", i] + qey[j])
        Icdt <- c(if(ID) B[i + 1, j, "D"] + A["ID", i] else -Inf,
                  B[i + 1, j + 1, "M"] + A["IM", i] + sij,
                  B[i, j + 1, "I"] + A["II", i] + qey[j])
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
                                      class = 'fullprob'))
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


