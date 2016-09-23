#' Forward algorithm.
#'
#' Forward algorithm for calculating the full log-odds
#' of a sequence given a HMM or PHMM and its random equivalent.
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
#' @param type character string indicating whether insert and delete states
#' at the beginning and end of the path should count towards the final score
#' ('global'; default). Note that semiglobal and local models
#' are not currently supported in this version.
#' @name forward

forward <- function(x, y, logspace = "autodetect", odds = FALSE,
                    type = "global"){
  UseMethod("forward")
}

#' @rdname forward
forward.PHMM <- function (x, y, logspace = "autodetect", odds = FALSE,
                          type = "global"){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
  if(type %in% c('semiglobal','local'))
    stop("semiglobal and local models are not supported in this version")
  pp <- inherits(y, 'PHMM')
  pd <- mode(y) == "raw"
  if(pd){
    x$E <- x$E[order(rownames(x$E)),]
    if(!(identical(rownames(x$E), c("a", "c", "g", "t")) |
         identical(rownames(x$E), c("A", "C", "G", "T")))){
      stop("invalid model, residue alphabet does not correspond to
           nucleotide alphabet")
    }
    if(is.matrix(y)){
      if(nrow(y) == 1) {
        y <- as.vector(y)
      } else {
        stop("forward algorithm can only process one sequence at a time")
      }
    }
  }
  n <- ncol(x$E)
  m <- if(pp) ncol(y$E) else length(y)
  states <- if(pp) c("MI", "DG", "MM", "GD", "IM") else c("D", "M", "I")
  R <- array(-Inf, dim = c(n + 1, m + 1, length(states)))
  dimnames(R) <- list(x = 0:n, y = 0:m, state = states)
  if(pp){
    stop("PHMM/PHMM forward algorithm not implemented in this version")
    #
    #
  }else{
    if(!is.null(x$qe)){
      qe <- if(logspace) x$qe else log(x$qe)
    }else{
      qe <- log(rep(1/nrow(x$E), nrow(x$E)))
      names(qe) <- rownames(x$E)
    }
    qey <- if(odds) rep(0, m) else if(pd) sapply(y, DNAprobC, qe) else qe[y]
    A <- if(logspace) x$A else log(x$A)
    E <- if(logspace) x$E else log(x$E)
    if(odds) E <- E - qe
    R[1, 1, 2] <- 0
    if(type == "global"){ #key: D = 1, M = 2, I = 3
      # R[-1, 1, 1] <- cumsum(c(0, A["D", 2:n, "D"])) + A["M", 1, "D"]
      # R[1, -1, 3] <- seq(from = A["M", 1, "I"], by = A["I", 1, "I"], length.out = m)
      R[2, 1, "D"] <- A["MD", 1]
      for(i in 2:n) R[i + 1, 1, "D"] <- R[i, 1, "D"] + A["DD", i]
      R[1, 2, "I"] <- A["MI", 1] + qey[1]
      for(j in 2:m) {
        R[1, j + 1, "I"] <- R[1, j, "I"] + A["II", 1] + qey[j]
      }
    }else{
      R[-1, 1, 1] <- R[1, -1, 3] <- 0 ### needs checking
    }
    for(i in 1:n){
      for(j in 1:m){
        sij <- if(pd) DNAprobC(y[j], E[, i]) else E[y[j], i]
        Dcdt <- c(R[i, j + 1, "D"] + A["DD", i],
                  R[i, j + 1, "M"] + A["MD", i],
                  R[i, j + 1, "I"] + A["ID", i])
        Mcdt <- c(R[i, j, "D"] + A["DM", i],
                  R[i, j, "M"] + A["MM", i],
                  R[i, j, "I"] + A["IM", i])
        Icdt <- c(R[i + 1, j, "D"] + A["DI", i + 1],
                  R[i + 1, j, "M"] + A["MI", i + 1],
                  R[i + 1, j, "I"] + A["II", i + 1])
        R[i + 1, j + 1, "D"] <- logsum(Dcdt)
        R[i + 1, j + 1, "M"] <- logsum(Mcdt) + sij
        R[i + 1, j + 1, "I"] <- logsum(Icdt) + qey[j]
      }
    }
    LLcdt <- c(R[n + 1, m + 1, "M"] + A["MM", n + 1],
               R[n + 1, m + 1, "I"] + A["IM", n + 1],
               R[n + 1, m + 1, "D"] + A["DM", n + 1])
    score <- logsum(LLcdt)
    res <- structure(list(score = score, array = R, odds = odds),
                     class = 'fullprob')
    return(res)
  }
}

#' @rdname forward
forward.HMM <- function (x, y, logspace = "autodetect"){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  n <- length(y)
  states <- rownames(x$E)
  H <- length(states)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  if(length(y) == 0) structure(list(score = A[1, 1], array = NULL), class = 'fullprob')
  R <- array(NA, dim = c(H, n), dimnames = list(state = states, rolls = 1:n))
  R[,1] <- E[, y[1]] + A[1, -1]
  fun <- Vectorize(function(k, l) R[k, i - 1] + A[k + 1, l + 1])
  for (i in 2:n) R[, i] <- E[, y[i]] + apply(outer(1:H, 1:H, fun), 2, logsum)
  ak0 <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H)
  score <- logsum(R[, n] + ak0)
  res <- structure(list(score = score, array = R, odds = FALSE),
                   class = 'fullprob')
  return(res)
}

#' @rdname forward
forward.default <- function(x, y, type = 'semiglobal', d = 8, e = 2,
                           S = NULL, itertab = NULL, offset = 0){
  ###check that x is a character vector###
  # this is incorrect, d and e have no probabilistic interpretation
  if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
  n <- length(x) + 1
  m <- length(y) + 1
  if (is.null(S)) {
    residues <- unique(c(x, y))
    S <- diag(2, nrow = length(residues))
    S <- S - 1
    dimnames(S) <- list(residues, residues)
  }
  # initialize scoring and pointer arrays (M and P)
  M <- array(-Inf, dim = c(n, m, 3))
  M[1, 1, 2] <- 0
  if(type == 'global'){
    M[, 1, 1] <- c(0, seq(from = -d, to = (- d + (n - 2) * -e), by = -e))
    M[1, , 3] <- c(0, seq(from = -d, to = (- d + (m - 2) * -e), by = -e))
  }else{
    M[2:n, 1, 1] <- M[1, 2:m, 3] <- 0 ### check this - should fill dim2 instead?
  }
  # recursion step
  for(i in if(is.null(itertab)) 2:n else 1:nrow(itertab)){
    for(j in if(is.null(itertab)) 2:m else 1){
      if(!is.null(itertab)){
        j <- itertab[i, 2]
        i <- itertab[i, 1]
      }
      sij <-  S[x[i - 1], y[j - 1]] + offset
      M[i, j, 1] <- logsum(c(M[i - 1, j, 1] - e, M[i - 1, j, 2] - (d + e)))#x alig to gap in y
      M[i, j, 2] <- logsum(c(M[i - 1, j - 1, 1],
                             M[i - 1, j - 1, 2],
                             M[i - 1, j - 1, 3])) + sij
      if(type == 'local' & M[i, j, 2] < 0) M[i, j, 2] <- 0
      M[i, j, 3] <- logsum(c(-Inf, M[i, j - 1, 2] - (d + e), M[i, j - 1, 3] - e))
    }
  }
  if(type == 'global'){
    # find highest score in bottom right corner of scoring array M
    score <- logsum(M[n, m, ])
  }else if(type == 'semiglobal'){
    # find highest score on bottom row or right column of scoring array M
    border <- rbind(M[n, , ], M[-n, m, ])
    ind <- which(border == max(border), arr.ind = TRUE)
    if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
    z <- if(ind[1] <= m) c(n, ind) else c(ind[1] - m, m, ind[2])
    score <- logsum(M[z[1], z[2], ])
  }else if(type == 'local'){
    # find highest score in scoring array M
    ind <- which(M[, , 2] == max(M[, , 2]), arr.ind = TRUE)
    score <- logsum(M[ind[1], ind[2], ])
  }
  res <- structure(list(score = score, array = M), class = 'forward')
  return(res)
}

