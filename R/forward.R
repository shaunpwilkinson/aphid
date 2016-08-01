#' Full log likelihood.
#'
#' Forward and backward algorithms for calculating the full (log) probability
#' of a sequence given a HMM or PHMM.
#'
#' @param x an object of class \code{HMM} or \code{PHMM}, or a character vector.
#' @name forward

forward <- function(x, y, logspace = FALSE, global = TRUE){
  UseMethod("forward")
}

#' @rdname forward
forward.PHMM <- function (x, y, logspace = FALSE, global = TRUE){
  E <- if(logspace) x$E - x$qe else log(x$E/x$qe)
  A <- if(logspace) x$A else log(x$A)
  n <- length(y)
  l <- ncol(x$E)
  states <- c("M", "I", "D")
  R <- array(NA, dim = c(n + 1, l + 1, 3))
  dimnames(R) <- list(y = 0:n, mod = 0:l, state = states)
  R["0", , "M"] <- -Inf
  R[, "0", "M"] <- -Inf
  R["0", "0", "M"] <- 0
  R["0", , "I"] <- -Inf
  R["1", "0", "I"] <- R["0", "0", "M"] + A["M", "0", "I"]
  for(i in 2:n) {
    R[i + 1, "0", "I"] <- R[i, "0", "I"] + A["I", "0", "I"]
  }
  R[, "0", "D"] <- -Inf
  R["0", "1", "D"] <- R["0", "0", "M"] + A["M", "0", "D"]
  for(j in 2:l){
    R["0", j + 1, "D"] <- R["0", j, "D"] + A["D", j, "D"]
  }
  for(i in 1:n){
    for(j in 1:l){
      Mcandidates <- c(R[i, j, "M"] + A["M", j, "M"],
                       R[i, j, "I"] + A["I", j, "M"],
                       R[i, j, "D"] + A["D", j, "M"])
      Icandidates <- c(R[i, j + 1, "M"] + A["M", j + 1, "I"],
                       R[i, j + 1, "I"] + A["I", j + 1, "I"],
                       R[i, j + 1, "D"] + A["D", j + 1, "I"])
      Dcandidates <- c(R[i + 1, j, "M"] + A["M", j, "D"],
                       R[i + 1, j, "I"] + A["I", j, "D"],
                       R[i + 1, j, "D"] + A["D", j, "D"])
      R[i + 1, j + 1, "M"] <- E[y[i], j] + logsum(Mcandidates)
      R[i + 1, j + 1, "I"] <- logsum(Icandidates)
      R[i + 1, j + 1, "D"] <- logsum(Dcandidates)
    }
  }
  LLcandidates <- c(R[n + 1, l + 1, "M"] + A["M", l + 1, "M"],
                    R[n + 1, l + 1, "I"] + A["I", l + 1, "M"],
                    R[n + 1, l + 1, "D"] + A["D", l + 1, "M"])
  res <- structure(list(logFullProb = logsum(LLcandidates),
                        forwardArray = R), class = 'forward')
  return(res)
}

#' @rdname forward
forward.HMM <- function (x, y, logspace = FALSE){
  n <- length(y)
  states <- rownames(x$E)
  H <- length(states)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  #s <- if(logspace) x$s else log(x$s)
  if(length(y) == 0) structure(list(score = A[1, 1], array = NULL), class = 'forward')
  R <- array(NA, dim = c(H, n), dimnames = list(state = states, rolls = 1:n))
  R[,1] <- E[, y[1]] + A[1, -1]
  fun <- Vectorize(function(k, l) R[k, i - 1] + A[k + 1, l + 1])
  for (i in 2:n) R[, i] <- E[, y[i]] + apply(outer(1:H, 1:H, fun), 2, logsum)
  ak0 <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H)
  score <- logsum(R[, n] + ak0)
  res <- structure(list(score = score, array = R), class = 'forward')
  return(res)
}

#' @rdname forward
forward.default <- function(x, y, type = 'semiglobal', d = 8, e = 2,
                           S = NULL, itertab = NULL, offset = 0){
  ###check that x is a character vector###
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
    #alig <- trackback(z, condition = "z[1] > 1 | z[2] > 1", P = P)
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

