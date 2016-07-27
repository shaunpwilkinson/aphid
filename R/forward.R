#' Full log probability of sequence given model.
#'
forward.HMM <- function (x, y, logspace = FALSE){
  n <- length(y)
  states <- rownames(x$E)
  H <- length(states)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  #s <- if(logspace) x$s else log(x$s)
  if(length(y) == 0) return(structure(list(score = A[1, 1], array = NULL),
                                      class = 'forward'))
  R <- array(NA, dim = c(H, n), dimnames = list(state = states, rolls = 1:n))
  R[,1] <- E[, y[1]] + A[1, -1] # formerly s
  fun <- Vectorize(function(k, l) R[k, i - 1] + A[k + 1, l + 1])
  for (i in 2:n) R[, i] <- E[, y[i]] + apply(outer(1:H, 1:H, fun), 2, logsum)
  ak0 <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H)
  score <- logsum(R[, n] + ak0)
  res <- structure(list(score = score, array = R), class = 'forward')
  return(res)
}
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



