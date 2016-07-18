#' Implement the forward algorithm
#' 
forward.HMM <- function (x, obs, logspace = FALSE){
  n <- length(obs)
  states <- names(x$s)
  H <- length(states)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  s <- if(logspace) x$s else log(x$s)
  R <- array(NA, dim = c(H, n), dimnames = list(states = states, index = 1:n))
  R[,1] <- E[, obs[1]] + s
  fun <- Vectorize(function(k, l) R[k, i - 1] + A[k, l])
  for (i in 2:n) R[, i] <- E[, obs[i]] + apply(outer(1:H, 1:H, fun), 2, logsum)
  res <- structure(list(logFullProb = logsum(R[, n]),
                        forwardArray = R), class = 'forward')
  return(res)
}
forward.PHMM <- function (x, obs, logspace = FALSE, global = TRUE){
  E <- if(logspace) x$E - x$qe else log(x$E/x$qe)
  A <- if(logspace) x$A else log(x$A)
  n <- length(obs)
  l <- ncol(x$E)
  states <- c("M", "I", "D")
  R <- array(NA, dim = c(n + 1, l + 1, 3)) 
  dimnames(R) <- list(obs = 0:n, mod = 0:l, state = states)
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
      R[i + 1, j + 1, "M"] <- E[obs[i], j] + logsum(Mcandidates)
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


 
