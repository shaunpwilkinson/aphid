#' Iterative estimation of model parameters.
#'
#' \code{BaumWelch} uses a special case of the expectation-maximization (EM) algorithm
#' to find the locally-optimal parameters of a HMM or PHMM.
#' @param x an object of class \code{'HMM'} or \code{'PHMM'} specifying the
#' starting parameter values.
#' @param y a list of training sequences (character vectors) whose hidden states are unknown.
#' @param maxiter the maximum number of EM iterations before the process is terminated
#' @param deltaLL the change in log likelihood specified as the convergence threshold
#' @param logspace logical argument indicating whether the emission and transition
#' probabilities of x are logged (base e; TRUE) or raw (FALSE). Alternatively, if
#' \code{logspace = "autodetect"} (default), the function will automatically detect
#' if the probabilities are in log space, returning an error if inconsistencies are found.
#' Note that choosing the latter option increases the computational
#' overhead; therefore specifying \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param quiet logical argument indicating whether the iteration progress should
#' be printed to the console (\code{quiet = FALSE}; default) or suppressed
#' (\code{quiet = TRUE}).
#' @param pseudocounts used to account for the possible absence of certain transition
#' and/or emission types in the training dataset.
#' either \code{'Laplace'} (adds one of each possible transition and emission type to the
#' training dataset; default), \code{'none'}, or a two-element list containing a matrix of
#' transition pseudocounts as its first element and a matrix of emission pseudocounts
#' as its second. If a list is supplied both matrices must have row and column names
#' according to the residues (column names of emission matrix) and states
#' (row and column names of the transition matrix and row names of the emission matrix).
#' The first row and column of the transition matrix must be 'BeginEnd'. Background
#' pseudocounts are recommended for small training sets,
#' since Laplacian counts can overinflate insert and delete transition probabilities
#' leading to convergence at suboptimal local maxima.
#' @seealso \code{\link{deriveHMM}} and \code{\link{derivePHMM}} for
#' maximum-likelihood parameter estimation when training sequence states are
#' known.
#' @name BaumWelch
#'
BaumWelch <- function(x, y, maxiter = 100, deltaLL = 1E-07,
                          logspace = "autodetect", quiet = FALSE,
                          pseudocounts = "background", DI = TRUE){
  UseMethod("BaumWelch")
}


#' @rdname BaumWelch
BaumWelch.PHMM <- function(x, y, maxiter = 100, deltaLL = 1E-07,
                           logspace = "autodetect", quiet = FALSE,
                           pseudocounts = "background", DI = TRUE){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  if(is.list(y)){
  }else if(is.vector(y, mode = "character")){
    y <- list(y)
  }else{
    stop("invalid y argument")
  }
  ### need qa, qe etc check function(valid?)

  n <- length(y)
  residues <- rownames(x$E)
  states <- c("D", "M", "I")

  if(!logspace){
    x$E <- log(x$E)
    x$A <- log(x$A)
  }

  if(!is.null(x$qe)){
    if(!logspace) x$qe <- log(x$qe)
  }else{
    allecs <- tab(unlist(y), residues = residues) + 1
    x$qe <- log(allecs/sum(allecs))
  }

  if(!is.null(x$qa)){
    if(!logspace) x$qa <- log(x$qa)
  } else{
    alignment <- align2phmm(y, x, gapchar = "-")
    gaps <- alignment == "-"
    inserts <- apply(gaps, 2, sum) > 0.5 * n
    xtr <- matrix(nrow = n, ncol = ncol(alignment))
    insertsn <- matrix(rep(inserts, n), nrow = n, byrow = T)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- tab9C(xtr, modules = sum(!inserts) + 2)
    transtotals <- apply(tcs, 1, sum) + 1 # forced addition of Laplacian pseudos
    if(!DI) transtotals[c(3, 7)] <- 0
    qa <- matrix(transtotals, nrow = 3, byrow = TRUE)
    dimnames(qa) <- list(from = c("D", "M", "I"), to = c("D", "M", "I"))
    # qa <- qa/apply(qa, 1, sum)
    x$qa <- log(qa/sum(qa)) ### needs tidying up
  }


  # these just provide preformatted structures for the pseudocounts
  Apseudocounts <- x$A
  Epseudocounts <- x$E
  qepseudocounts <- x$qe
  # don't need qa since not updated

  if(identical(pseudocounts, "background")){
    qepseudocounts <- exp(x$qe) * length(x$qe)
    Epseudocounts[] <- rep(qepseudocounts, x$size)
    qacounts <- exp(x$qa) * if(DI) 9 else 7
    Apseudocounts[,,1] <- c(rep(qacounts[,1], x$size), rep(0, 3))
    Apseudocounts[,,2] <- rep(qacounts[,2], x$size + 1)
    Apseudocounts[,,3] <- rep(qacounts[,3], x$size + 1)
    Apseudocounts[1, 1, ] <- 0
  } else if(identical(pseudocounts, "Laplace")){
    Apseudocounts[] <- Epseudocounts[] <- qepseudocounts[] <- 1
    Apseudocounts["D", "0", ] <- Apseudocounts[ , ncol(x$A), "D"] <- rep(0, 3)
  } else if(identical(pseudocounts, "none")){
    Apseudocounts[] <- Epseudocounts[] <- qepseudocounts[] <- 0
  } else if(is.list(pseudocounts)){
    stopifnot(length(pseudocounts == 3))
    Apseudocounts[] <- pseudocounts[[1]]
    ### check for DI conflict
    Epseudocounts[] <- pseudocounts[[2]]
    qepseudocounts[] <- pseudocounts[[3]]
  } else stop("invalid 'pseudocounts' argument")
  if(!DI) Apseudocounts["D", , "I"] <- Apseudocounts["I", , "D"] <- 0

  out <- x
  E <- out$E
  A <- out$A
  qe <- out$qe
  LL <- -1E06
  for(i in 1:maxiter){
    tmpA <- Apseudocounts
    tmpE <- Epseudocounts
    tmpqe <- qepseudocounts
    tmplogPx <- rep(NA, n)
    for(j in 1:n){
      yj <- y[[j]]
      nj <- length(yj)
      if(nj == 0){
        tmpA["D", 2:(ncol(tmpA) - 1), "D"] <- tmpA["D", 2:(ncol(tmpA) - 1), "D"] + 1
        tmplogPx[j] <- sum(c(A["M", 1, "D"],
                             A["D", 2:(ncol(A[, , "D"]) - 1) , "D"],
                             A["D", ncol(A[, , "D"]), "M"]))
      }else{
        forwj <- forward(out, yj, logspace = TRUE, odds = FALSE)
        Rj <- forwj$array
        logPxj <- forwj$score
        tmplogPx[j] <- logPxj
        backj <- backward(out, yj, logspace = TRUE, odds = FALSE)
        Bj <- backj$array
        for(k in 1:ncol(x$E)){ #modules
          for(a in seq_along(residues)){  #residue alphabet
            yjequalsa <- c(FALSE, yj == residues[a])
            if(any(yjequalsa)){
              tmpE[a, k] <- tmpE[a, k] + exp(logsum(Rj[k + 1, yjequalsa, "M"] +
                                                      Bj[k + 1, yjequalsa, "M"]) -
                                               logPxj)
              tmpqe[a] <- tmpqe[a] + exp(logsum(Rj[k + 1, yjequalsa, "I"] +
                                                  Bj[k + 1, yjequalsa, "I"]) -
                                           logPxj)
            }
          }
        }
        for(k in 1:(ncol(A[,,"M"]) - 1)){
          # all vectors length 1 or nj + 1 ( i = 0... L)
          tmpA["D", k, "D"] <- tmpA["D", k, "D"] +
            exp(logsum(Rj[k, , "D"] + A["D", k, "D"] + Bj[k + 1, ,"D"]) - logPxj) #d to d
          tmpA["M", k, "D"] <- tmpA["M", k, "D"] +
            exp(logsum(Rj[k, , "M"] + A["M", k, "D"] + Bj[k + 1, ,"D"]) - logPxj) #m to d
          tmpA["I", k, "D"] <- tmpA["I", k, "D"] +
            exp(logsum(Rj[k, , "I"] + A["I", k, "D"] + Bj[k + 1, ,"D"]) - logPxj) #i to d
          tmpA["D", k, "M"] <- tmpA["D", k, "M"] +
            exp(logsum(Rj[k, , "D"] + A["D", k, "M"] + c(E[yj, k], -Inf) +
                         c(Bj[k + 1, -1,"M"], -Inf)) - logPxj) #d to m
          tmpA["M", k, "M"] <- tmpA["M", k, "M"] +
            exp(logsum(Rj[k, , "M"] + A["M", k, "M"] + c(E[yj, k], -Inf) +
                         c(Bj[k + 1, -1,"M"], -Inf)) - logPxj) #m to m
          tmpA["I", k, "M"] <- tmpA["I", k, "M"] +
            exp(logsum(Rj[k, , "I"] + A["I", k, "M"] + c(E[yj, k], -Inf) +
                         c(Bj[k + 1, -1,"M"], -Inf)) - logPxj) #i to m
          tmpA["D", k, "I"] <- tmpA["D", k, "I"] +
            exp(logsum(Rj[k, , "D"] + A["D", k, "I"] + c(qe[yj], -Inf) +
                         c(Bj[k, -1, "I"], -Inf)) - logPxj) #d to i
          tmpA["M", k, "I"] <- tmpA["M", k, "I"] +
            exp(logsum(Rj[k, , "M"] + A["M", k, "I"] + c(qe[yj], -Inf) +
                         c(Bj[k, -1, "I"], -Inf)) - logPxj) #m to i
          tmpA["I", k, "I"] <- tmpA["I", k, "I"] +
            exp(logsum(Rj[k, , "I"] + A["I", k, "I"] + c(qe[yj], -Inf) +
                         c(Bj[k, -1, "I"], -Inf)) - logPxj) #i to i
        }
        k <- ncol(A[,,"M"])
        tmpA["D", k, "M"] <- tmpA["D", k, "M"] +
          exp(Rj[k, nj + 1,"D"] + A["D", k, "M"] - logPxj)
        tmpA["M", k, "M"] <- tmpA["M", k, "M"] +
          exp(Rj[k, nj + 1,"M"] + A["M", k, "M"] - logPxj)
        tmpA["I", k, "M"] <- tmpA["I", k, "M"] +
          exp(Rj[k, nj + 1,"I"] + A["I", k, "M"] - logPxj)
        tmpA["D", k, "I"] <- tmpA["D", k, "I"] +
          exp(logsum(Rj[k, , "D"] + A["D", k, "I"] + c(qe[yj], -Inf) +
                       c(Bj[k, -1, "I"], -Inf)) - logPxj) #d to i
        tmpA["M", k, "I"] <- tmpA["M", k, "I"] +
          exp(logsum(Rj[k, , "M"] + A["M", k, "I"] + c(qe[yj], -Inf) +
                       c(Bj[k, -1, "I"], -Inf)) - logPxj) #m to i
        tmpA["I", k, "I"] <- tmpA["I", k, "I"] +
          exp(logsum(Rj[k, , "I"] + A["I", k, "I"] + c(qe[yj], -Inf) +
                       c(Bj[k, -1, "I"], -Inf)) - logPxj)
      }
    }
    tmpE <- log(tmpE/apply(tmpE, 1, sum))
    tmpqe <- log(tmpqe/sum(tmpqe))
    for(X in states) tmpA[X, , ] <- log(tmpA[X, , ]/apply(tmpA[X, , ], 1, sum))
    tmpA["D", "0", ] <- rep(-Inf, 3) # replace NaNs generated when dividing by 0
    A <- tmpA
    E <- tmpE
    qe <- tmpqe
    out$A <- A
    out$E <- E
    out$qe <- qe
    logPx <- sum(tmplogPx) # page 62 eq 3.17
    if(!quiet) cat("Iteration", i, "log likelihood = ", logPx, "\n")
    if(abs(LL - logPx) < deltaLL){
      if(!logspace){
        out$A <- exp(out$A)
        out$E <- exp(out$E)
        out$qe <- exp(out$qe)
      }
      if(!quiet) cat("Convergence threshold reached after", i, "EM iterations\n")
      return(out)
    }
    LL <- logPx
  }
}


#' @rdname BaumWelch
BaumWelch.HMM <- function(x, y, maxiter = 100, deltaLL = 1E-07,
                      logspace = "autodetect", quiet = FALSE,
                      modelend = FALSE, pseudocounts = "Laplace"){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  if(is.list(y)){
  } else if(is.vector(y, mode = "character")){
    y <- list(y)
  } else stop("invalid y argument")
  n <- length(y)
  states <- rownames(x$A)
  nstates <- length(states)
  residues <- colnames(x$E)
  nres <- length(residues)
  Apseudocounts <- matrix(nrow = nstates, ncol = nstates)
  Epseudocounts <- matrix(nrow = nstates - 1, ncol = nres)
  dimnames(Apseudocounts) <- list(from = states, to =  states)
  dimnames(Epseudocounts) <- list(state = states[-1], residue = residues)
  if(identical(pseudocounts, "Laplace")){
    Apseudocounts[] <- Epseudocounts[] <- 1
    if(!modelend) Apseudocounts[, 1] <- 0
  } else if(identical(pseudocounts, "none")){
    AApseudocounts[] <- Epseudocounts[] <- 0
  } else if(is.list(pseudocounts)){
    stopifnot(length(pseudocounts == 2))
    Apseudocounts[] <- pseudocounts[[1]]
    Epseudocounts[] <- pseudocounts[[2]]
  } else stop("invalid 'pseudocounts' argument")
  out <- x
  if(!logspace){
    out$E <- log(out$E)
    out$A <- log(out$A)
  }
  E <- out$E
  A <- out$A
  LL <- -1E06
  for(i in 1:maxiter){
    tmpA <- Apseudocounts
    tmpE <- Epseudocounts
    tmplogPx <- rep(NA, n)
    for(j in 1:n){
      yj <- y[[j]]
      nj <- length(yj)
      if(nj == 0){
        tmpA[1, 1] <- tmpA[1, 1] + if(modelend) 1 else 0
      }else{
        forwj <- forward(out, yj, logspace = TRUE)
        Rj <- forwj$array
        logPxj <- forwj$score
        tmplogPx[j] <- logPxj
        backj <- backward(out, yj, logspace = TRUE)
        Bj <- backj$array
        for(k in states[-1]){
          tmpA[1, -1] <- tmpA[1, -1] + exp(A[1, -1] + E[, yj[1]] + Bj[, 1] - logPxj)
          tmpA[-1, 1] <- tmpA[-1, 1] + exp(Rj[, nj] + A[-1, 1] - logPxj)
          for(l in states[-1]){
            tmpA[k, l] <- tmpA[k, l] + exp(logsum(Rj[k, -nj] + A[k, l] +
                                                    E[l, yj[-1]] + Bj[l, -1])
                                           - logPxj)
          }
          for(b in residues){
            cond <- yj == b
            tmpE[k, b] <- tmpE[k, b] + exp(logsum(Rj[k, cond] + Bj[k, cond]) - logPxj)
          }
        }
      }
    }
    A[] <- log(tmpA/apply(tmpA, 1, sum))
    E[] <- log(tmpE/apply(tmpE, 1, sum))
    out$A <- A
    out$E <- E
    logPx <- sum(tmplogPx) # page 62 eq 3.17
    if(!quiet) cat("Iteration", i, "log likelihood = ", logPx, "\n")
    if(abs(LL - logPx) < deltaLL){
      if(!logspace){
        out$A <- exp(out$A)
        out$E <- exp(out$E)
      }
      if(!quiet) cat("Convergence threshold reached after", i, "EM iterations\n")
      return(out)
    }
    LL <- logPx
  }
  stop("Failed to converge on a local maximum. Try increasing 'maxiter',
        decreasing 'deltaLL' or modifying start parameters")
}


