#' Iterative refinement of model parameters.
#'
#' Update model parameters using a list of training sequences. Methods available include
#' Viterbi training (also known as the segmental K-means algorithm (Juang & Rabiner 1990)),
#' and the Baum Welch algorithm, a special case of the expectation-maximization (EM) algorithm
#' that iteratively finds the local (but not necessarily global) optimal parameters of a HMM or PHMM.
#'
#' @param x an object of class \code{'HMM'} or \code{'PHMM'} specifying the
#' starting parameter values.
#' @param y a list of training sequences (character vectors) whose hidden states are unknown.
#' @param method a character string specifying the iterative model training method to use.
#' Must be set to either \code{method = "Viterbi"} (default) or \code{method = "BaumWelch"}.
#' @param maxiter the maximum number of EM iterations before the cycling process is terminated
#' with an error.
#' @param deltaLL a numeric value giving the change in log likelihood to specify as the convergence
#' threshold. Only applicable if \code{method = "BaumWelch"}
#' @param logspace logical argument indicating whether the emission and transition
#' probabilities of x are logged (base e; TRUE) or raw (FALSE). Alternatively, if
#' \code{logspace = "autodetect"} (default), the function will automatically detect
#' if the probabilities are in log space, returning an error if inconsistencies are found.
#' Note that choosing the latter option increases the computational
#' overhead; therefore specifying \code{TRUE} or \code{FALSE} can reduce the running time for
#' large models with many parameters.
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
#' @references Juang B-H & Rabiner L R (1990) The segmental K-means
#' algorithm for estimating parameters of hidden Markov models.
#' IEEE transactions on Acoustics,
#' Speech...
#' @references Durbin...
#' @name train
#'
train <- function(x, y, method = "Viterbi", seqweights = NULL, logspace = "autodetect",
                  maxiter = if(method == "Viterbi") 10 else 100,
                  quiet = FALSE, deltaLL = 1E-07, modelend = FALSE, pseudocounts = "Laplace",
                  gapchar = "-", fixqa = FALSE, fixqe = FALSE, inserts = "map",
                  threshold = 0.5, lambda = 0, DI = TRUE, ID = TRUE, cpp = TRUE){
  UseMethod("train")
}

#' @rdname train
train.PHMM <- function(x, y, method = "Viterbi", seqweights = NULL,
                       maxiter = if(method == "Viterbi") 10 else 100,
                       logspace = "autodetect", quiet = FALSE, deltaLL = 1E-07,
                       pseudocounts = "background", gapchar = "-",
                       fixqa = FALSE, fixqe = FALSE,
                       inserts = "map", threshold = 0.5, lambda = 0,
                       DI = TRUE, ID = TRUE, cpp = TRUE){
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  #note any changes below also need apply to align2phmm
  DNA <- is.DNA(y)
  AA <- is.AA(y)
  # DNA <- inherits(y, "DNAbin")
  gapchar <- if(DNA) as.raw(4) else if(AA) as.raw(45) else gapchar
  if(!is.list(y)){
    if(DNA | AA){
      if(is.matrix(y)){
        nseq <- nrow(y)
        seqnames <- rownames(y)
        tmp <- structure(vector(mode = "list", length = nseq), class = if(DNA) "DNAbin" else "AAbin")
        for(i in 1: nseq){
          seqi <- as.vector(y[i, ])
          tmp[[i]] <- seqi[seqi != gapchar]
        }
        names(tmp) <- seqnames
      }else{
        tmp <- structure(list(y), class = if(DNA) "DNAbin" else "AAbin")
        names(tmp) <- deparse(substitute(y))
      }
      y <- tmp
    }else if(is.matrix(y)){
      if(mode(y[1, 1]) != "character") stop("invalid mode")
      nseq <- nrow(y)
      seqnames <- rownames(y)
      tmp <- structure(vector(mode = "list", length = nseq), class = "DNAbin")
      for(i in 1: nseq){
        seqi <- as.vector(y[i, ])
        tmp[[i]] <- seqi[seqi != gapchar]
      }
      names(tmp) <- seqnames
      y <- tmp
    }else if(is.null(dim(y))){
      if(mode(y) == "character"){
        yname <- deparse(substitute(y))
        #y <- list(matrix(y, nrow = 1, dimnames = list(yname, NULL)))
        y <- list(y)
        names(y) <- yname
      }else stop("invalid mode")
    }else stop("invalid 'y' argument")
  }
  n <- length(y)
  if(is.null(seqweights)) seqweights <- rep(1, n)
  states <- c("D", "M", "I")
  residues <- rownames(x$E)
  l <- x$size
  if(!logspace){
    x$E <- log(x$E)
    x$A <- log(x$A)
  }
  if(!is.null(x$qe)){
    if(!logspace) x$qe <- log(x$qe)
  }else{
    allecs <- if(DNA){
      apply(t(sapply(y, tabulate.DNA, ambiguities = TRUE)) * seqweights, 2, sum) + 1
    }else if(AA){
      apply(t(sapply(y, tabulate.AA, ambiguities = TRUE)) * seqweights, 2, sum) + 1
    }else{
      apply(t(sapply(y, tabulate.char, residues = residues)) * seqweights, 2, sum) + 1
    } #tab(unlist(y), residues = residues) + 1 ### needs fixing
    x$qe <- log(allecs/sum(allecs))
  }
  if(!is.null(x$qa)){
    if(!logspace) x$qa <- log(x$qa)
  }else{
    alignment <- align2phmm(y, x, cpp = cpp)
    gaps <- alignment == gapchar
    inserts <- apply(gaps, 2, sum) > 0.5 * n
    xtr <- matrix(nrow = n, ncol = ncol(alignment))
    insertsn <- matrix(rep(inserts, n), nrow = n, byrow = T)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- tab9C(xtr, modules = sum(!inserts) + 2)
    transtotals <- apply(tcs, 1, sum) + 1 # force addition of Laplacian pseudos
    if(!DI) transtotals[3] <- 0
    if(!ID) transtotals[7] <- 0
    x$qa <- log(transtotals/sum(transtotals)) ### needs tidying up
  }
  out <- x
  if(method  == "Viterbi"){
    alignment <- align2phmm(y, out, logspace = TRUE, cpp = cpp)
    #scores <- attr(alignment, "score")
    #maxscore <- attr(alignment, "score")
    for(i in 1:maxiter){
      if(!quiet) cat("iteration", i, "\n")
      out <- derivePHMM(alignment, seqweights = seqweights, residues = residues,
                        inserts = inserts,
                        pseudocounts = pseudocounts, logspace = TRUE,
                        qa = if(fixqa) out$qa else NULL,
                        qe = if(fixqe) out$qe else NULL,
                        DI = DI, ID = ID) ### add others too, could just use ... = ...?
      newalig <- align2phmm(y, out, logspace = TRUE, cpp = cpp)
      #score <- attr(newalig, "score")
      #newscore <- attr(newalig, "score")
      #if(!identical(alignment, newalig) & !(score %in% scores)){
      if(!identical(alignment, newalig)){
        alignment <- newalig
        #scores <- c(scores, score)
        #maxscore <- newscore
      }else{
        if(!logspace){
          out$A <- exp(out$A)
          out$E <- exp(out$E)
          out$qa <- exp(out$qa)
          out$qe <- exp(out$qe)
        }
        if(!quiet) cat("converged after", i, "iterations\n")
        return(out)
      }
    }
    if(!quiet) cat("note: sequential alignments were not identical after", i, "iterations\n")
    return(out)
    #stop("Failed to converge. Try increasing 'maxiter' or modifying start parameters")
  }else if(method == "BaumWelch"){
    if(DNA){
      rownames(x$E) <- toupper(rownames(x$E))
      NUCorder <- sapply(rownames(x$E), match, c("A", "T", "G", "C"))
      x$E <- x$E[NUCorder, ]
      if(!(identical(rownames(x$E), c("A", "T", "G", "C")))){
        stop("Invalid model for DNA, residue alphabet does not correspond to
              nucleotide alphabet")
      }
      y <- DNA2quaternary(y, random = FALSE)
    }else if(AA){
      rownames(x$E) <- toupper(rownames(x$E))
      PFAMorder <- sapply(rownames(x$E), match, LETTERS[-c(2, 10, 15, 21, 24, 26)])
      x$E <- x$E[PFAMorder, ]
      if(!(identical(rownames(x$E), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
        stop("Invalid model residue alphabet does not correspond to
              20-letter amino acid alphabet")
      }
      y <- AA2vigesimal(y, random = FALSE)
    }else{
      y <- lapply(y, function(e) match(e, residues) - 1)
      if(any(is.na(y))) stop("Residues in sequence(s) are missing from the model")
    }
    # these just provide preformatted containers for the pseudocounts
    Apseudocounts <- x$A
    Epseudocounts <- x$E
    qepseudocounts <- x$qe
    # don't need x$qa since not updated during iteration
    if(identical(pseudocounts, "background")){
      qepseudocounts <- exp(x$qe) * length(x$qe)
      Epseudocounts[] <- rep(qepseudocounts, l)
      qacounts <- exp(x$qa) * if(DI & ID) 9 else if (DI | ID) 8 else 7
      Apseudocounts[] <- rep(qacounts, l + 1)
    }else if(identical(pseudocounts, "Laplace")){
      Apseudocounts[] <- Epseudocounts[] <- qepseudocounts[] <- 1
    }else if(identical(pseudocounts, "none")){
      Apseudocounts[] <- Epseudocounts[] <- qepseudocounts[] <- 0
    }else if(is.list(pseudocounts)){
      stopifnot(length(pseudocounts) == 3)
      stopifnot(identical(dim(pseudocounts[[1]]), dim(x$A)))
      stopifnot(identical(dim(pseudocounts[[2]]), dim(x$E)))
      stopifnot(identical(length(pseudocounts[[3]]), length(x$qe)))
      Apseudocounts[] <- pseudocounts[[1]]
      Epseudocounts[] <- pseudocounts[[2]]
      qepseudocounts[] <- pseudocounts[[3]]
    }else stop("Invalid 'pseudocounts' argument")
    Apseudocounts[1:3, 1] <- Apseudocounts[c(1, 4, 7), l + 1] <- 0
    if(!DI) Apseudocounts["DI", ] <- 0
    if(!ID) Apseudocounts["ID", ] <- 0
    E <- out$E
    A <- out$A
    qe <- out$qe
    LL <- -1E12
    for(i in 1:maxiter){
      tmpA <- Apseudocounts
      tmpE <- Epseudocounts
      tmpqe <- qepseudocounts
      tmplogPx <- rep(NA, n)
      for(j in 1:n){
        yj <- y[[j]]
        nj <- length(yj)
        if(nj == 0){
          tmpA["DD", 2:(ncol(tmpA) - 1)] <- tmpA["DD", 2:(ncol(tmpA) - 1)] + seqweights[j]
          tmplogPx[j] <- sum(c(A["MD", 1], A["DD", 2:l], A["DM", l + 1]))
        }else{
          forwj <- forward(out, yj, logspace = TRUE, odds = FALSE)
          Rj <- forwj$array
          logPxj <- forwj$score
          tmplogPx[j] <- logPxj
          backj <- backward(out, yj, logspace = TRUE, odds = FALSE)
          Bj <- backj$array
          tmpEj <- tmpE
          tmpAj <- tmpA
          tmpqej <- tmpqe
          tmpEj[] <- tmpAj[] <- tmpqej[] <- 0
          yj <- yj + 1 # R indexing style
          yjea <- cbind(FALSE, sapply(yj, function(r) 1:length(residues) == r))
          for(k in 1:l){
            # Emission and background emission counts
            for(a in seq_along(residues)){
              #yjea <- c(FALSE, yj == residues[a])
              if(any(yjea[a, ])){
                tmpEj[a, k] <- exp(logsum(Rj[k + 1, yjea[a, ], "M"] + Bj[k + 1, yjea[a, ], "M"]) - logPxj)
                tmpqej[a] <- exp(logsum(Rj[k + 1, yjea[a, ], "I"] + Bj[k + 1, yjea[a, ], "I"]) - logPxj)
              }
            }
            # Transition counts - all vectors length 1 or nj + 1 (i = 0... L)
            tmpAj["DD", k] <- exp(logsum(Rj[k, , "D"] + A["DD", k] + Bj[k + 1, , "D"]) - logPxj)
            tmpAj["MD", k] <- exp(logsum(Rj[k, , "M"] + A["MD", k] + Bj[k + 1, , "D"]) - logPxj)
            if(ID) tmpAj["ID", k] <- exp(logsum(Rj[k, , "I"] + A["ID", k] + Bj[k + 1, , "D"]) - logPxj)
            tmpAj["DM", k] <- exp(logsum(Rj[k, , "D"] + A["DM", k] + c(E[yj, k], -Inf) + c(Bj[k + 1, -1, "M"], -Inf)) - logPxj)
            tmpAj["MM", k] <- exp(logsum(Rj[k, , "M"] + A["MM", k] + c(E[yj, k], -Inf) + c(Bj[k + 1, -1, "M"], -Inf)) - logPxj)
            tmpAj["IM", k] <- exp(logsum(Rj[k, , "I"] + A["IM", k] + c(E[yj, k], -Inf) + c(Bj[k + 1, -1, "M"], -Inf)) - logPxj)
            if(DI) tmpAj["DI", k] <- exp(logsum(Rj[k, , "D"] + A["DI", k] + c(qe[yj], -Inf) + c(Bj[k, -1, "I"], -Inf)) - logPxj)
            tmpAj["MI", k] <- exp(logsum(Rj[k, , "M"] + A["MI", k] + c(qe[yj], -Inf) + c(Bj[k, -1, "I"], -Inf)) - logPxj)
            tmpAj["II", k] <- exp(logsum(Rj[k, , "I"] + A["II", k] + c(qe[yj], -Inf) + c(Bj[k, -1, "I"], -Inf)) - logPxj)
          }
          k <- l + 1
          tmpAj["DM", k] <- exp(Rj[k, nj + 1,"D"] + A["DM", k] - logPxj)
          tmpAj["MM", k] <- exp(Rj[k, nj + 1,"M"] + A["MM", k] - logPxj)
          tmpAj["IM", k] <- exp(Rj[k, nj + 1,"I"] + A["IM", k] - logPxj)
          if(DI) tmpAj["DI", k] <- exp(logsum(Rj[k, , "D"] + A["DI", k] + c(qe[yj], -Inf) + c(Bj[k, -1, "I"], -Inf)) - logPxj)
          tmpAj["MI", k] <- exp(logsum(Rj[k, , "M"] + A["MI", k] + c(qe[yj], -Inf) + c(Bj[k, -1, "I"], -Inf)) - logPxj)
          tmpAj["II", k] <- exp(logsum(Rj[k, , "I"] + A["II", k] + c(qe[yj], -Inf) + c(Bj[k, -1, "I"], -Inf)) - logPxj)
          # correct for sequence weight
          tmpA <- tmpA + tmpAj * seqweights[j]
          tmpE <- tmpE + tmpEj * seqweights[j]
          tmpqe <- tmpqe + tmpqej * seqweights[j]
        }
      }
      tmpE <- log(tmpE/apply(tmpE, 1, sum))
      tmpqe <- log(tmpqe/sum(tmpqe))
      tmpA <- t(tmpA)
      for(X in c(1, 4, 7)) tmpA[, X:(X + 2)] <- log(tmpA[, X:(X + 2)]/apply(tmpA[, X:(X + 2)], 1, sum))
      tmpA <- t(tmpA)
      tmpA[1:3, 1] <- -Inf # replace NaNs generated when dividing by 0
      A <- tmpA
      E <- tmpE
      if(!fixqe) qe <- tmpqe
      out$A <- A
      out$E <- E
      if(!fixqe) out$qe <- qe
      ### some cleanup required above
      logPx <- sum(tmplogPx) # page 62 eq 3.17
      if(!quiet) cat("Iteration", i, "log likelihood =", logPx, "\n")
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
    stop("Failed to converge. Try increasing 'maxiter' or modifying start parameters")
  }else stop("Invalid argument given for 'method'")

}

#' @rdname train
train.HMM <- function(x, y, method = "Viterbi", seqweights = NULL,
                      maxiter = if(method == "Viterbi") 10 else 100,
                      deltaLL = 1E-07,
                      logspace = "autodetect", quiet = FALSE, modelend = FALSE,
                      pseudocounts = "Laplace", cpp = TRUE){
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  if(is.list(y)){
  } else if(is.vector(y, mode = "character")){
    y <- list(y)
  }else stop("invalid y argument")
  n <- length(y)
  if(is.null(seqweights)) seqweights <- rep(1, n)
  stopifnot(sum(seqweights) == n)
  states <- rownames(x$A)
  nstates <- length(states)
  residues <- colnames(x$E)
  nres <- length(residues)
  out <- x
  if(!logspace){
    out$E <- log(out$E)
    out$A <- log(out$A)
  }
  if(method == "Viterbi"){
    for(i in 1:maxiter){
      samename <- logical(n)
      for(j in 1:n){
        vitj <- Viterbi(out, y[[j]], logspace = TRUE, cpp = cpp)
        pathchar <- states[-1][vitj$path + 1]
        if(identical(pathchar, names(y[[j]]))) samename[j] <- TRUE
        # need to be named to feed into deriveHMM
        names(y[[j]]) <- pathchar
      }
      if(all(samename)){
        if(!logspace){
          out$A <- exp(out$A)
          out$E <- exp(out$E)
        }
        if(!quiet) cat("iteration", i, "\nconverged after", i, "iterations\n")
        return(out)
      }else{
        if(!quiet) cat("iteration", i, "\n")
        out <- deriveHMM(y, seqweights = seqweights, residues = residues,
                         states = states, modelend = modelend,
                         pseudocounts = pseudocounts, logspace = TRUE)
      }
    }
    stop("Failed to converge. Try increasing 'maxiter' or modifying start parameters")
  }else if(method == "BaumWelch"){
    Apseudocounts <- matrix(0, nrow = nstates, ncol = nstates)
    Epseudocounts <- matrix(0, nrow = nstates - 1, ncol = nres)
    dimnames(Apseudocounts) <- list(from = states, to =  states)
    dimnames(Epseudocounts) <- list(state = states[-1], residue = residues)
    if(identical(pseudocounts, "Laplace")){
      Apseudocounts[] <- Epseudocounts[] <- 1
      if(!modelend) Apseudocounts[, 1] <- 0
    } else if(is.list(pseudocounts)){
      stopifnot(length(pseudocounts) == 2)
      stopifnot(identical(dim(pseudocounts[[1]]), dim(x$A)))
      stopifnot(identical(dim(pseudocounts[[2]]), dim(x$E)))
      Apseudocounts[] <- pseudocounts[[1]]
      Epseudocounts[] <- pseudocounts[[2]]
    } else if(!identical(pseudocounts, "none")) stop("invalid 'pseudocounts' argument")
    E <- out$E
    A <- out$A
    LL <- -1E12
    for(i in 1:maxiter){
      tmpA <- Apseudocounts
      tmpE <- Epseudocounts
      tmplogPx <- rep(NA, n)
      for(j in 1:n){
        yj <- y[[j]]
        nj <- length(yj)
        if(nj == 0){
          tmpA[1, 1] <- tmpA[1, 1] + if(modelend) seqweights[j] else 0
        }else{
          forwj <- forward(out, yj, logspace = TRUE)
          Rj <- forwj$array
          logPxj <- forwj$score
          tmplogPx[j] <- logPxj
          backj <- backward(out, yj, logspace = TRUE)
          Bj <- backj$array
          tmpAj <- tmpA
          tmpEj <- tmpE
          tmpAj[] <- tmpAj[] <- 0
          for(k in states[-1]){
            tmpAj[1, -1] <- exp(A[1, -1] + E[, yj[1]] + Bj[, 1] - logPxj)
            tmpAj[-1, 1] <- exp(Rj[, nj] + A[-1, 1] - logPxj)
            for(l in states[-1]){
              tmpAj[k, l] <- exp(logsum(Rj[k, -nj] + A[k, l] + E[l, yj[-1]] + Bj[l, -1]) - logPxj)
            }
            for(b in residues){
              cond <- yj == b
              tmpEj[k, b] <- exp(logsum(Rj[k, cond] + Bj[k, cond]) - logPxj)
            }
          }
          # correct for sequence weight
          tmpA <- tmpA + tmpAj * seqweights[j]
          tmpE <- tmpE + tmpEj * seqweights[j]
        }
      }
      A[] <- log(tmpA/apply(tmpA, 1, sum))
      E[] <- log(tmpE/apply(tmpE, 1, sum))
      out$A <- A
      out$E <- E
      logPx <- sum(tmplogPx) # page 62 eq 3.17
      if(!quiet) cat("iteration", i, "log likelihood = ", logPx, "\n")
      if(abs(LL - logPx) < deltaLL){
        if(!logspace){
          out$A <- exp(out$A)
          out$E <- exp(out$E)
        }
        if(!quiet) cat("convergence threshold reached after", i, "EM iterations\n")
        return(out)
      }
      LL <- logPx
    }
    stop("Failed to converge on a local maximum. Try increasing 'maxiter',
        decreasing 'deltaLL' or modifying start parameters")
  }else stop("invalid argument given for 'method'")
}
