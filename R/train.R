#' Iterative model refinement.
#'
#' Update model parameters using a list of training sequences, with either
#'   the Viterbi training or Baum-Welch algorithm.
#'
#' @param x an object of class \code{'HMM'} or \code{'PHMM'} specifying the
#'   initial parameter values.
#' @param y a list of training sequences whose hidden states are unknown.
#'   Accepted modes are "character" and "raw" (for "DNAbin" and "AAbin"
#'   objects).
#' @param method a character string specifying the iterative model training
#'   method to use. Accepted methods are \code{"Viterbi"} (the default)
#'   and \code{"BaumWelch"}.
#' @param seqweights either NULL (default; all sequences are given
#'   weights of 1) or a numeric vector the same length as \code{x}
#'   representing the sequence weights used to derive the model.
#' @param logspace logical indicating whether the emission and transition
#'   probabilities of x are logged. If \code{logspace = "autodetect"}
#'   (default setting), the function will automatically detect
#'   if the probabilities are logged, returning an error if
#'   inconsistencies are found. Note that choosing the latter option
#'   increases the computational overhead; therefore specifying
#'   \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param maxiter the maximum number of EM iterations or Viterbi training
#'   iterations to carry out before the cycling process is terminated and
#'   the partially trained model is returned. Defaults to 100.
#' @param deltaLL numeric, the maximum change in log likelihood between EM
#'   iterations before the cycling procedure is terminated (signifying model
#'   convergence). Defaults to 1E-07. Only applicable if
#'   \code{method = "BaumWelch"}.
#' @param modelend logical indicating whether transition probabilites
#'   to the end state of the standard hidden Markov model should be
#'   modeled (if applicable). Defaults to FALSE.
#' @param pseudocounts character string, either "background", Laplace"
#'   or "none". Used to account for the possible absence of certain
#'   transition and/or emission types in the input sequences.
#'   If \code{pseudocounts = "background"} (default), pseudocounts
#'   are calculated from the background transition and emission
#'   frequencies in the training dataset.
#'   If \code{pseudocounts = "Laplace"} one of each possible transition
#'   and emission type is added to the training dataset (default).
#'   If \code{pseudocounts = "none"} no pseudocounts are added (not
#'   usually recommended, since low frequency transition/emission types
#'   may be excluded from the model).
#'   Alternatively this argument can be a two-element list containing
#'   a matrix of transition pseudocounts
#'   as its first element and a matrix of emission pseudocounts as its
#'   second. If this option is selected, both matrices must have row and column
#'   names corresponding with the residues (column names of emission matrix)
#'   and states (row and column names of the transition matrix and
#'   row names of the emission matrix). For standard HMMs
#'   the first row and column of the transition matrix should be named
#'   "Begin".
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param fixqa logical. Should the background transition probabilities
#'   be fixed (TRUE), or allowed to vary between iterations (FALSE)?
#'   Defaults to FALSE. Only applicable if \code{method = "Viterbi"}.
#' @param fixqe logical. Should the background emission probabilities
#'   be fixed (TRUE), or allowed to vary between iterations (FALSE)?
#'   Defaults to FALSE. Only applicable if \code{method = "Viterbi"}.
#' @param maxsize integer giving the upper bound on the number of modules
#'   in the PHMM. If NULL (default) no maximum size is enforced.
#' @param inserts character string giving the model construction method
#'   by which alignment columns
#'   are marked as either match or insert states. Accepted methods include
#'   \code{"threshold"} (only columns with fewer than a specified
#'   proportion of gaps form match states in the model), \code{"map"} (default;
#'   match and insert columns are found using the maximum \emph{a posteriori}
#'   method outlined in Durbin et al. (1998) chapter 5.7), \code{"inherited"}
#'   (match and insert columns are inherited from the "inserts" attribute
#'   of the input alignment), and \code{"none"} (all columns are assigned
#'   match states in the model). Alternatively, insert columns can be
#'   specified manually by providing a logical vector the same length
#'   as the number of columns in the alignment, with \code{TRUE} for insert
#'   columns and \code{FALSE} for match states.
#' @param threshold the maximum proportion of gaps for an alignment column
#'   to be considered for a match state in the PHMM (defaults to 0.5).
#'   Only applicable when \code{inserts = "threshold"}.
#'   Note that the maximum \emph{a posteriori}
#'   method works poorly for alignments with few sequences,
#'   so the 'threshold' method is
#'   automatically used when the number of sequences is less than 5.
#' @param lambda penalty parameter used to favour models with fewer match
#'   states. Equivalent to the log of the prior probability of marking each
#'   column (Durbin et al. 1998, chapter 5.7). Only applicable when
#'   \code{inserts = "map"}.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @param ... aditional arguments to be passed to \code{"Viterbi"} (if
#'   \code{method = "Viterbi"}) or \code{"forward"} (if
#'   \code{method = "BaumWelch"}).
#' @return an object of class \code{"HMM"} or \code{"PHMM"}, depending
#'   on the input model \code{x}.
#' @details
#'   This function offers a choice of two model training methods,
#'   Viterbi training (also known as the segmental
#'   K-means algorithm (Juang & Rabiner 1990)), and the Baum Welch algorithm,
#'   a special case of the expectation-maximization (EM) algorithm that
#'   iteratively finds the locally (but not necessarily globally) optimal
#'   parameters of a HMM or PHMM.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Juang B-H, Rabiner LR (1990) The segmental K-means
#'   algorithm for estimating parameters of hidden Markov models.
#'   \emph{IEEE Transactions on Acoustics, Speech, and Signal Processing},
#'   \strong{38}, 1639-1641.
#'
#' @seealso \code{\link{deriveHMM}} and \code{\link{derivePHMM}} for
#'   maximum-likelihood parameter estimation when the training sequence states
#'   are known.
#' @examples
#'   ## Baum Welch training for standard HMMs
#'   ## The dishonest casino example from Durbin et al. (1998) chapter 3.2
#'   A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9),
#'               nrow = 3) # transition probability matrix
#'   dimnames(A) <- list(from = c("Begin", "Fair", "Loaded"),
#'                       to = c("Begin", "Fair", "Loaded"))
#'   E <- matrix(c((1/6), (1/6), (1/6), (1/6), (1/6), (1/6),
#'                 (1/10), (1/10), (1/10), (1/10), (1/10), (1/2)),
#'               nrow = 2, byrow = TRUE) # emission probability matrix
#'   dimnames(E) <- list(states = c('Fair', 'Loaded'), residues = paste(1:6))
#'   x <- structure(list(A = A, E = E), class = "HMM") # create hidden Markov model
#'   par(mfrow = c(2, 1))
#'   plot(x, main = "Dishonest casino HMM before training")
#'   data(casino)
#'   x <- train(x, list(casino), method = "BaumWelch", deltaLL = 0.001)
#'   plot(x, main = "Dishonest casino HMM after training")
#' @name train
################################################################################
# train <- function(x, y, method = "Viterbi", seqweights = NULL,
#                   logspace = "autodetect", maxiter = 100,
#                   deltaLL = 1E-07, modelend = FALSE,
#                   pseudocounts = "background", gap = "-",
#                   fixqa = FALSE, fixqe = FALSE, maxsize = NULL,
#                   inserts = "map", threshold = 0.5, lambda = 0,
#                   quiet = FALSE, ...){
#   UseMethod("train")
# }
train <- function(x, y, ...){
  UseMethod("train")
}
################################################################################
#' @rdname train
################################################################################
train.PHMM <- function(x, y, method = "Viterbi", seqweights = NULL,
                       logspace = "autodetect", maxiter = 100, deltaLL = 1E-07,
                       pseudocounts = "background", gap = "-",
                       fixqa = FALSE, fixqe = FALSE, maxsize = NULL,
                       inserts = "map", threshold = 0.5, lambda = 0,
                       quiet = FALSE, ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  #note any changes below also need apply to align2phmm
  DNA <- .isDNA(y)
  AA <- .isAA(y)
  DI <- !all(x$A["DI", ] == if(logspace) -Inf else 0)
  ID <- !all(x$A["ID", ] == if(logspace) -Inf else 0)
  maxiter <- maxiter
  #method <- toupper(method)
  gap <- if(DNA) as.raw(4) else if(AA) as.raw(45) else gap
  if(!is.list(y)){
    if(DNA | AA){
      if(is.matrix(y)){
        nseq <- nrow(y)
        seqnames <- rownames(y)
        tmp <- structure(vector(mode = "list", length = nseq), class = if(DNA) "DNAbin" else "AAbin")
        for(i in 1: nseq){
          seqi <- as.vector(y[i, ])
          tmp[[i]] <- seqi[seqi != gap]
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
        tmp[[i]] <- seqi[seqi != gap]
      }
      names(tmp) <- seqnames
      y <- tmp
    }else if(is.null(dim(y))){
      if(mode(y) == "character"){
        yname <- deparse(substitute(y))
        y <- list(y)
        names(y) <- yname
      }else stop("invalid mode")
    }else stop("invalid 'y' argument")
  }
  n <- length(y)
  if(is.null(seqweights)) seqweights <- rep(1, n)
  ##### TODO else (seqweights used as numeric below)
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
      apply(t(sapply(y, .tabulateDNA, ambiguities = TRUE)) * seqweights, 2, sum) + 1
    }else if(AA){
      apply(t(sapply(y, .tabulateAA, ambiguities = TRUE)) * seqweights, 2, sum) + 1
    }else{
      apply(t(sapply(y, .tabulateCH, residues = residues)) * seqweights, 2, sum) + 1
    }
    x$qe <- log(allecs/sum(allecs))
  }
  if(!is.null(x$qa)){
    if(!logspace) x$qa <- log(x$qa)
  }else{
    alignment <- align(y, x, ... = ...)
    gaps <- alignment == gap
    inserts <- apply(gaps, 2, sum) > 0.5 * n
    xtr <- matrix(nrow = n, ncol = ncol(alignment))
    insertsn <- matrix(rep(inserts, n), nrow = n, byrow = T)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- .atab(xtr, seqweights = seqweights)
    transtotals <- apply(tcs, 1, sum) + 1 # force addition of Laplacian pseudos
    if(!DI) transtotals[3] <- 0
    if(!ID) transtotals[7] <- 0
    x$qa <- log(transtotals/sum(transtotals)) ### needs tidying up
  }
  if(method  == "Viterbi"){
    alignment <- align(y, model = x, logspace = TRUE, ... = ...)
    #scores <- attr(alignment, "score")
    #maxscore <- attr(alignment, "score")
    alig_cache <- character(maxiter)
    alig_cache[1] <- paste(openssl::md5(as.vector(alignment)))
    # alig_cache <- list()
    # alig_cache[[1]] <- as.vector(alignment)
    for(i in 1:maxiter){
      out <- derivePHMM.default(alignment, seqweights = seqweights, residues = residues,
                                gap = gap, DI = DI, ID = ID, maxsize = maxsize,
                        inserts = inserts, lambda = lambda, threshold = threshold,
                        pseudocounts = pseudocounts, logspace = TRUE,
                        qa = if(fixqa) x$qa else NULL,
                        qe = if(fixqe) x$qe else NULL)
      if(!quiet) cat("Iteration", i, "PHMM with", out$size, "internal modules\n")
      ### what about DI and ID
      newalig <- align(y, model = out, logspace = TRUE, ... = ...)
      # newaligv <- as.vector(newalig)
      newaligv <- paste(openssl::md5(as.vector(newalig)))
      if(!any(sapply(alig_cache, identical, newaligv))){
        # alig_cache[[i + 1]] <- newaligv
        alig_cache[i + 1] <- newaligv
      #score <- attr(newalig, "score")
      #newscore <- attr(newalig, "score")
      #if(!identical(alignment, newalig) & !(score %in% scores)){
      #if(!identical(as.vector(alignment), as.vector(newalig))){
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
        if(!quiet) cat("Sequential alignments were identical after", i, "iterations\n")
        return(out)
      }
    }
    if(!quiet) cat("Note: sequential alignments were not identical after", i, "iterations\n")
    return(out)
    #stop("Failed to converge. Try increasing 'maxiter' or modifying start parameters")
  }else if(method == "BaumWelch"){
    if(DNA){
      NUCorder <- sapply(toupper(rownames(x$E)), match, c("A", "T", "G", "C"))
      x$E <- x$E[NUCorder, ]
      x$qe <- x$qe[NUCorder]
      if(!(identical(toupper(rownames(x$E)), c("A", "T", "G", "C")))){
        stop("Invalid model for DNA, residue alphabet does not correspond to
              nucleotide alphabet")
      }
      y <- .encodeDNA(y, arity = 4, probs = exp(x$qe), random = FALSE, na.rm = TRUE)
    }else if(AA){
      #rownames(x$E) <- toupper(rownames(x$E))
      PFAMorder <- sapply(toupper(rownames(x$E)), match, LETTERS[-c(2, 10, 15, 21, 24, 26)])
      x$E <- x$E[PFAMorder, ]
      x$qe <- x$qe[PFAMorder]
      if(!(identical(toupper(rownames(x$E)), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
        stop("Invalid model for AA, residue alphabet does not correspond to
              20-letter amino acid alphabet")
      }
      y <- .encodeAA(y, arity = 20, probs = exp(x$qe), random = FALSE, na.rm = TRUE)
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
    E <- x$E
    A <- x$A
    qe <- x$qe
    LL <- -1E12
    out <- x
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
          forwj <- forward(out, yj, logspace = TRUE, odds = FALSE, ... = ...)
          Rj <- forwj$array
          logPxj <- forwj$score
          tmplogPx[j] <- logPxj
          backj <- backward(out, yj, logspace = TRUE, odds = FALSE, ... = ...)
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
      tmpE <- t(tmpE)
      tmpE <- log(tmpE/apply(tmpE, 1, sum))
      tmpE <- t(tmpE)
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
        if(DNA){
          out$E <- out$E[NUCorder, ]
          out$qe <- out$qe[NUCorder]
        }else if(AA){
          out$E <- out$E[PFAMorder, ]
          out$qe <- out$qe[PFAMorder]
        }
        if(!quiet) cat("Convergence threshold reached after", i, "EM iterations\n")
        return(out)
      }
      LL <- logPx
    }
    warning("Failed to converge. Try increasing 'maxiter' or modifying start parameters")
    if(!logspace){
      out$A <- exp(out$A)
      out$E <- exp(out$E)
      out$qe <- exp(out$qe)
    }
    if(DNA){
      out$E <- out$E[NUCorder, ]
      out$qe <- out$qe[NUCorder]
    }else if(AA){
      out$E <- out$E[PFAMorder, ]
      out$qe <- out$qe[PFAMorder]
    }
    return(out)
  }else stop("Invalid argument given for 'method'")
}
################################################################################
#' @rdname train
################################################################################
train.HMM <- function(x, y, method = "Viterbi", seqweights = NULL,
                      maxiter = 100,
                      deltaLL = 1E-07,
                      logspace = "autodetect", quiet = FALSE, modelend = FALSE,
                      pseudocounts = "Laplace", ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  if(identical(pseudocounts, "background")) pseudocounts <- "Laplace"
  # no method for background pseudocounts yet
  #method <- toupper(method)
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
        vitj <- Viterbi(out, y[[j]], logspace = TRUE, ... = ...)
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
        if(!quiet) cat("Iteration", i, "\nPaths were identical after",
                       i, "iterations\n")
        return(out)
      }else{
        if(!quiet) cat("Iteration", i, "\n")
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
    } else if(!identical(pseudocounts, "none")) stop("Invalid 'pseudocounts' argument")
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
          forwj <- forward(out, yj, logspace = TRUE, ... = ...)
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
    if(!logspace){
      out$A <- exp(out$A)
      out$E <- exp(out$E)
    }
    warning("Failed to converge on a local maximum. Try increasing 'maxiter',
        decreasing 'deltaLL' or modifying start parameters")
    return(out)
  }else stop("Invalid argument given for 'method'")
}
################################################################################
