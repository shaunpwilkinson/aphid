#' Iterative model refinement.
#'
#' Update model parameters using a list of training sequences, with either
#'   the Viterbi training or Baum-Welch algorithm.
#'
#' @param x an object of class \code{"HMM"} or \code{"PHMM"} specifying the
#'   initial parameter values.
#' @param y a list of training sequences whose hidden states are unknown.
#'   Accepted modes are "character" and "raw" (for "DNAbin" and "AAbin"
#'   objects).
#' @param method a character string specifying the iterative model training
#'   method to use. Accepted methods are \code{"Viterbi"} (the default)
#'   and \code{"BaumWelch"}.
#' @param seqweights either NULL (all sequences are given weights
#'   of 1), a numeric vector the same length as \code{y} representing
#'   the sequence weights used to derive the model, or a character string giving
#'   the method to derive the weights from the sequences
#'   (see \code{\link{weight}}).
#' @param wfactor numeric. The factor to multiply the sequence weights by.
#'   Defaults to 1.
#' @param k integer representing the k-mer size to be used in tree-based
#'   sequence weighting (if applicable). Defaults to 5. Note that higher
#'   values of k may be slow to compute and use excessive memory due to
#'   the large numbers of calculations required.
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
#' @param limit the proportion of alignment rows that must be identical
#'   between subsequent iterations for the Viterbi training algorithm
#'   to terminate. Defaults to 1.
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
#'   method outlined in Durbin et al (1998) chapter 5.7), \code{"inherited"}
#'   (match and insert columns are inherited from the input alignment),
#'   and \code{"none"} (all columns are assigned match states in the model).
#'   Alternatively, insert columns can be
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
#'   column (Durbin et al 1998, chapter 5.7). Only applicable when
#'   \code{inserts = "map"}.
#' @param alignment logical indicating whether the alignment used to
#'   derive the final model (if applicable) should be included as an element of
#'   the returned PHMM object. Defaults to FALSE.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if x is not a list.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string "autodetect" is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores
#'   available. Note that in this case there
#'   may be a tradeoff in terms of speed depending on the number and size
#'   of sequences to be aligned, due to the extra time required to initialize
#'   the cluster.
#'   Only applicable if x is an object of class \code{"PHMM"}.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @param ... aditional arguments to be passed to \code{"Viterbi"} (if
#'   \code{method = "Viterbi"}) or \code{"forward"} (if
#'   \code{method = "BaumWelch"}).
#' @return an object of class \code{"HMM"} or \code{"PHMM"}, depending
#'   on the input model \code{x}.
#' @details
#'   This function optimizes the parameters of a hidden Markov model
#'   (object class: \code{"HMM"}) or profile hidden Markov model
#'   (object class: \code{"PHMM"}) using the methods described in
#'   Durbin et al (1998) chapters 3.3 and 6.5, respectively.
#'   For standard HMMs, the function assumes the state sequence is unknown
#'   (as opposed to the \code{\link{deriveHMM}} function, which is used
#'   when the state sequence is known).
#'   For profile HMMs, the input object is generally a list of non-aligned
#'   sequences rather than an alignment (for which the \code{\link{derivePHMM}}
#'   function may be more suitable).
#'
#'   This function offers a choice of two model training methods,
#'   Viterbi training (also known as the segmental
#'   K-means algorithm (Juang & Rabiner 1990)), and the Baum Welch algorithm,
#'   a special case of the expectation-maximization (EM) algorithm that
#'   iteratively finds the locally (but not necessarily globally) optimal
#'   parameters of a HMM or PHMM.
#'
#'   The Viterbi training method is generally much faster, particularly for
#'   profile HMMs and when the multi-threading option is used
#'   (see the \code{"cores"} argument). The comparison in accuracy will depend
#'   on the nature of the problem, but personal experience suggests that
#'   the methods are comparable for training profile HMMs for DNA and
#'   amino acid sequences.
#'
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
#'   ## Baum Welch training for standard HMMs:
#'   ## The dishonest casino example from Durbin et al (1998) chapter 3.2
#'   states <- c("Begin", "Fair", "Loaded")
#'   residues <- paste(1:6)
#'   ### Define the transition probability matrix
#'   A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9), nrow = 3)
#'   dimnames(A) <- list(from = states, to = states)
#'   ### Define the emission probability matrix
#'   E <- matrix(c(rep(1/6, 6), rep(1/10, 5), 1/2), nrow = 2, byrow = TRUE)
#'   dimnames(E) <- list(states = states[-1], residues = residues)
#'   ### Build and plot the HMM object
#'   x <- structure(list(A = A, E = E), class = "HMM")
#'   op <- par(no.readonly = TRUE)
#'   par(mfrow = c(2, 1))
#'   plot(x, main = "Dishonest casino HMM before training")
#'   data(casino)
#'   x <- train(x, list(casino), method = "BaumWelch", deltaLL = 0.001)
#'   plot(x, main = "Dishonest casino HMM after training")
#'   par(op)
#' @name train
################################################################################
train <- function(x, y, ...){
  UseMethod("train")
}
################################################################################
#' @rdname train
################################################################################
train.PHMM <- function(x, y, method = "Viterbi", seqweights = "Henikoff",
                       wfactor = 1, k = 5, logspace = "autodetect",
                       maxiter = 100, limit = 1, deltaLL = 1E-07,
                       pseudocounts = "background", gap = "-",
                       fixqa = FALSE, fixqe = FALSE, maxsize = NULL,
                       inserts = "map", threshold = 0.5, lambda = 0,
                       alignment = FALSE, cores = 1, quiet = FALSE, ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  DNA <- .isDNA(y)
  AA <- .isAA(y)
  DI <- !all(x$A["DI", ] == if(logspace) -Inf else 0)
  ID <- !all(x$A["ID", ] == if(logspace) -Inf else 0)
  gap <- if(DNA) as.raw(4) else if(AA) as.raw(45) else gap
  if(!is.list(y)) y <- if(is.null(dim(y))) list(y) else unalign(y, gap = gap)
  n <- length(y)
  if(is.null(seqweights)){
    seqweights <- rep(1, n)
  }else if(identical(seqweights, "Henikoff")){
    seqweights <- weight(y, k = k, gap = gap, method = "Henikoff")
  }else if(identical(seqweights, "Gerstein")){
    seqweights <- weight(y, k = k, gap = gap, method = "Gerstein")
  }else{
    stopifnot(
      length(seqweights) == n,
      !any(is.na(seqweights)),
      mode(seqweights) %in% c("numeric", "integer")
    )
  }
  seqweights <- seqweights * wfactor
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
      apply(t(sapply(y, .tabulateDNA, ambiguities = TRUE)) *
              seqweights, 2, sum) + 1
    }else if(AA){
      apply(t(sapply(y, .tabulateAA, ambiguities = TRUE)) *
              seqweights, 2, sum) + 1
    }else{
      apply(t(sapply(y, .tabulateCH, residues = residues)) *
              seqweights, 2, sum) + 1
    }
    x$qe <- log(allecs/sum(allecs))
  }
  if(!is.null(x$qa)){
    if(!logspace) x$qa <- log(x$qa)
  }else{
    alig <- align(y, x, ... = ...)
    gaps <- alig == gap
    insrts <- apply(gaps, 2, sum) > 0.5 * n
    xtr <- matrix(nrow = n, ncol = ncol(alig))
    insertsn <- matrix(rep(insrts, n), nrow = n, byrow = T)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- .atab(xtr, seqweights = seqweights)
    transtotals <- apply(tcs, 1, sum) + 1 # force addition of Lapl pseudos
    if(!DI) transtotals[3] <- 0
    if(!ID) transtotals[7] <- 0
    x$qa <- log(transtotals/sum(transtotals))
  }
  ## set up multithread
  if(inherits(cores, "cluster")){
    para <- TRUE
    stopclustr <- FALSE
  }else if(cores == 1){
    para <- FALSE
    stopclustr <- FALSE
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(cores > 1){
      # if(cores > navailcores) stop("No. cores is more than number available")
      if(!quiet) cat("Multithreading over", cores, "cores\n")
      cores <- parallel::makeCluster(cores)
      para <- TRUE
      stopclustr <- TRUE
    }else{
      para <- FALSE
      stopclustr <- FALSE
    }
  }
  if(method  == "Viterbi"){
    alig <- align.list(y, model = x, logspace = TRUE, cores = cores, ... = ...)
    alig_cache <- character(maxiter + 1)
    hashis <- apply(alig, 1, .digest)
    hashi <- .digest(paste0(hashis, collapse = ""))
    stopifnot(length(hashi) == 1)
    alig_cache[1] <- hashi
    for(i in 1:maxiter){
      model <- derivePHMM.default(alig, seqweights = seqweights,
                                residues = residues, gap = gap,
                                DI = DI, ID = ID, maxsize = maxsize,
                                inserts = inserts, lambda = lambda,
                                threshold = threshold,
                                pseudocounts = pseudocounts,
                                logspace = TRUE, alignment = alignment,
                                qa = if(fixqa) x$qa else NULL,
                                qe = if(fixqe) x$qe else NULL)
      if(!quiet){
        cat("Iteration", i)
        cat(": alignment with", nrow(alig), "rows &", ncol(alig), "columns, ")
        cat("PHMM with", model$size, "modules\n")
      }
      rm(alig) ## free up space for next alignment
      gc()
      alig <- align(y, model = model, logspace = TRUE, cores = cores, ... = ...)
      tmp <- apply(alig, 1, .digest)
      breakme <- sum(tmp == hashis)/length(hashis) > limit
      hashis <- tmp
      hashi <- .digest(paste0(hashis, collapse = ""))
      stopifnot(length(hashi) == 1)
      #newhash <- .digest(alig)
      if(!hashi %in% alig_cache & !breakme){
      #if(!newhash %in% alig_cache){
        #alig_cache[i + 1] <- newhash
        alig_cache[i + 1] <- hashi
      }else{
        if(!logspace){
          model$A <- exp(model$A)
          model$E <- exp(model$E)
          model$qa <- exp(model$qa)
          model$qe <- exp(model$qe)
        }
        if(!quiet) cat("Sequential alignments were identical after",
                       i, "iterations\n")
        if(para & stopclustr) parallel::stopCluster(cores)
        gc()
        return(model)
      }
    }
    if(!quiet) cat("Sequential alignments were not identical after",
                   i, "iterations\n")
    if(para & stopclustr) parallel::stopCluster(cores)
    model <- derivePHMM.default(alig, seqweights = seqweights,
                                residues = residues, gap = gap,
                                DI = DI, ID = ID, maxsize = maxsize,
                                inserts = inserts, lambda = lambda,
                                threshold = threshold,
                                pseudocounts = pseudocounts,
                                logspace = TRUE, alignment = alignment,
                                qa = if(fixqa) x$qa else NULL,
                                qe = if(fixqe) x$qe else NULL)
    return(model)
  }else if(method == "BaumWelch"){
    if(DNA){
      NUCorder <- sapply(toupper(rownames(x$E)), match,
                         c("A", "T", "G", "C"))
      x$E <- x$E[NUCorder, ]
      x$qe <- x$qe[NUCorder]
      if(!(identical(toupper(rownames(x$E)), c("A", "T", "G", "C")))){
        stop("Invalid model for DNA, residue alphabet does not correspond to
              nucleotide alphabet")
      }
      y <- .encodeDNA(y, arity = 4, probs = exp(x$qe),
                      random = FALSE, na.rm = TRUE)
    }else if(AA){
      PFAMorder <- sapply(toupper(rownames(x$E)), match,
                          LETTERS[-c(2, 10, 15, 21, 24, 26)])
      x$E <- x$E[PFAMorder, ]
      x$qe <- x$qe[PFAMorder]
      if(!(identical(toupper(rownames(x$E)), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
        stop("Invalid model for AA, residue alphabet does not correspond to
              20-letter amino acid alphabet")
      }
      y <- .encodeAA(y, arity = 20, probs = exp(x$qe), random = FALSE, na.rm = TRUE)
    }else{
      y <- lapply(y, function(s) match(s[s != gap], residues) - 1)
      if(any(is.na(unlist(y, use.names = FALSE)))) {
        stop("Residues in sequence(s) are missing from the model")
      }
    }
    # these just provide preformatted containers for the pseudocounts
    Apseudocounts <- x$A
    Epseudocounts <- x$E
    qepseudocounts <- x$qe
    # don't need x$qa since not updated during iteration
    if(mode(pseudocounts) %in% c("numeric", "integer") & length(pseudocounts) == 1L){
      Apseudocounts[] <- Epseudocounts[] <- qepseudocounts[] <- pseudocounts
    }else if(identical(pseudocounts, "background")){
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
    model <- x
    LL <- -1E12
    count <- function(y, model, DI = FALSE, ID = FALSE, ...){
      # y is sequences with weights attr
      ## output list including A, E, qe (counts)
      n <- length(y)
      l <- model$size
      seqweights <- attr(y, "weights")
      A <- model$A
      E <- model$E
      qe <- model$qe
      A[] <- E[] <- qe[] <- 0 # counts
      A0 <- A
      E0 <- E
      qe0 <- qe
      A0[] <- E0[] <- qe0[] <- -Inf # log probs
      logPx <- 0
      for(j in 1:n){
        Aj <- A0
        Ej <- E0
        qej <- qe0  ### note position specific
        seqweight <- seqweights[j]
        yj <- y[[j]]
        nj <- length(yj)
        if(nj == 0){ # can occasionally simulate zero length seqs
          A["DD", 2:l] <- A["DD", 2:l] + seqweight
          logPx <- logPx + sum(c(model$A["MD", 1], model$A["DD", 2:l],
                                 model$A["DM", l + 1]))
        }else{
          forwj <- forward.PHMM(model, yj, logspace = TRUE,
                                odds = FALSE, ... = ...)
          Rj <- forwj$array
          logPxj <- forwj$score
          logPx <- logPx + logPxj
          backj <- backward.PHMM(model, yj, logspace = TRUE, odds = FALSE, ... = ...)
          Bj <- backj$array
          yj <- yj + 1 # R indexing style
          ## EMISSIONS
          eM <- Rj[-1, -1, "M"] + Bj[-1, -1, "M"]
          eI <- Rj[, -1, "I"] + Bj[, -1, "I"]
          f <- factor(yj, levels = seq_len(nrow(E)))
          Ej <- qej <- indices <- split(seq_along(yj), f = f)
          indlens <- sapply(indices, length)
          for(i in seq_len(nrow(E))){
            if(indlens[i] == 0){
              Ej[[i]] <- rep(-Inf, l)
              qej[[i]] <- rep(-Inf, l + 1)
            }else if(indlens[i] == 1){
              Ej[[i]] <- eM[, indices[[i]]]
              qej[[i]] <- eI[, indices[[i]]] ## I still has state 0 column
            }else{
              Ej[[i]] <- apply(eM[, indices[[i]]], 1, logsum)
              qej[[i]] <- apply(eI[, indices[[i]]], 1, logsum)
            }
          }
          Ej <- exp(do.call("rbind", Ej) - logPxj)
          qej <- exp(do.call("rbind", qej) - logPxj)
          rownames(Ej) <- rownames(qej) <- rownames(E)
          E <- E + (Ej * seqweight)
          ####qej <- apply(qej, 1, sum) ############
          qe <- qe + (qej * seqweight)
          ## TRANSITIONS
          Ey <- t(model$E[yj, ])
          qey <- model$qe[yj]
          qey <- matrix(rep(qey, l + 1), nrow = l + 1, byrow = TRUE)
          colnames(qey) <- colnames(Ey)
          rownames(qey) <- c("0", rownames(Ey))
          Aj["DD", 1:l] <- apply(Rj[1:l, , "D"] + model$A["DD", 1:l] +
                                   Bj[-1, , "D"] , 1, logsum)
          Aj["MD", 1:l] <- apply(Rj[1:l, , "M"] + model$A["MD", 1:l] +
                                   Bj[-1, , "D"] , 1, logsum)
          if(ID) Aj["ID", 1:l] <- apply(Rj[1:l, , "I"] + model$A["ID", 1:l] +
                                          Bj[-1, , "D"] , 1, logsum)
          Aj["DM", 1:l] <- apply(Rj[1:l, 1:nj, "D"] + model$A["DM", 1:l] +
                                   Ey + Bj[-1, -1, "M"], 1, logsum)
          Aj["DM", l+1] <- Rj[l + 1, nj + 1, "D"] + model$A["DM", l + 1]
          Aj["MM", 1:l] <- apply(Rj[1:l, 1:nj, "M"] + model$A["MM", 1:l] +
                                   Ey + Bj[-1, -1, "M"], 1, logsum)
          Aj["MM", l+1] <- Rj[l + 1, nj + 1, "M"] + model$A["MM", l + 1]
          Aj["IM", 1:l] <- apply(Rj[1:l, 1:nj, "I"] + model$A["IM", 1:l] +
                                   Ey + Bj[-1, -1, "M"], 1, logsum)
          Aj["IM", l+1] <- Rj[l + 1, nj + 1, "I"] + model$A["IM", l + 1]
          if(DI) Aj["DI", 1:l] <- apply(Rj[ ,1:nj, "D"] + model$A["DI", ] +
                                          qey + Bj[, -1, "I"], 1, logsum)
          Aj["MI", ] <- apply(Rj[, 1:nj, "M"] + model$A["MI", ] + qey +
                                Bj[, -1, "I"], 1, logsum)
          Aj["II", ] <- apply(Rj[, 1:nj, "I"] + model$A["II", ] + qey +
                                Bj[, -1, "I"], 1, logsum)
          Aj <- exp(Aj - logPxj)
          A <- A + (Aj * seqweight)
        }
      }
      return(list(A = A, E = E, qe = qe, logPx = logPx))
    }
    do_multi <- para & length(y) > length(cores)
    if(do_multi){
      f <- rep(seq(1, length(cores)), each = ceiling(n/length(cores)))[1:n]
      y <- split(y, f)
      swl <- split(seqweights, f)
      for(i in seq_along(y)) attr(y[[i]], "weights") <- swl[[i]]
    }else{
      attr(y, "weights") <- seqweights
    }
    for(i in 1:maxiter){
      if(do_multi){
        counts <- parallel::parLapply(cores, y, count, model = model,
                                      DI = DI, ID = ID, ... = ...)
        A <- Reduce("+", lapply(counts, "[[", 1))
        E <- Reduce("+", lapply(counts, "[[", 2))
        qe <- Reduce("+", lapply(counts, "[[", 3))
        logPx <- Reduce("+", lapply(counts, "[[", 4))
      }else{
        counts <- count(y, model = model, DI = DI, ID = ID, ... = ...)
        A <- counts$A
        E <- counts$E
        qe <- counts$qe
        logPx <- counts$logPx
      }
      ## add pseudocounts
      qe <- apply(qe, 1, sum)
      A <- A + Apseudocounts
      E <- E + Epseudocounts
      qe <- qe + qepseudocounts
      E <- t(E)
      E <- log(E/apply(E, 1, sum))
      E <- t(E)
      qe <- log(qe/sum(qe))
      A <- t(A)
      M <- matrix(1:9, 3)
      for(m in 1:3) A[, M[, m]] <- log(A[, M[, m]]/apply(A[, M[, m]], 1, sum))
      A <- t(A)
      A[1:3, 1] <- -Inf # replace NaNs generated when dividing by 0
      model$A <- A
      model$E <- E
      if(!fixqe) model$qe <- qe
      # logPx <- sum(logPxs) # page 62 eq 3.17
      if(!quiet) cat("Iteration", i, "log likelihood =", logPx, "\n")
      if(abs(LL - logPx) < deltaLL){
        if(!logspace){
          model$A <- exp(model$A)
          model$E <- exp(model$E)
          model$qe <- exp(model$qe)
        }
        if(DNA){
          model$E <- model$E[NUCorder, ]
          model$qe <- model$qe[NUCorder]
        }else if(AA){
          model$E <- model$E[PFAMorder, ]
          model$qe <- model$qe[PFAMorder]
        }
        if(para & stopclustr) parallel::stopCluster(cores)
        if(!quiet) cat("Convergence threshold reached after", i, "EM iterations\n")
        return(model)
      }
      LL <- logPx
      gc()
    }
    if(!quiet) cat("Warning: failed to converge. Try increasing 'maxiter' or modifying start parameters")
    if(!logspace){
      model$A <- exp(model$A)
      model$E <- exp(model$E)
      model$qe <- exp(model$qe)
    }
    if(DNA){
      model$E <- model$E[NUCorder, ]
      model$qe <- model$qe[NUCorder]
    }else if(AA){
      model$E <- model$E[PFAMorder, ]
      model$qe <- model$qe[PFAMorder]
    }
    if(para & stopclustr) parallel::stopCluster(cores)
    gc()
    return(model)
  }else stop("Invalid argument given for 'method'")
}
################################################################################
#' @rdname train
################################################################################
train.HMM <- function(x, y, method = "Viterbi", seqweights = NULL, wfactor = 1,
                      maxiter = 100, deltaLL = 1E-07,
                      logspace = "autodetect", quiet = FALSE, modelend = FALSE,
                      pseudocounts = "Laplace", ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  if(identical(pseudocounts, "background")) pseudocounts <- "Laplace"
  # no method for background pseudocounts yet
  if(is.list(y)){
  }else if(is.vector(y, mode = "character")){
    y <- list(y)
  }else stop("Invalid y argument")
  n <- length(y)
  if(is.null(seqweights)) seqweights <- rep(1, n)
  stopifnot(sum(seqweights) == n)
  seqweights <- seqweights * wfactor
  states <- rownames(x$A)
  nstates <- length(states)
  residues <- colnames(x$E)
  nres <- length(residues)
  model <- x
  if(!logspace){
    model$E <- log(model$E)
    model$A <- log(model$A)
  }
  if(method == "Viterbi"){
    for(i in 1:maxiter){
      samename <- logical(n)
      for(j in 1:n){
        vitj <- Viterbi(model, y[[j]], logspace = TRUE, ... = ...)
        pathchar <- states[-1][vitj$path + 1]
        if(identical(pathchar, names(y[[j]]))) samename[j] <- TRUE
        # need to be named to feed into deriveHMM
        names(y[[j]]) <- pathchar
      }
      if(all(samename)){
        if(!logspace){
          model$A <- exp(model$A)
          model$E <- exp(model$E)
        }
        if(!quiet) cat("Iteration", i, "\nPaths were identical after",
                       i, "iterations\n")
        return(model)
      }else{
        if(!quiet) cat("Iteration", i, "\n")
        model <- deriveHMM(y, seqweights = seqweights, residues = residues,
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
    E <- model$E
    A <- model$A
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
          forwj <- forward(model, yj, logspace = TRUE, ... = ...)
          Rj <- forwj$array
          logPxj <- forwj$score
          tmplogPx[j] <- logPxj
          backj <- backward(model, yj, logspace = TRUE)
          Bj <- backj$array
          tmpAj <- tmpA
          tmpEj <- tmpE
          tmpAj[] <- tmpEj[] <- 0
          for(k in states[-1]){
            tmpAj[1, -1] <- exp(A[1, -1] + E[, yj[1]] + Bj[, 1] - logPxj)
            tmpAj[-1, 1] <- exp(Rj[, nj] + A[-1, 1] - logPxj)
            for(l in states[-1]){
              tmpAj[k, l] <- exp(logsum(Rj[k, -nj] + A[k, l] + E[l, yj[-1]] + Bj[l, -1]) - logPxj)
            }
            for(b in residues){
              cond <- yj == b
              if(any(cond)) tmpEj[k, b] <- exp(logsum(Rj[k, cond] + Bj[k, cond]) - logPxj)
            }
          }
          # correct for sequence weight
          tmpA <- tmpA + tmpAj * seqweights[j]
          tmpE <- tmpE + tmpEj * seqweights[j]
        }
      }
      A[] <- log(tmpA/apply(tmpA, 1, sum))
      E[] <- log(tmpE/apply(tmpE, 1, sum))
      model$A <- A
      model$E <- E
      logPx <- sum(tmplogPx) # page 62 eq 3.17
      if(!quiet) cat("Iteration", i, "log likelihood =", logPx, "\n")
      if(abs(LL - logPx) < deltaLL){
        if(!logspace){
          model$A <- exp(model$A)
          model$E <- exp(model$E)
        }
        if(!quiet) cat("Convergence threshold reached after", i, "EM iterations\n")
        return(model)
      }
      LL <- logPx
    }
    if(!logspace){
      model$A <- exp(model$A)
      model$E <- exp(model$E)
    }
    warning("Failed to converge on a local maximum. Try increasing 'maxiter',
        decreasing 'deltaLL' or modifying start parameters")
    return(model)
  }else stop("Invalid argument given for 'method'")
}
################################################################################
