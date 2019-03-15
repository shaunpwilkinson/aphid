#' The backward algorithm.
#'
#' This function calculates the full (log) probability or odds
#'   of a sequence given a hidden Markov model or profile HMM using the
#'   backward dynamic programming algorithm.
#'
#' @param x an object of class \code{PHMM} or \code{HMM}.
#' @param y a vector of mode "character" or "raw" (a "DNAbin" or "AAbin"
#'   object) representing a single sequence hypothetically emitted by
#'   the model in \code{x}.
#' @inheritParams Viterbi
#' @return an object of class \code{"DPA"}, which is a list
#'   containing the score and dynamic programming array.
#' @details
#'   This function is a wrapper for a compiled C++ function that recursively
#'   fills a dynamic programming matrix with logged probabilities, and
#'   calculates the full (logged) probability of a sequence given a HMM or
#'   PHMM.
#'   For a thorough explanation of the backward, forward and Viterbi
#'   algorithms, see Durbin et al (1998) chapters 3.2 (HMMs) and 5.4 (PHMMs).
#'
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Wilbur WJ, Lipman DJ (1983) Rapid similarity searches of nucleic acid and
#'   protein data banks. \emph{Proc Natl Acad Sci USA}, \strong{80}, 726-730.
#'
#' @seealso \code{\link{forward}}, \code{\link{Viterbi}}.
#' @examples
#'   ## Backward algorithm for standard HMMs:
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
#'   plot(x, main = "Dishonest casino HMM")
#'   data(casino)
#'   backward(x, casino)
#'   ##
#'   ## Backward algorithm for profile HMMs:
#'   ## Small globin alignment data from Durbin et al (1998) Figure 5.3
#'   data(globins)
#'   ### Derive a profile hidden Markov model from the alignment
#'   globins.PHMM <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
#'   plot(globins.PHMM, main = "Profile hidden Markov model for globins")
#'   ### Simulate a random sequence from the model
#'   suppressWarnings(RNGversion("3.5.0"))
#'   set.seed(999)
#'   simulation <- generate(globins.PHMM, size = 20)
#'   simulation ## "F" "S" "A" "N" "N" "D" "W" "E"
#'   ### Calculate the full (log) probability of the sequence given the model
#'   x <- backward(globins.PHMM, simulation, odds = FALSE)
#'   x # -23.0586
#'   ### Show dynammic programming array
#'   x$array
#' @name backward
################################################################################
backward <- function(x, y, ...){
  UseMethod("backward")
}
################################################################################
#' @rdname backward
################################################################################
backward.PHMM <- function(x, y, qe = NULL, logspace = "autodetect",
                          odds = TRUE, windowspace = "all",
                          DI = FALSE, ID = FALSE, cpp = TRUE, ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  pp <- inherits(y, "PHMM")
  if(pp) stop("PHMM vs PHMM backward comparison is not supported")
  pd <- .isDNA(y)
  pa <- .isAA(y)
  pc <- !pp & !pd & !pa
  if(!pp & is.list(y)) if(length(y) == 1) y <- y[[1]] else stop("y is invalid")
  if(pd){
    rownames(x$E) <- toupper(rownames(x$E))
    if("U" %in% rownames(x$E)) rownames(x$E)[rownames(x$E) == "U"] <- "T"
    NUCorder <- sapply(rownames(x$E), match, c("A", "T", "G", "C"))
    x$E <- x$E[NUCorder, ]
    if(!(identical(rownames(x$E), c("A", "T", "G", "C")))){
      stop("Invalid model for DNA, residue alphabet does not correspond to
           nucleotide alphabet")
    }
    class(y) <- "DNAbin"
    y.DNAbin <- y
    y <- .encodeDNA(y, arity = 15, na.rm = TRUE)
  }else if(pa){
    rownames(x$E) <- toupper(rownames(x$E))
    PFAMorder <- sapply(rownames(x$E), match, LETTERS[-c(2, 10, 15, 21, 24, 26)])
    x$E <- x$E[PFAMorder, ]
    if(!(identical(rownames(x$E), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
      stop("Invalid model for AA, residue alphabet does not correspond to
             20-letter amino acid alphabet")
    }
    class(y) <- "AAbin"
    y.AAbin <- y
    #y <- AA2heptovigesimal(y, na.rm = TRUE)
    y <- .encodeAA(y, arity = 27, na.rm = TRUE)
  }else if(pc){
    if(mode(y) == "character"){
      y <- match(y, rownames(x$E)) - 1
      if(any(is.na(y))){
        warning("Residues in sequence(s) are missing from the model")
        y <- y[!is.na(y)]
      }
    }
  }
  n <- ncol(x$E) + 1
  m <- if(pp) ncol(y$E) + 1 else length(y) + 1
  states <- if(pp) c("MI", "DG", "MM", "GD", "IM") else c("D", "M", "I")
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
  B <- array(-Inf, dim = c(n, m, length(states)))
  dimnames(B) <- list(x = 0:(n - 1), y = 0:(m - 1), state = states)
  if(pp){
    ### placeholder
  }else{
    if(identical(windowspace, "WilburLipman")){
      xseq <- generate.PHMM(x, size = 10 * ncol(x$A), random = FALSE, AA = pa, DNA = pd)
      if(pd){
        xqt  <- match(xseq, as.raw(c(136, 24, 72, 40))) - 1
        yqt <- .encodeDNA(y.DNAbin, arity = 4, na.rm = TRUE)
        windowspace <- .streak(xqt, yqt, arity = 4, k = 5)
      }else if(pa){
        y.comp <- .encodeAA(y.AAbin, arity = 6, na.rm = TRUE)
        xseq.comp <- .encodeAA(xseq, arity = 6, na.rm = TRUE)
        windowspace <- .streak(xseq.comp, y.comp, arity = 6, k = 5)
      }else{
        xseq <- match(xseq, rownames(x$E)) - 1
        windowspace <- .streak(xseq, y, arity = nrow(x$E), k = 3)
      }
    }else if(identical(windowspace, "all")){
      windowspace <- c(-x$size, length(y))
    }else if(length(windowspace) != 2) stop("Invalid windowspace argument")
    qey <- if(odds) {
      rep(0, m - 1)
    }else if(pd){
      sapply(y, .probDNA, qe)
    }else if(pa){
      sapply(y, .probAA, qe)
    }else qe[y + 1]
    A <- if(logspace) x$A else log(x$A)
    E <- if(logspace) x$E else log(x$E)
    B[n, m, ] <- A[c("DM", "MM", "IM"), n]
    if(odds) E <- E - qe
    if(cpp){
      res <- .backwardP(y, A, E, qe, qey, windowspace, DI, ID, DNA = pd, AA = pa)
      B[, , 1] <- res$Dmatrix
      B[, , 2] <- res$Mmatrix
      B[, , 3] <- res$Imatrix
      res$array <- B
      res$Dmatrix <- NULL
      res$Mmatrix <- NULL
      res$Imatrix <- NULL
    }else{
      for(i in n:2) {
        B[i - 1, m, "D"] <- B[i, m, "D"] + A["DD", i - 1]
        B[i - 1, m, "M"] <- B[i, m, "D"] + A["MD", i - 1]
        B[i - 1, m, "I"] <- if(ID) B[i, m, "D"] + A["ID", i - 1] else -Inf
      }
      for(j in m:2) {
        B[n, j - 1, "D"] <- if(DI) B[n, j, "I"] + A["DI", n] + qey[j - 1] else -Inf
        B[n, j - 1, "M"] <- B[n, j, "I"] + A["MI", n] + qey[j - 1]
        B[n, j - 1, "I"] <- B[n, j, "I"] + A["II", n] + qey[j - 1]
      }
      for(i in (n - 1):1){
        for(j in (m - 1):1){
          if(j - i >= windowspace[1] & j - i <= windowspace[2]){
            if(pd){
              sij <- .probDNA(y[j], E[, i])
            }else if(pa){
              sij <- .probAA(y[j], E[, i])
            }else{
              sij <- E[y[j] + 1, i]
            }
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
      }
      score <- B[1, 1, "M"]
      res <- structure(list(score = score,
                            odds = odds,
                            array = B),
                       class = "DPA")
    }
  }
  return(res)
  }
################################################################################
#' @rdname backward
################################################################################
backward.HMM <- function (x, y, logspace = "autodetect", cpp = TRUE, ...){
  if(identical(logspace, 'autodetect')) logspace <- .logdetect(x)
  DNA <- .isDNA(y)
  AA <- .isAA(y)
  if(DNA){
    colnames(x$E) <- toupper(colnames(x$E))
    NUCorder <- sapply(colnames(x$E), match, c("A", "T", "G", "C"))
    x$E <- x$E[, NUCorder]
    if(!(identical(colnames(x$E), c("A", "T", "G", "C")))){
      stop("invalid model for DNA, residue alphabet does not correspond to
           nucleotide alphabet")
    }
    if(is.list(y)){
      if(length(y) == 1){
        y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "DNAbin"
      }else stop("Invalid input object y: multi-sequence list")
    }
    y <- .encodeDNA(y, arity = 15, na.rm = TRUE)
  }else if(AA){
    colnames(x$E) <- toupper(colnames(x$E))
    PFAMorder <- sapply(colnames(x$E), match, LETTERS[-c(2, 10, 15, 21, 24, 26)])
    x$E <- x$E[, PFAMorder]
    if(!(identical(colnames(x$E), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
      stop("invalid model for AA, residue alphabet does not correspond to
           20-letter amino acid alphabet")
    }
    if(is.list(y)){
      if(length(y) == 1){
        y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "AAbin"
      }else stop("Invalid input object y: multi-sequence list")
    }
    y <- .encodeAA(y, arity = 27, na.rm = TRUE)
  }else{
    if(is.list(y)){
      if(length(y) == 1){
        y <- y[[1]]
      }else stop("Invalid input object y: multi-sequence list")
    }
    if(mode(y) == "character") y <- match(y, colnames(x$E)) - 1
  }
  n <- length(y)
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  states <- rownames(E)
  H <- length(states)
  if(length(y) == 0) structure(list(score = A[1, 1], array = NULL), class = "DPA")
  if(cpp){
    res <- .backwardH(y, A, E)
    rownames(res$array) <- states
  }else{
    B <- array(NA, dim = c(H, n))
    rownames(B) = states
    B[, n] <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H) #ak0
    tmp <- matrix(nrow = H, ncol = H)
    if(DNA){
      for (i in n:2){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- A[k + 1, l + 1] + .probDNA(y[i], E[l, ]) + B[l, i]
          }
        }
        B[, i - 1] <- apply(tmp, 1, logsum)
      }
    }else if(AA){
      for (i in n:2){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- A[k + 1, l + 1] + .probAA(y[i], E[l, ]) + B[l, i]
          }
        }
        B[, i - 1] <- apply(tmp, 1, logsum)
      }
    }else{
      y <- y + 1
      for (i in n:2){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- A[k + 1, l + 1] + E[l, y[i]] + B[l, i]
          }
        }
        B[, i - 1] <- apply(tmp, 1, logsum)
      }
    }
    logprobs <- numeric(H)
    for(l in 1:H){
      p <- if(DNA) .probDNA(y[1], E[l, ]) else if(AA) .probAA(y[1], E[l, ]) else E[l, y[1]]
      logprobs[l] <- A[1, l + 1] + p + B[l, 1] #termination
    }
    score <- logsum(logprobs)
    res <- structure(list(score = score, array = B, odds = FALSE),
                     class = "DPA")
  }
  return(res)
}
################################################################################
