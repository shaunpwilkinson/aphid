#' The forward algorithm.
#'
#' This function calculates the full (log) probability or odds
#'   of a sequence given a hidden Markov model or profile HMM using the
#'   forward dynamic programming algorithm.
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
#'
#'   For a thorough explanation of the backward, forward and Viterbi
#'   algorithms, see Durbin et al (1998) chapters 3.2 (HMMs) and 5.4 (PHMMs).
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Wilbur WJ, Lipman DJ (1983) Rapid similarity searches of nucleic acid and
#'   protein data banks. \emph{Proc Natl Acad Sci USA}, \strong{80}, 726-730.
#'
#' @seealso \code{\link{backward}}, \code{\link{Viterbi}}.
#' @examples
#'   ## Forward algorithm for standard HMMs:
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
#'   ### Find full probability of the sequence given the model
#'   data(casino)
#'   forward(x, casino)
#'   ###
#'   ## Forward algorithm for profile HMMs:
#'   ## Small globin alignment data from Durbin et al (1998) Figure 5.3
#'   data(globins)
#'   ### Derive a profile HMM from the alignment
#'   globins.PHMM <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
#'   plot(globins.PHMM, main = "Profile hidden Markov model for globins")
#'   ### Simulate a random sequence from the model
#'   suppressWarnings(RNGversion("3.5.0"))
#'   set.seed(999)
#'   simulation <- generate(globins.PHMM, size = 20)
#'   simulation ## "F" "S" "A" "N" "N" "D" "W" "E"
#'   ### Calculate the full (log) probability of the sequence given the model
#'   x <- forward(globins.PHMM, simulation, odds = FALSE)
#'   x # -23.0586
#'   ### Show the dynammic programming array
#'   x$array
#' @name forward
################################################################################
forward <- function(x, y, ...){
  UseMethod("forward")
}
################################################################################
#' @rdname forward
################################################################################
forward.PHMM <- function(x, y, qe = NULL, logspace = "autodetect",
                         odds = TRUE, windowspace = "all",
                         DI = FALSE, ID = FALSE, cpp = TRUE, ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  pp <- inherits(y, "PHMM")
  if(pp) stop("PHMM vs PHMM forward comparison is not supported")
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
  R <- array(-Inf, dim = c(n, m, length(states)))
  dimnames(R) <- list(x = 0:(n - 1), y = 0:(m - 1), state = states)
  R[1, 1, 2] <- 0
  if(pp){
    ### placeholder
  }else{
    if(identical(windowspace, "WilburLipman")){
      xseq <- generate.PHMM(x, size = 10 * ncol(x$A), random = FALSE,
                            AA = pa, DNA = pd)
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
    if(n == 1){
      R[1, 1, "M"] <- 0
      if(m > 1) R[1, 2, "I"] <- A["MI", 1]
      if(m > 2) for(j in 3:m) R[1, j, "I"] <- R[1, j - 1, "I"] + A["II", 1] + qey[j - 1]
      score <- if(m > 1) R[1, m, "I"] + A["IM", 1] else 0
      res <- structure(list(score = score, array = R), class = "DPA")
      return(res)
    }
    if(m == 1){
      R[1, 1, "M"] <- 0
      if(n > 1) R[2, 1, "D"] <- A["MD", 1]
      if(n > 2) for(i in 3:n) R[i, 1, "D"] <- R[i - 1, 1, "D"] + A["DD", i - 1]
      score <- if(n > 1) R[n, 1, "D"] + A["DM", n] else 0
      res <- structure(list(score = score, array = R), class = "DPA")
      return(res)
    }
    if(odds) E <- E - qe
    if(cpp){
      res <- .forwardP(y, A, E, qe, qey, windowspace, DI, ID, DNA = pd, AA = pa)
      R[, , 1] <- res$Dmatrix
      R[, , 2] <- res$Mmatrix
      R[, , 3] <- res$Imatrix
      res$array <- R
      res$Dmatrix <- NULL
      res$Mmatrix <- NULL
      res$Imatrix <- NULL
    }else{
      R[2, 1, "D"] <- A["MD", 1]
      for(i in 3:n) R[i, 1, "D"] <- R[i - 1, 1, "D"] + A["DD", i - 1]
      R[1, 2, "I"] <- A["MI", 1] + qey[1]
      for(j in 3:m) R[1, j, "I"] <- R[1, j - 1, "I"] + A["II", 1] + qey[j - 1]
      for(i in 2:n){
        for(j in 2:m){
          if(j - i >= windowspace[1] & j - i <= windowspace[2]){
            if(pd){
              sij <- .probDNA(y[j - 1], E[, i - 1])
            }else if(pa){
              sij <- .probAA(y[j - 1], E[, i - 1])
            }else{
              sij <- E[y[j - 1] + 1, i - 1]
            }
            Dcdt <- c(R[i - 1, j, "D"] + A["DD", i - 1],
                      R[i - 1, j, "M"] + A["MD", i - 1],
                      if(ID) R[i - 1, j, "I"] + A["ID", i - 1] else -Inf)
            Mcdt <- c(R[i - 1, j - 1, "D"] + A["DM", i - 1],
                      R[i - 1, j - 1, "M"] + A["MM", i - 1],
                      R[i - 1, j - 1, "I"] + A["IM", i - 1])
            Icdt <- c(if(DI) R[i, j - 1, "D"] + A["DI", i] else -Inf,
                      R[i, j - 1, "M"] + A["MI", i],
                      R[i, j - 1, "I"] + A["II", i])
            R[i, j, "D"] <- logsum(Dcdt)
            R[i, j, "M"] <- logsum(Mcdt) + sij
            R[i, j, "I"] <- logsum(Icdt) + qey[j - 1]
          }
        }
      }
      LLcdt <- c(R[n, m, "M"] + A["MM", n],
                 R[n, m, "I"] + A["IM", n],
                 R[n, m, "D"] + A["DM", n])
      score <- logsum(LLcdt)
      res <- structure(list(score = score,
                            odds = odds,
                            array = R),
                       class = "DPA")
    }
  }
  return(res)
}
################################################################################
#' @rdname forward
################################################################################
forward.HMM <- function (x, y, logspace = "autodetect", cpp = TRUE, ...){
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
  residues <- colnames(E)
  H <- length(states)
  if(length(y) == 0) structure(list(score = A[1, 1], array = NULL), class = "DPA")
  if(cpp){
    res <- .forwardH(y, A, E, DNA, AA)
    rownames(res$array) <- states
  }else{
    R <- array(NA, dim = c(H, n))
    rownames(R) <- states
    tmp <- matrix(nrow = H, ncol = H)
    if(DNA){
      for(k in 1:H) R[k, 1] <- .probDNA(y[1], E[k, ]) + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- R[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H) R[k, i] <- logsum(tmp[, k]) + .probDNA(y[i], E[k, ])
      }
    }else if(AA){
      for(k in 1:H) R[k, 1] <- .probAA(y[1], E[k, ]) + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- R[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H) R[k, i] <- logsum(tmp[, k]) + .probAA(y[i], E[k, ])
      }
    }else{
      y <- y + 1
      for(k in 1:H) R[k, 1] <- E[k, y[1]] + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- R[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H) R[k, i] <- logsum(tmp[, k]) + E[k, y[i]]
      }
    }
    ak0 <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H)
    score <- logsum(R[, n] + ak0)
    res <- structure(list(score = score, array = R, odds = FALSE),
                     class = "DPA")
  }
  return(res)
}
################################################################################
