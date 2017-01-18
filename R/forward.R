#' Forward algorithm.
#'
#' Forward algorithm for calculating the full (log) probability or odds
#' of a sequence given a hidden Markov model or profile HMM.
#'
#' @param x an object of class \code{PHMM} or \code{HMM}.
#' @param y a character vector representing a single instance of a sequence
#' hypothetically emitted by the model.
#' @param type character string indicating whether insert and delete states
#' at the beginning and end of the path should count towards the final score
#' ('global'; default). Note that semiglobal and local models
#' are not currently supported in this version.
#' @inheritParams Viterbi
#' @return a list containing the score and forward dynamic programming arrays.
#' @name forward
#' @export
#'
forward <- function(x, y, qe = NULL, logspace = "autodetect", odds = TRUE,
                    windowspace = "all", type = "global", eta = 0.01,
                    DI = FALSE, ID = FALSE, cpp = TRUE){
  UseMethod("forward")
}


#' @rdname forward
#' @export
#'
forward.PHMM <- function(x, y, qe = NULL, logspace = "autodetect",
                         type = "global", odds = TRUE,
                         windowspace = "all", DI = FALSE, ID = FALSE, cpp = TRUE){
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  pp <- inherits(y, "PHMM")
  if(pp) stop("PHMM vs PHMM forward comparison is not supported")
  pd <- is.DNA(y)
  pa <- is.AA(y)
  pc <- !pp & !pd & !pa
  if(pd){
    rownames(x$E) <- toupper(rownames(x$E))
    if("U" %in% rownames(x$E)) rownames(x$E)[rownames(x$E) == "U"] <- "T"
    NUCorder <- sapply(rownames(x$E), match, c("A", "T", "G", "C"))
    x$E <- x$E[NUCorder, ]
    if(!(identical(rownames(x$E), c("A", "T", "G", "C")))){
      stop("invalid model for DNA, residue alphabet does not correspond to
           nucleotide alphabet")
    }
    if(is.list(y)){
      if(length(y) == 1){
        y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "DNAbin"
      }else stop("Invalid input object y: multi-sequence list")
    }
    y.DNAbin <- y
    y <- encode.DNA(y, arity = 15, na.rm = TRUE)
  }else if(pa){
    rownames(x$E) <- toupper(rownames(x$E))
    PFAMorder <- sapply(rownames(x$E), match, LETTERS[-c(2, 10, 15, 21, 24, 26)])
    x$E <- x$E[PFAMorder, ]
    if(!(identical(rownames(x$E), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
      stop("invalid model for AA, residue alphabet does not correspond to
           20-letter amino acid alphabet")
    }
    if(is.list(y)){
      if(length(y) == 1){
        y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "AAbin"
      }else stop("Invalid input object y: multi-sequence list")
    }
    y.AAbin <- y
    #y <- AA2heptovigesimal(y, na.rm = TRUE)
    y <- encode.AA(y, arity = 27, na.rm = TRUE)
  }else if(pc){
    if(is.list(y)){
      if(length(y) == 1){
        y <- y[[1]]
      }else stop("Invalid input object y: multi-sequence list")
    }
    #y <- setNames(seq_along(colnames(x$E)) - 1, colnames(x$E))[y]
    if(mode(y) == "character"){
      y <- match(y, rownames(x$E)) - 1
      if(any(is.na(y))) {
        warning("residues in sequence(s) are missing from the model")
        y <- y[!is.na(y)]
      }
    }#else if length(unique(y)) > nrow(x$E) stop("")
  }
  n <- ncol(x$E) + 1
  m <- if(pp) ncol(y$E) + 1 else length(y) + 1
  # if(identical(windowspace, "WilburLipman") | identical(windowspace, "all")){
  #   windowspace <- c(-x$size, if(pp) y$size else length(y)) ### placeholder
  # }else if(length(windowspace) != 2) stop("invalid windowspace argument")
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
  type = switch(type, "global" = 0L, "semiglobal" = 1L, stop("invalid type"))
  R <- array(-Inf, dim = c(n, m, length(states)))
  dimnames(R) <- list(x = 0:(n - 1), y = 0:(m - 1), state = states)
  R[1, 1, 2] <- 0
  if(pp){
    ### placeholder
  }else{
    if(identical(windowspace, "WilburLipman")){
      xseq <- generate.PHMM(x, size = 10 * ncol(x$A), random = FALSE, AA = pa, DNA = pd)
      if(pd){
        xqt  <- match(xseq, as.raw(c(136, 24, 72, 40))) - 1
        #yqt <- DNA2quaternary(y.DNAbin, na.rm = TRUE)
        yqt <- encode.DNA(y.DNAbin, arity = 4, na.rm = TRUE)
        windowspace <- WilburLipman(xqt, yqt, arity = 4, k = 5)
      }else if(pa){
        y.comp <- encode.AA(y.AAbin, arity = 6, na.rm = TRUE)
        xseq.comp <- encode.AA(xseq, arity = 6, na.rm = TRUE)
        windowspace <- WilburLipman(xseq.comp, y.comp, arity = 6, k = 5)
      }else{
        xseq <- match(xseq, rownames(x$E)) - 1
        windowspace <- WilburLipman(xseq, y, arity = nrow(x$E), k = 3)
      }
    }else if(identical(windowspace, "all")){
      windowspace <- c(-x$size, length(y))
    }else if(length(windowspace) != 2) stop("invalid windowspace argument")
    qey <- if(odds) {
      rep(0, m - 1)
    }else if(pd){
      sapply(y, DNAprobC2, qe)
    }else if(pa){
      sapply(y, AAprobC2, qe)
    }else qe[y + 1]
    A <- if(logspace) x$A else log(x$A)
    E <- if(logspace) x$E else log(x$E)
    if(n == 1){
      R[1, 1, "M"] <- 0
      if(m > 1) R[1, 2, "I"] <- A["MI", 1]
      if(m > 2) for(j in 3:m) R[1, j, "I"] <- R[1, j - 1, "I"] + A["II", 1] + qey[j - 1]
      score <- if(m > 1) R[1, m, "I"] + A["IM", 1] else 0
      res <- structure(list(score = score, array = R), class = 'fullprob')
      return(res)
    }
    if(m == 1){
      R[1, 1, "M"] <- 0
      if(n > 1) R[2, 1, "D"] <- A["MD", 1]
      if(n > 2) for(i in 3:n) R[i, 1, "D"] <- R[i - 1, 1, "D"] + A["DD", i - 1]
      score <- if(n > 1) R[n, 1, "D"] + A["DM", n] else 0
      res <- structure(list(score = score, array = R), class = 'fullprob')
      return(res)
    }
    if(odds) E <- E - qe
    if(cpp){
      res <- forward_PHMM(y, A, E, qe, qey, type, windowspace, DI, ID, DNA = pd, AA = pa)
      R[, , 1] <- res$Dmatrix
      R[, , 2] <- res$Mmatrix
      R[, , 3] <- res$Imatrix
      res$array <- R
    }else{
      if(type == 0){
        R[2, 1, "D"] <- A["MD", 1]
        for(i in 3:n) R[i, 1, "D"] <- R[i - 1, 1, "D"] + A["DD", i - 1]
        R[1, 2, "I"] <- A["MI", 1] + qey[1]
        for(j in 3:m) R[1, j, "I"] <- R[1, j - 1, "I"] + A["II", 1] + qey[j - 1]
      }else{
        R[-1, 1, 1] <- R[1, -1, 3] <- 0 ### needs checking, in viterbi too
      }
      for(i in 2:n){
        for(j in 2:m){
          if(j - i >= windowspace[1] & j - i <= windowspace[2]){
            if(pd){
              sij <- DNAprobC2(y[j - 1], E[, i - 1])
            }else if(pa){
              sij <- AAprobC2(y[j - 1], E[, i - 1])
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
      if(type == 0){
        LLcdt <- c(R[n, m, "M"] + A["MM", n],
                   R[n, m, "I"] + A["IM", n],
                   R[n, m, "D"] + A["DM", n])
        score <- logsum(LLcdt)
      }else{
        stop("semiglobal type not available for forward.PHMM yet")
      }
      res <- structure(list(score = score,
                            odds = odds,
                            array = R),
                            # Dmatrix = R[, , 1],
                            # Mmatrix = R[, , 2],
                            # Imatrix = R[, , 3]),
                       class = 'fullprob')
    }
  }
  return(res)
}


#' @rdname forward
#' @export
#'
forward.HMM <- function (x, y, logspace = "autodetect", cpp = TRUE){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  DNA <- is.DNA(y)
  AA <- is.AA(y)
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
    #y <- DNA2pentadecimal(y, na.rm = TRUE)
    y <- encode.DNA(y, arity = 15, na.rm = TRUE)
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
    #y <- AA2heptovigesimal(y, na.rm = TRUE)
    y <- encode.AA(y, arity = 27, na.rm = TRUE)
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
  if(length(y) == 0) structure(list(score = A[1, 1], array = NULL), class = 'fullprob')
  if(cpp){
    res <- forward_HMM(y, A, E, DNA, AA)
    rownames(res$array) <- states
  }else{
    R <- array(NA, dim = c(H, n))#, dimnames = list(state = states, step = 1:n))
    rownames(R) <- states
    tmp <- matrix(nrow = H, ncol = H)
    if(DNA){
      for(k in 1:H) R[k, 1] <- DNAprobC2(y[1], E[k, ]) + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- R[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H) R[k, i] <- logsum(tmp[, k]) + DNAprobC2(y[i], E[k, ])
      }
    }else if(AA){
      for(k in 1:H) R[k, 1] <- AAprobC2(y[1], E[k, ]) + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- R[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H) R[k, i] <- logsum(tmp[, k]) + AAprobC2(y[i], E[k, ])
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
                     class = 'fullprob')
  }
  return(res)
}

#
# forward.default <- function(x, y, type = 'global', d = 8, e = 2,
#                            S = NULL, windowspace = "all"){
#   ###check that x is a character vector###
#   # this is invalid, d and e have no probabilistic interpretation
#   if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
#   n <- length(x) + 1
#   m <- length(y) + 1
#   if (is.null(S)) {
#     residues <- unique(c(x, y))
#     S <- diag(2, nrow = length(residues))
#     S <- S - 1
#     dimnames(S) <- list(residues, residues)
#   }
#   # initialize scoring and pointer arrays (M and P)
#   M <- array(-Inf, dim = c(n, m, 3))
#   M[1, 1, 2] <- 0
#   if(type == 'global'){
#     M[, 1, 1] <- c(0, seq(from = -d, to = (- d + (n - 2) * -e), by = -e))
#     M[1, , 3] <- c(0, seq(from = -d, to = (- d + (m - 2) * -e), by = -e))
#   }else{
#     M[2:n, 1, 1] <- M[1, 2:m, 3] <- 0 ### check this - should fill dim2 instead?
#   }
#   # recursion step
#   for(i in if(is.null(itertab)) 2:n else 1:nrow(itertab)){
#     for(j in if(is.null(itertab)) 2:m else 1){
#       if(!is.null(itertab)){
#         j <- itertab[i, 2]
#         i <- itertab[i, 1]
#       }
#       sij <-  S[x[i - 1], y[j - 1]]
#       M[i, j, 1] <- logsum(c(M[i - 1, j, 1] - e, M[i - 1, j, 2] - (d + e)))#x alig to gap in y
#       M[i, j, 2] <- logsum(c(M[i - 1, j - 1, 1],
#                              M[i - 1, j - 1, 2],
#                              M[i - 1, j - 1, 3])) + sij
#       if(type == 'local' & M[i, j, 2] < 0) M[i, j, 2] <- 0
#       M[i, j, 3] <- logsum(c(-Inf, M[i, j - 1, 2] - (d + e), M[i, j - 1, 3] - e))
#     }
#   }
#   if(type == 'global'){
#     # find highest score in bottom right corner of scoring array M
#     score <- logsum(M[n, m, ])
#   }else if(type == 'semiglobal'){
#     # find highest score on bottom row or right column of scoring array M
#     border <- rbind(M[n, , ], M[-n, m, ])
#     ind <- which(border == max(border), arr.ind = TRUE)
#     if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
#     z <- if(ind[1] <= m) c(n, ind) else c(ind[1] - m, m, ind[2])
#     score <- logsum(M[z[1], z[2], ])
#   }else if(type == 'local'){
#     # find highest score in scoring array M
#     ind <- which(M[, , 2] == max(M[, , 2]), arr.ind = TRUE)
#     score <- logsum(M[ind[1], ind[2], ])
#   }
#   res <- structure(list(score = score, array = M), class = 'forward')
#   return(res)
# }
#
