#' Backward algorithm.
#'
#' Backward algorithm for calculating the full (log) probability or odds
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
#' @name backward
#'
backward <- function(x, y, qe = NULL, logspace = "autodetect",  odds = TRUE,
                     windowspace = "all",
                     type = "global", DI = FALSE, ID = FALSE, cpp = TRUE){
  UseMethod("backward")
}


#' @rdname backward
backward.PHMM <- function(x, y, qe = NULL, logspace = "autodetect",
                          type = "global", odds = TRUE,
                          windowspace = "all", DI = FALSE, ID = FALSE, cpp = TRUE){
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  pp <- inherits(y, "PHMM")
  if(pp) stop("PHMM vs PHMM back comparison is not supported")
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
    #y <- DNA2pentadecimal(y, na.rm = TRUE)
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
  B <- array(-Inf, dim = c(n, m, length(states)))
  dimnames(B) <- list(x = 0:(n - 1), y = 0:(m - 1), state = states)
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
    B[n, m, ] <- A[c("DM", "MM", "IM"), n]
    if(odds) E <- E - qe
    if(cpp){
      res <- backward_PHMM(y, A, E, qe, qey, type, windowspace, DI, ID, DNA = pd, AA = pa)
      B[, , 1] <- res$Dmatrix
      B[, , 2] <- res$Mmatrix
      B[, , 3] <- res$Imatrix
      res$array <- B
    }else{
      if(type == 0){
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
      }else{
        B[n, , 2] <- B[, m, 2] <- 0 ### needs checking, in viterbi, forward too
      }
      for(i in (n - 1):1){
        for(j in (m - 1):1){
          if(j - i >= windowspace[1] & j - i <= windowspace[2]){
            #sij <- if(pd) DNAprobC2(y[j], E[, i]) else E[y[j] + 1, i]
            if(pd){
              sij <- DNAprobC2(y[j], E[, i])
            }else if(pa){
              sij <- AAprobC2(y[j], E[, i])
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
      if(type == 0){
        score <- B[1, 1, "M"]
      }else{
        stop("semiglobal type not available for forward.PHMM yet")
      }
      res <- structure(list(score = score,
                            odds = odds,
                            array = B),
                       # Dmatrix = B[, , 1],
                       # Mmatrix = B[, , 2],
                       # Imatrix = B[, , 3]),
                       class = 'fullprob')
    }
  }
  return(res)
  }


#' @rdname backward
backward.HMM <- function (x, y, logspace = "autodetect", cpp = TRUE){
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
  H <- length(states)
  if(length(y) == 0) structure(list(score = A[1, 1], array = NULL), class = 'fullprob')
  if(cpp){
    res <- backward_HMM(y, A, E)
    rownames(res$array) <- states
  }else{
    B <- array(NA, dim = c(H, n))#, dimnames = list(state = states, rolls = 1:n))
    rownames(B) = states
    B[, n] <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H) #ak0
    tmp <- matrix(nrow = H, ncol = H)
    if(DNA){
      for (i in n:2){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- A[k + 1, l + 1] + DNAprobC2(y[i], E[l, ]) + B[l, i]
          }
        }
        B[, i - 1] <- apply(tmp, 1, logsum)
      }
    }else if(AA){
      for (i in n:2){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- A[k + 1, l + 1] + AAprobC2(y[i], E[l, ]) + B[l, i]
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
    #names(logprobs) <- states
    for(l in 1:H){
      p <- if(DNA) DNAprobC2(y[1], E[l, ]) else if(AA) AAprobC2(y[1], E[l, ]) else E[l, y[1]]
      logprobs[l] <- A[1, l + 1] + p + B[l, 1] #termination
    }
    score <- logsum(logprobs)
    res <- structure(list(score = score, array = B, odds = FALSE),
                     class = 'fullprob')
  }
  return(res)
}

