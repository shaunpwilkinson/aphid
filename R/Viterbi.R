#' Find the optimal sequence path recursively
#'
#' \code{Viterbi} finds the optimal path of a sequence through a HMM
#' or PHMM and returns its log-odds score.
#'
#' @param x an object of class \code{HMM} or \code{PHMM}, or a character vector.
#' @param y a character vector consisting of residues emitted by the
#' HMM or PHMM.
#' @param qe an optional named vector of background residue frequencies. If NULL
#' background residue frequencies from the PHMMs are used. If these are not available
#' equal background residue frequencies are assumed.
#' @param logspace logical argument indicating whether the emission and transition
#' probabilities of x are logged (base e; TRUE) or raw (FALSE). Alternatively, if
#' \code{logspace = "autodetect"} (default), the function will automatically detect
#' if the probabilities are in log space, returning an error if inconsistencies are found.
#' Note that choosing the latter option increases the computational
#' overhead; therefore specifying \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param type character string indicating whether insert and delete states
#' at the beginning and end of the path should count towards the final score
#' ('global'; default), or not ('semiglobal'), or whether the highest scoring
#' sub-path should be returned ('local').
#' @param S an optional scoring matrix with rows and columns named according
#' to the residue alphabet. Note that for local alignments scores for
#' mismatches should generally take negative values or else spurious
#' alignments could occur. If NULL matches are scored as 1 and
#' mismatches scored as -1.
#' @param offset column score offset to specify level of greediness. Defaults to
#' -0.1 bits as recommended by Soding (2005).
#' @param itertab an optional two column matrix of row and column indices used
#' to iterate through a subset of the viterbi array.
#' @param d gap opening penalty for sequence vs. sequence alignment
#' @param e gap extension penalty for sequence vs. sequence alignment
#' @name Viterbi
#'
#' @examples
#' x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#' y <- c("P", "A", "W", "H", "E", "A", "E")
#' Viterbi(x, y,  d = 8, e = 2)
#'
Viterbi <- function(x, y, qe = NULL, logspace = "autodetect", type = "semiglobal",
                    offset = -0.1, d = 8, e = 2, S = NULL, itertab = NULL){
  UseMethod("Viterbi")
}

#' @rdname Viterbi
Viterbi.PHMM <- function(x, y, qe = NULL, logspace = "autodetect",
                         type = "semiglobal", offset = -0.1,
                         itertab = NULL){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
  if(identical(itertab, "WilburLipman")) itertab <- NULL ### placeholder
  pp <- inherits(y, "PHMM")
  #pd <- mode(y) == "raw"
  pd <- inherits(y, "DNAbin")
  if(pd){
    x$E <- x$E[order(rownames(x$E)),]
    if(!(identical(rownames(x$E), c("a", "c", "g", "t")) |
         identical(rownames(x$E), c("A", "C", "G", "T")))){
      stop("invalid model, residue alphabet does not correspond to
           nucleotide alphabet")
    }
    if(is.list(y)) y <- y[[1]]
    if(is.matrix(y)){
      if(nrow(y) == 1) {
        y <- as.vector(y)
      } else {
        stop("Viterbi can only process one sequence at a time")
      }
    }
  }
  n <- ncol(x$E)
  m <- if(pp) ncol(y$E) else length(y)
  states <- if(pp) c("MI", "DG", "MM", "GD", "IM") else c("D", "M", "I")
  V <- array(-Inf, dim = c(n + 1, m + 1, length(states)),
             dimnames = list(x = 0:n, y = 0:m, state = states))
  P <- V + NA # pointer array
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
  if(pp){
    if(logdetect(y) != logspace) stop("Both models must be in the same logspace format")
    Ax <- if(logspace) x$A else log(x$A)
    Ay <- if(logspace) y$A else log(y$A)
    Ex <- if(logspace) x$E else log(x$E)
    Ey <- if(logspace) y$E else log(y$E)
    fun <- Vectorize(function(a, b) logsum((Ex[, a] + Ey[, b]) - qe))
    Saa <- outer(1:n, 1:m, fun)
    if(type == "global"){
      V[-1, 1, "MI"] <- cumsum(c(Ax["MM", 1] + Ay["MI", 1], Ax["MM", 2:n] + Ay["II", 1]))
      P[-1, 1, "MI"] <- c(3, rep(1, n - 1))
      V[-1, 1, "DG"] <- cumsum(c(Ax["MD", 1], Ax["DD", 2:n]))
      P[-1, 1, "DG"] <- c(3, rep(2, n - 1))
      V[1, 1, "MM"] <- 0
      V[1, -1, "GD"] <- cumsum(c(Ay["MD", 1], Ay["DD", 2:m]))
      P[1, -1, "GD"] <- c(3, rep(4, m - 1))
      V[1, -1, "IM"] <- cumsum(c(Ax["MI", 1] + Ay["MM", 1], Ax["II", 1] + Ay["MM", 2:m]))
      P[1, -1, "IM"] <- c(3, rep(5, m - 1))
    }else{
      V[1, , "MM"] <- V[, 1, "MM"] <- 0
    }
    for(i in if(is.null(itertab)) 1:n else 1:nrow(itertab)){
      for(j in if(is.null(itertab)) 1:m else 1){
        if(!is.null(itertab)){
          j <- itertab[i, 2]
          i <- itertab[i, 1]
        }
        sij <- Saa[i, j] + offset
        MIcdt <- c(V[i, j + 1, "MM"] + Ax["MM", i] + Ay["MI", j + 1],
                   V[i, j + 1, "MI"] + Ax["MM", i] + Ay["II", j + 1])
        DGcdt <- c(V[i, j + 1, "MM"] + Ax["MD", i],
                   V[i, j + 1, "DG"] + Ax["DD", i])
        MMcdt <- c(V[i, j, "MI"] + Ax["MM", i] + Ay["IM", j] + sij,
                   V[i, j, "DG"] + Ax["DM", i] + Ay["MM", j] + sij,
                   V[i, j, "MM"] + Ax["MM", i] + Ay["MM", j] + sij,
                   V[i, j, "GD"] + Ax["MM", i] + Ay["DM", j] + sij,
                   V[i, j, "IM"] + Ax["IM", i] + Ay["MM", j] + sij,
                   if(type == "local") 0 else NULL)
        GDcdt <- c(V[i + 1, j, "MM"] + Ay["MD", j],
                   V[i + 1, j, "GD"] + Ay["DD", j])
        IMcdt <- c(V[i + 1, j, "MM"] + Ax["MI", i + 1] + Ay["MM", j],
                   V[i + 1, j, "IM"] + Ax["II", i + 1] + Ay["MM", j])
        MImax <- whichismax(MIcdt)
        DGmax <- whichismax(DGcdt)
        MMmax <- whichismax(MMcdt)
        GDmax <- whichismax(GDcdt)
        IMmax <- whichismax(IMcdt)
        V[i + 1, j + 1, "MI"] <- MIcdt[MImax]
        V[i + 1, j + 1, "DG"] <- DGcdt[DGmax]
        V[i + 1, j + 1, "MM"] <- MMcdt[MMmax]
        V[i + 1, j + 1, "GD"] <- GDcdt[GDmax]
        V[i + 1, j + 1, "IM"] <- IMcdt[IMmax]
        P[i + 1, j + 1, "MI"] <- c(3, 1)[MImax]
        P[i + 1, j + 1, "DG"] <- c(3, 2)[DGmax]
        P[i + 1, j + 1, "MM"] <- c(1:6)[MMmax]
        P[i + 1, j + 1, "GD"] <- c(3, 4)[GDmax]
        P[i + 1, j + 1, "IM"] <- c(3, 5)[IMmax]
      }
    }
    path <- c()
    progression <- matrix(nrow = 2, ncol = 0)
    if(type == 'global'){
      LLcdt <- c(V[n + 1, m + 1, "MI"] + Ax["MM", n + 1] + Ay["IM", m + 1],
                 V[n + 1, m + 1, "DG"] + Ax["DM", n + 1] + Ay["MM", m + 1],
                 V[n + 1, m + 1, "MM"] + Ax["MM", n + 1] + Ay["MM", m + 1],
                 V[n + 1, m + 1, "GD"] + Ax["MM", n + 1] + Ay["DM", m + 1],
                 V[n + 1, m + 1, "IM"] + Ax["IM", n + 1] + Ay["MM", m + 1])
      LLptr <- whichismax(LLcdt)
      score <- LLcdt[LLptr]
      z <- c(n + 1, m + 1, LLptr)
      while(z[1] > 1 | z[2] > 1){
        path <- c(z[3], path)
        progression <- cbind(z[1:2], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1,0,0), c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,0))
      }
    }else if (type == 'semiglobal'){
      tmp <- V[, , "MM"]
      tmp[1:n, 1:m] <- -Inf
      ind <- which(tmp == max(tmp), arr.ind = TRUE)
      if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
      z <- c(ind, 3)
      score <- V[z[1], z[2], z[3]]
      if(z[1] < n + 1){
        path <- rep(2, n + 1 - z[1])
        progression <- rbind(z[1]:n + 1, rep(z[2], n + 1 - z[1]))
      }else if(z[2] < m + 1){
        path <- rep(4, m + 1 - z[2])
        progression <- rbind(z[2]:m + 1, rep(z[1], m + 1 - z[2]))
      }
      while(z[1] > 1 & z[2] > 1){
        path <- c(z[3], path)
        progression <- cbind(z[1:2], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1,0,0), c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,0))
      }
      if(z[1] > 1){
        path <- c(rep(2, z[1] - 1), path)
        progression <- cbind(rbind(1:(z[1] - 1), rep(1, z[1] - 1)), progression)
      }else if(z[2] > 1){
        path <- c(rep(4, z[2] - 1), path)
        progression <- cbind(rbind(rep(1, z[2] - 1), 1:(z[2] - 1)), progression)
      }
    }else{
      ind <- which(V[, , "MM"] == max(V[, , "MM"]), arr.ind = TRUE)
      if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
      z <- c(ind, 3)
      score <- V[z[1], z[2], z[3]]
      P[1, 1, "MM"] <- 6
      while(P[z[1], z[2], z[3]] != 6){
        path <- c(z[3], path)
        progression <- cbind(z[1:2], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1,0,0), c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,0))
      }
    }
    key <- "MI = 1, DG = 2, MM = 3, GD = 4, IM = 5"
    rownames(progression) <- c(deparse(substitute(x)), deparse(substitute(y)))
  }else{
    A <- if(logspace) x$A else log(x$A)
    E <- if(logspace) x$E else log(x$E)
    Saa <- E - qe
    P[, 1, 1] <- c(NA, 2, rep(1, n - 1))
    P[1, , 3] <- c(NA, 2, rep(3, m - 1))
    V[1, 1, 2] <- 0
    if(type == "global"){ #key: D = 1, M = 2, I = 3
      V[-1, 1, 1] <- cumsum(c(0, A["DD", 2:n])) + A["MD", 1]
      V[1, -1, 3] <- seq(from = A["MI", 1], by = A["II", 1], length.out = m)
    }else{
      V[-1, 1, 1] <- V[1, -1, 3] <- 0 ### check this
    }
    for(i in if(is.null(itertab)) 1:n else 1:nrow(itertab)){
      for(j in if(is.null(itertab)) 1:m else 1){
        if(!is.null(itertab)){
          j <- itertab[i, 2]
          i <- itertab[i, 1]
        }
        if(pd){
          sij <- DNAprobC(y[j], Saa[, i]) + offset
        } else{
          sij <- Saa[y[j], i] + offset
        }
        Dcdt <- c(V[i, j + 1, "D"] + A["DD", i],
                  V[i, j + 1, "M"] + A["MD", i],
                  V[i, j + 1, "I"] + A["ID", i])
        Mcdt <- c(V[i, j, "D"] + A["DM", i] + sij,
                  V[i, j, "M"] + A["MM", i] + sij,
                  V[i, j, "I"] + A["IM", i] + sij,
                  if(type == "local") 0 else NULL)
        Icdt <- c(V[i + 1, j, "D"] + A["DI", i + 1],
                  V[i + 1, j, "M"] + A["MI", i + 1],
                  V[i + 1, j, "I"] + A["II", i + 1])
        Dmax <- whichismax(Dcdt)
        Mmax <- whichismax(Mcdt)
        Imax <- whichismax(Icdt)
        V[i + 1, j + 1, "D"] <- Dcdt[Dmax]
        V[i + 1, j + 1, "M"] <- Mcdt[Mmax]
        V[i + 1, j + 1, "I"] <- Icdt[Imax]
        P[i + 1, j + 1, "D"] <- Dmax
        P[i + 1, j + 1, "M"] <- Mmax
        P[i + 1, j + 1, "I"] <- Imax
      }
    }
    path <- c()
    progression <- c()
    if(type == 'global'){
      LLcdt <- c(V[n + 1, m + 1, "D"] + A["DM", n + 1],
                 V[n + 1, m + 1, "M"] + A["MM", n + 1],
                 V[n + 1, m + 1, "I"] + A["IM", n + 1])
      LLptr <- whichismax(LLcdt)
      score <- LLcdt[LLptr]
      z <- c(n + 1, m + 1, LLptr)
      while(z[1] > 1 | z[2] > 1){
        path <- c(z[3], path)
        progression <- c(z[1], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
      }
    }else if (type == 'semiglobal'){
      tmp <- V[, , "M"]
      tmp[1:n, 1:m] <- -Inf
      ind <- which(tmp == max(tmp), arr.ind = TRUE)
      if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
      z <- c(ind, 2)
      score <- V[z[1], z[2], z[3]]
      if(z[1] < n + 1){
        path <- rep(1, n + 1 - z[1])
        progression <- z[1]:n + 1
      }else if(z[2] < m + 1){
        path <- rep(3, m + 1 - z[2])
        progression <- rep(n + 1, m + 1 - z[2])
      }
      while(z[1] > 1 & z[2] > 1){
        path <- c(z[3], path)
        progression <- c(z[1], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
      }
      if(z[1] > 1){
        path <- c(rep(1, z[1] - 1), path)
        progression <- c(1:(z[1] - 1), progression)
      }else if(z[2] > 1){
        path <- c(rep(3, z[2] - 1), path)
        progression <- c(rep(1, z[2] - 1), progression)
      }
    }else{
      ind <- which(V[, , "M"] == max(V[, , "M"]), arr.ind = TRUE)
      if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
      z <- c(ind, 2)
      score <- V[z[1], z[2], z[3]]
      P[1, 1, "M"] <- 4
      while(P[z[1], z[2], z[3]] != 4){
        path <- c(z[3], path)
        progression <- c(z[1], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
      }
    }
    key <- "1 = delete, 2 = match, 3 = insert"
  }
  progression <- progression - 1 #to account for 0 row
  res <- structure(list(score = score,
                        path = path,
                        progression = progression,
                        key = key,
                        V = V,
                        pointer = P),
                   class = 'Viterbi')
  return(res)
}

#' @rdname Viterbi
Viterbi.HMM <- function (x, y, logspace = "autodetect"){
  if(identical(logspace, 'autodetect')) logspace <- logdetect(x)
  n <- length(y)
  states <- rownames(x$E)
  H <- length(states) # not including BeginEnd state
  path <- rep(NA, n)
  V <- array(-Inf, dim = c(H, n), dimnames = list(state = states, roll = 1:n))
  P <- V + NA # pointer array
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  #s <- if(logspace) x$s else log(x$s)
  V[, 1] <- E[, y[1]] + A[1, -1] # formerly s
  fun <- Vectorize(function(k, l) V[k, i - 1] + A[k + 1, l + 1])
  for (i in 2:n){
    tmp <- outer(1:H, 1:H, fun)
    V[, i] <- E[, y[i]] + apply(tmp, 2, max)
    P[, i] <- states[apply(tmp, 2, which.max)]
  }
  ak0 <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H)
  endstate <- which.max(V[, n] + ak0)
  maxLL <- V[endstate, n] + ak0[endstate]
  path[n] <- states[endstate]
  tmp <- path[n]
  for(i in n:2){
    path[i - 1] <- P[tmp, i]
    tmp <- path[i - 1]
  }
  res <- structure(list(score = maxLL,
                        path = path,
                        progression = NULL,
                        V = V,
                        pointer = P),
                   class = 'Viterbi')
  return(res)
}


#' @rdname Viterbi
Viterbi.DNAbin <- function(x, y, type = 'semiglobal', d = 8, e = 2,
                            S = NULL, itertab = NULL, offset = 0){
  ###check that x is a character vector###
  if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
  if(!inherits(y, "DNAbin")) stop("second argument (y) must also be of class 'DNAbin'")
  if(is.list(x)){
    if(length(x) == 1){
      tmp <- attributes(x)
      x <- x[[1]]
      attributes(x) <- tmp
    }else stop("Invalid input: multi-sequence list")
  }
  if(is.list(y)){
    if(length(y) == 1){
      tmp <- attributes(y)
      y <- y[[1]]
      attributes(y) <- tmp
    }else stop("Invalid input: multi-sequence list")
  }
  #class(x) <- class(y) <- NULL
  n <- length(x) + 1
  m <- length(y) + 1
  if(identical(itertab, "WilburLipman")) itertab = WilburLipman(x, y)
  if(!any(itertab)) itertab <- matrix(TRUE, n, m)
  if(is.null(S)){
    S <- NUC4.4
    #S <- diag(2, nrow = 4) - 1
  }else{
    if(nrow(S) != ncol(S) & nrow(S) != 15) stop("Invalid scoring matrix (S)")
    ### also check names order
  }
  guide <- as.raw(c(136, 24, 72, 40, 96, 144, 192, 48, 80, 160, 112, 224, 176, 208, 240))
  # initialize scoring and pointer arrays (M and P)
  M <- array(-Inf, dim = c(n, m, 3))
  P <- M + NA
  M[1, 1, 2] <- 0
  if(type == 'global'){
    M[, 1, 1] <- c(0, seq(from = -d, to = (- d + (n - 2) * -e), by = -e))
    M[1, , 3] <- c(0, seq(from = -d, to = (- d + (m - 2) * -e), by = -e))
  }else{
    M[2:n, 1, 1] <- M[1, 2:m, 3] <- 0 ### check this - should fill dim2 instead?
  }
  P[2:n, 1, 1] <- 1
  P[1, 2:m, 3] <- 3
  # recursion step
  # for(i in if(is.null(itertab)) 2:n else 1:nrow(itertab)){
  #   for(j in if(is.null(itertab)) 2:m else 1){
  #     if(!is.null(itertab)){
  #       j <- itertab[i, 2]
  #       i <- itertab[i, 1]
  #     }
  for(i in 2:n){
    for(j in 2:m){
      if(itertab[i, j]){
        ###
        #sij <-  S[x[i - 1], y[j - 1]] + offset
        sij <-  S[match(x[i - 1], guide), match(y[j - 1], guide)] + offset
        ###
        Ixcdt <- c(M[i - 1, j, 1] - e, M[i - 1, j, 2] - (d + e))#x alig to gap in y
        Mcdt <- c(M[i - 1, j - 1, 1] + sij,
                  M[i - 1, j - 1, 2] + sij,
                  M[i - 1, j - 1, 3] + sij,
                  if(type == 'local') 0 else NULL)
        Iycdt <- c(-Inf, M[i, j - 1, 2] - (d + e), M[i, j - 1, 3] - e)
        Ixmax <- whichismax(Ixcdt)
        Mmax <- whichismax(Mcdt)
        Iymax <- whichismax(Iycdt)
        M[i, j, 1] <- Ixcdt[Ixmax]
        M[i, j, 2] <- Mcdt[Mmax]
        M[i, j, 3] <- Iycdt[Iymax]
        P[i, j, 1] <- Ixmax
        P[i, j, 2] <- Mmax
        P[i, j, 3] <- Iymax
      }
    }
  }
  path <- c()
  progression <- matrix(nrow = 2, ncol = 0)
  if(type == 'global'){
    # find highest score in bottom right corner of scoring array M
    z <- c(n, m, whichismax(M[n, m, ]))
    score <- M[z[1], z[2], z[3]]
    while(z[1] > 1 | z[2] > 1){
      path <- c(z[3], path)
      progression <- cbind(z[1:2], progression)
      z[3] <- P[z[1], z[2], z[3]]
      z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
    }
    #alig <- trackback(z, condition = "z[1] > 1 | z[2] > 1", P = P)
  }else if(type == 'semiglobal'){
    # find highest score on bottom row or right column of scoring array M
    border <- rbind(M[n, , ], M[-n, m, ])
    ind <- which(border == max(border), arr.ind = TRUE)
    if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
    z <- if(ind[1] <= m) c(n, ind) else c(ind[1] - m, m, ind[2])
    if(is.na(P[z[1], z[2], z[3]])) stop("alignment failed, try increasing offset")
    score <- M[z[1], z[2], z[3]]
    if(z[1] < n){
      path <- rep(1, n - z[1])
      progression <- rbind((z[1] + 1):n, rep(z[2], n - z[1]))
    }else if(z[2] < m){
      path <- rep(3, m - z[2])
      progression <- rbind((z[2] + 1):m, rep(z[1], m - z[2]))
    }
    while(z[1] > 1 & z[2] > 1){
      path <- c(z[3], path)
      progression <- cbind(z[1:2], progression)
      z[3] <- P[z[1], z[2], z[3]]
      z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
    }
    if(z[1] > 1){
      path <- c(rep(1, z[1] - 1), path)
      progression <- cbind(rbind(1:(z[1] - 1), rep(1, z[1] - 1)), progression)
    }else if(z[2] > 1){
      path <- c(rep(3, z[2] - 1), path)
      progression <- cbind(rbind(rep(1, z[2] - 1), 1:(z[2] - 1)), progression)
    }
  }else if(type == 'local'){
    # find highest score in scoring array M
    ind <- which(M[, , 2] == max(M[, , 2]), arr.ind = TRUE)
    z <- c(ind, 2)
    score <- M[z[1], z[2], z[3]]
    P[1, 1, 2] <- 4
    while(P[z[1], z[2], z[3]] != 4){
      path <- c(z[3], path)
      progression <- cbind(z[1:2], progression)
      z[3] <- P[z[1], z[2], z[3]]
      z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
    }
  }
  key <- "1: x aligns to gap in y, 2: match, 3: y aligns to gap in x"
  progression <- progression - 1
  #rownames(progression) <- c(deparse(substitute(x)), deparse(substitute(y)))
  rownames(progression) <- c("x", "y") ### fix
  res <- structure(list(score = score,
                        #alignment = alig,
                        path = path,
                        #firstmatch = firstmatch,
                        progression = progression,
                        key = key,
                        V = M,
                        pointer = P),
                   class = 'Viterbi')
  return(res)
}



#' @rdname Viterbi
Viterbi.default <- function(x, y, type = 'semiglobal', d = 8, e = 2,
                            S = NULL, itertab = NULL, offset = 0){
  ###check that x is a character vector###
  if(!(type %in% c('global','semiglobal','local'))) stop("invalid type")
  if(mode(x) != "character" | mode(y) != "character") stop("x and y must be of mode 'character'")
  n <- length(x) + 1
  m <- length(y) + 1
  if(identical(itertab, "WilburLipman")) itertab = WilburLipman(x, y)
  if(!any(itertab)) itertab <- matrix(TRUE, n, m)
  if (is.null(S)) {
    residues <- unique(c(x, y))
    S <- diag(2, nrow = length(residues)) - 1
    dimnames(S) <- list(residues, residues)
  }
  # initialize scoring and pointer arrays (M and P)
  M <- array(-Inf, dim = c(n, m, 3))
  P <- M + NA
  M[1, 1, 2] <- 0
  if(type == 'global'){
    M[, 1, 1] <- c(0, seq(from = -d, to = (- d + (n - 2) * -e), by = -e))
    M[1, , 3] <- c(0, seq(from = -d, to = (- d + (m - 2) * -e), by = -e))
  }else{
    M[2:n, 1, 1] <- M[1, 2:m, 3] <- 0 ### check this - should fill dim2 instead?
  }
  P[2:n, 1, 1] <- 1
  P[1, 2:m, 3] <- 3
  # recursion step
  # for(i in if(is.null(itertab)) 2:n else 1:nrow(itertab)){
  #   for(j in if(is.null(itertab)) 2:m else 1){
  #     if(!is.null(itertab)){
  #       j <- itertab[i, 2]
  #       i <- itertab[i, 1]
  #     }
  for(i in 2:n){
    for(j in 2:m){
      if(itertab[i, j]){
        sij <-  S[x[i - 1], y[j - 1]] + offset
        Ixcdt <- c(M[i - 1, j, 1] - e, M[i - 1, j, 2] - (d + e))#x alig to gap in y
        Mcdt <- c(M[i - 1, j - 1, 1] + sij,
                  M[i - 1, j - 1, 2] + sij,
                  M[i - 1, j - 1, 3] + sij,
                  if(type == 'local') 0 else NULL)
        Iycdt <- c(-Inf, M[i, j - 1, 2] - (d + e), M[i, j - 1, 3] - e)
        Ixmax <- whichismax(Ixcdt) ### could improve using nnet::which.is.max###
        Mmax <- whichismax(Mcdt)
        Iymax <- whichismax(Iycdt)
        M[i, j, 1] <- Ixcdt[Ixmax]
        M[i, j, 2] <- Mcdt[Mmax]
        M[i, j, 3] <- Iycdt[Iymax]
        P[i, j, 1] <- Ixmax
        P[i, j, 2] <- Mmax
        P[i, j, 3] <- Iymax
      }
    }
  }
  path <- c()
  progression <- matrix(nrow = 2, ncol = 0)
  if(type == 'global'){
    # find highest score in bottom right corner of scoring array M
    z <- c(n, m, whichismax(M[n, m, ]))
    score <- M[z[1], z[2], z[3]]
    while(z[1] > 1 | z[2] > 1){
      path <- c(z[3], path)
      progression <- cbind(z[1:2], progression)
      z[3] <- P[z[1], z[2], z[3]]
      z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
    }
    #alig <- trackback(z, condition = "z[1] > 1 | z[2] > 1", P = P)
  }else if(type == 'semiglobal'){
    # find highest score on bottom row or right column of scoring array M
    border <- rbind(M[n, , ], M[-n, m, ])
    ind <- which(border == max(border), arr.ind = TRUE)
    if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
    z <- if(ind[1] <= m) c(n, ind) else c(ind[1] - m, m, ind[2])
    if(is.na(P[z[1], z[2], z[3]])) stop("alignment failed, try increasing offset")
    score <- M[z[1], z[2], z[3]]
    if(z[1] < n){
      path <- rep(1, n - z[1])
      progression <- rbind((z[1] + 1):n, rep(z[2], n - z[1]))
    }else if(z[2] < m){
      path <- rep(3, m - z[2])
      progression <- rbind((z[2] + 1):m, rep(z[1], m - z[2]))
    }
    while(z[1] > 1 & z[2] > 1){
      path <- c(z[3], path)
      progression <- cbind(z[1:2], progression)
      z[3] <- P[z[1], z[2], z[3]]
      z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
    }
    if(z[1] > 1){
      path <- c(rep(1, z[1] - 1), path)
      progression <- cbind(rbind(1:(z[1] - 1), rep(1, z[1] - 1)), progression)
    }else if(z[2] > 1){
      path <- c(rep(3, z[2] - 1), path)
      progression <- cbind(rbind(rep(1, z[2] - 1), 1:(z[2] - 1)), progression)
    }
  }else if(type == 'local'){
    # find highest score in scoring array M
    ind <- which(M[, , 2] == max(M[, , 2]), arr.ind = TRUE)
    z <- c(ind, 2)
    score <- M[z[1], z[2], z[3]]
    P[1, 1, 2] <- 4
    while(P[z[1], z[2], z[3]] != 4){
      path <- c(z[3], path)
      progression <- cbind(z[1:2], progression)
      z[3] <- P[z[1], z[2], z[3]]
      z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
    }
  }
  key <- "1: x aligns to gap in y, 2: match, 3: y aligns to gap in x"
  progression <- progression - 1
  rownames(progression) <- c(deparse(substitute(x)), deparse(substitute(y)))
  res <- structure(list(score = score,
                        #alignment = alig,
                        path = path,
                        #firstmatch = firstmatch,
                        progression = progression,
                        key = key,
                        V = M,
                        pointer = P),
                   class = 'Viterbi')
  return(res)
}
