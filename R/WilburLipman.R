#' Define search space for dynammic programming.
#' 
#' Implement the algorithm of Wilbur & Lipman (1983) to define windowspace for 
#' downstream dynammic programming applications such as the Viterbi, forward and 
#' backward algorithms.
#' 
#' @param threshold the number of standard deviations above the mean
#' number of k-tuple matches 
#' for all diagonals with k-tuple matches for a diagonal to be be considered 
#' significant...


WilburLipman <- function(x, y, k = 4, w = 10, threshold = 5){
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
  N1 <- length(x)
  N2 <- length(y)
  n2t <- function(n) switch(n, 'a' = 0, 'c' = 1, 'g'= 2, 't' = 3, sample(0:3, 1)) ### hack
  xtr <- if(inherits(x, "DNAbin")) as.ternary.DNAbin(x) else unname(sapply(x, n2t))
  ytr <- if(inherits(y, "DNAbin")) as.ternary.DNAbin(y) else unname(sapply(y, n2t))
  xkm <- matrix(nrow = k, ncol = N1 - k + 1) # kmers for x
  ykm <- matrix(nrow = k, ncol = N2 - k + 1) # kmers for y
  for(i in 1:k){
    xkm[i, ] <- xtr[i:(N1 - (k - i))]
    ykm[i, ] <- ytr[i:(N2 - (k - i))]
  }
  tmp <- apply(ykm, 2, decimal, from = 4) + 1
  pointer <- lapply(1:k^4, function(x) which(x == tmp))
  S1 <- apply(xkm, 2, function(x) pointer[[decimal(x, from = 4) + 1]])
  diags <- table(unlist(mapply("-", S1, 1:ncol(xkm)))) # could speed this up a lot with cpp
  if(length(diags) == 0) return(NULL)
  sigdiags <- as.numeric(names(diags[diags > mean(diags) + threshold * sd(diags)]))
  if(length(sigdiags) > 0){
    windowspace <- c()
    for(i in sigdiags) windowspace <- c(windowspace, (i - w):(i + w))
    return(outer(1:(N1 + 1), 1:(N2 + 1), function(x, y) (y - x) %in% windowspace))
  }else{
    return(matrix(FALSE, nrow = N1 + 1, ncol = N2 + 1))
  } 
}


# WilburLipman <- function(x, y, k = 4, d = 6, e = 0, w = 10, S = NULL){
#   N1 <- length(x)
#   N2 <- length(y)
#   if (is.null(S)) {
#     residues <- unique(c(x, y))
#     S <- diag(2, nrow = length(residues))
#     S <- S - 1
#     dimnames(S) <- list(residues, residues)
#   }
#   n2t <- function(x) switch(x, 'a' = 0, 'c' = 1, 'g'= 2, 't' = 3)
#   xtr <- unname(sapply(x, n2t)) #x expressed as a ternary vector
#   ytr <- unname(sapply(y, n2t)) #y expressed as a ternary vector
#   xnc <- N1 - k + 1 # number of columns for the k-row x matrix
#   ync <- N2 - k + 1 # number of columns for the k-row y matrix
#   xsq <- 1:xnc
#   ysq <- 1:ync
#   xdf <- xtr[xsq]
#   ydf <- ytr[ysq]
#   for(i in 2:k){
#     xsq <- xsq + 1
#     ysq <- ysq + 1
#     xdf <- rbind(xdf, xtr[xsq])
#     ydf <- rbind(ydf, ytr[ysq])
#   }
#   tmp <- apply(ydf, 2, decimal, from = 4) + 1
#   pointer <- apply(matrix(1:k^4, nrow = 1), 2, function(x) which(x == tmp))
#   S1 <- apply(xdf, 2, function(x) pointer[[decimal(x, from = 4) + 1]])
#   diags <- table(unlist(mapply("-", S1, 1:xnc)))
#   if(length(diags) == 0) return(NULL)
#   possdiags <- xnc + ync + 1
#   poismean <- mean(c(rep(0, possdiags - length(diags)), diags))
#   pvals <- sapply(diags, function(x) poisson.test(x,  r = poismean)$p.value)
#   sigdiags <- as.integer(names(pvals[pvals < 0.001/possdiags]))
#   if(length(sigdiags) == 0) return(NULL)
#   windowspace <- c()
#   for(i in sigdiags) windowspace <- c(windowspace, (i - w):(i + w))
#   tmp <- outer(1:N1, 1:N2, function(x, y) (y - x) %in% windowspace)
#   tmp <- cbind(rep(FALSE, N1), tmp)
#   tmp <- rbind(rep(FALSE, N2 + 1), tmp)
#   itertab <- which(tmp, arr.ind = TRUE)
#   itertab <- itertab[order(itertab[, 1], itertab[, 2]), ]
#   # initialize scoring and pointer arrays (M and P)
#   n <- N1 + 1
#   m <- N2 + 1
#   M <- array(-Inf, dim = c(n, m, 3))
#   M[1, 1, 2] <- 0
#   P <- array(NA, dim = c(n, m, 3))
#   M[2:n, 1, 1] <- 0
#   M[1, 2:m, 3] <- 0
#   P[2:n, 1, 1] <- 1
#   P[1, 2:m, 3] <- 3
#   # recursion step
#   for(h in 1:nrow(itertab)){
#     i <- itertab[h, 1]
#     j <- itertab[h, 2]
#     # match state
#     sij <- S[x[i - 1], y[j - 1]]
#     Mcandidates <- c(M[i - 1, j - 1, 1] + sij,
#                      M[i - 1, j - 1, 2] + sij,
#                      M[i - 1, j - 1, 3] + sij)
#     M[i, j, 2] <- max(Mcandidates)
#     Mpointer <- which(Mcandidates == max(Mcandidates))
#     if(length(Mpointer) > 1) Mpointer <- sample(Mpointer, size = 1)
#     P[i, j, 2] <- Mpointer
#     # x aligned to gap in y
#     Ixcandidates <- c(M[i - 1, j, 1] - e, M[i - 1, j, 2] - (d + e))
#     M[i, j, 1] <- max(Ixcandidates)
#     if(!(all(Ixcandidates == -Inf))){
#       Ixpointer <- which(Ixcandidates == max(Ixcandidates))
#       if(length(Ixpointer) > 1) Ixpointer <- sample(Ixpointer, size = 1)
#       P[i, j, 1] <- Ixpointer
#     }
#     # y aligned to a gap in x
#     Iycandidates <- c(-Inf, M[i, j - 1, 2] - (d + e), M[i, j - 1, 3] - e)
#     M[i, j, 3] <- max(Iycandidates)
#     if(!(all(Iycandidates == -Inf))){
#       Iypointer <- which(Iycandidates == max(Iycandidates))
#       if(length(Iypointer) > 1) Iypointer <- sample(Iypointer, size = 1)
#       P[i, j, 3] <- Iypointer
#     }
#   }
#   # initialize result alignment (matrix)
#   alig <- matrix(c(NA, NA))
#   # find highest score
#   bottomright <- rbind(M[n, , ], M[-n, m, ])
#   max.ind <- which(bottomright == max(bottomright), arr.ind = TRUE)
#   tbm <- max.ind[2]
#   if(max.ind[1] <= m){
#     tbr <- n
#     tbc <- max.ind[1]
#     if(tbc != m){
#       rightend <- rbind(rep("-", m - tbc), y[tbc:(m - 1)])
#     }else{
#       rightend <- NULL
#     }
#   }else{
#     tbr <- max.ind[1] - m
#     tbc <- m
#     if(tbr != n){
#       rightend <- rbind(x[tbr:(n - 1)], rep("-", n - tbr))
#     }else{
#       rightend <- NULL
#     }
#   }
#   score <- M[tbr, tbc, tbm] # score of optimal alignment
#   # traceback
#   while(!(tbr == 1 | tbc == 1)){
#     app <- switch(tbm,
#                   c(x[tbr - 1], "-"),
#                   c(x[tbr - 1], y[tbc - 1]),
#                   c("-", y[tbc - 1]))
#     alig <- cbind(alig, unname(app))
#     temp <- tbm
#     tbm <- P[tbr, tbc, tbm]
#     if(temp == 1){
#       tbr <- tbr - 1
#     } else if (temp == 2){
#       tbr <- tbr - 1
#       tbc <- tbc - 1
#     } else if(temp == 3){
#       tbc <- tbc - 1
#     }
#   }
#   if(tbr == 1 & tbc == 1){
#     leftend <- NULL
#   }else if(tbr == 1){
#     leftend <- rbind(rep("-", tbc - 1), y[1:(tbc - 1)])
#   }else{
#     leftend <- rbind(x[1:(tbr - 1)], rep("-", tbr - 1))
#   }
#   # reverse alignment matrix
#   alig <- alig[, ncol(alig):2]
#   rownames(alig) <- c(deparse(substitute(x)), deparse(substitute(y)))
#   globalalig <- alig
#   alig <- cbind(leftend, alig, rightend)
#   res <- structure(list(score = score, globalAlignment = globalalig,
#                         alignment = alig, ViterbiPath = NULL,
#                         ViterbiArray = M, PointerArray = P),
#                    class = 'alignment')
#   res
# }
