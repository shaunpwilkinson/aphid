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

