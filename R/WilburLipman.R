#' Define search space for dynammic programming.
#'
#' Implement the algorithm of Wilbur & Lipman (1983) to define windowspace for
#' downstream dynammic programming applications such as the Viterbi, forward and
#' backward algorithms.
#' @param x,y integer vectors.
#' @param arity integer indicating the numbering system that x and y are coded in,
#' for example 4 for DNA, 20 for proteins
#' @param threshold the number of standard deviations above the mean
#' number of k-tuple matches
#' for all diagonals with k-tuple matches for a diagonal to be be considered
#' significant...


WilburLipman <- function(x, y, arity = "autodetect", k = 4, w = 20, threshold = 5){
  # x and y coded as integers starting from 0
  if(arity == "autodetect") arity <- max(c(x, y)) + 1
  N1 <- length(x)
  N2 <- length(y)
  xkm <- matrix(nrow = k, ncol = N1 - k + 1) # kmers for x
  ykm <- matrix(nrow = k, ncol = N2 - k + 1) # kmers for y
  for(i in 1:k){
    xkm[i, ] <- x[i:(N1 - (k - i))]
    ykm[i, ] <- y[i:(N2 - (k - i))]
  }
  tmp <- apply(ykm, 2, decimal, from = arity) + 1 ### poss bug here with auto-simplify
  pointer <- lapply(1:arity^k, function(x) which(tmp == x))
  S1 <- apply(xkm, 2, function(x) pointer[[decimal(x, from = arity) + 1]])
  diagbin <- unlist(mapply("-", S1, 1:ncol(xkm)))
  diags <- tabulate(diagbin + (N1 + 1), nbins = N1 + N2 + 1)
  #names(diags) <- -N1:N2
  #diags <- table(unlist(mapply("-", S1, 1:ncol(xkm)))) # could prob speed this up a lot
  if(length(diags) == 0) return(c(-length(x), length(y)))
  #sigdiags <- as.numeric(names(diags[diags > mean(diags) + threshold * sd(diags)]))
  sigdiags <- which(diags > mean(diags) + threshold * sd(diags)) - (N1 + 1)
  if(length(sigdiags) > 0){
    #windowspace <- c()
    windowspace <- c(sigdiags[1] - w, sigdiags[length(sigdiags)] + w)
    #for(i in sigdiags) windowspace <- c(windowspace, (i - w):(i + w))
    #return(outer(1:(N1 + 1), 1:(N2 + 1), function(x, y) (y - x) %in% windowspace))
    return(windowspace)
  }else{
    return(c(-length(x), length(y)))
    #return(matrix(FALSE, nrow = N1 + 1, ncol = N2 + 1))
    #return(NULL)
  }
}

