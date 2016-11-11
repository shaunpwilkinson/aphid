#' K-mer distance matrix calculation.
#'
#' Computes the distance matrix ...
#'
kdist <- function(x, k = 5, compress = "Dayhoff6", asmatrix = FALSE){
  UseMethod("kdist")
}

kdist.AAbin <- function(x, k = 5, compress = "Dayhoff6", ...){
  x <- lapply(x, compress.AA, alphabet = compress)
  arity <- switch(compress, "Dayhoff6" = 6) #placeholder
  tuplecount <- function(y, k, arity){
    tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
    for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
    res <- apply(tuplemat, 2, decimal, from = arity) + 1
    res <- tabulate(res, nbins = arity^k)
    return(res)
  }
  counts <- sapply(x, tuplecount, k, arity)
  denoms <- sapply(x, length) - k + 1
  freqs <- t(counts)/denoms
  return(dist(freqs, ... = ...))
}
