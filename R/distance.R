#' K-mer distance matrix calculation.
#'
#' Computes the distance matrix ...
#'
#' @param x a list of sequences, can be an object of class \code{"DNAbin"} or \code{"AAbin"}.
#' @param k an integer specifying the size of the k-tuples.
#' @param alpha the residue alphabet.
#' @param ... further arguments to be passed to \code{"dist"}.
#'
kdistance <- function(x, k = 5, alpha = "Dayhoff6", ...){
  UseMethod("kdistance")
}

kdistance.AAbin <- function(x, k = 5, alpha = "Dayhoff6", ...){
  #x <- lapply(x, compress.AA, alpha = compress)
  x <- compress.AA(x, na.rm = TRUE)
  arity <- switch(alpha, "Dayhoff6" = 6) #placeholder
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

kdistance.DNAbin <- function(x, k = 5, alpha = NULL, ...){
  counts <- kcount_DNA(x, k = k)
  denoms <- sapply(x, length) - k + 1
  freqs <- counts/denoms
  return(dist(freqs, ... = ...))
}

kdistance.default <- function(x, k = 5, alpha = "autodetect", ...){
  if(identical(alpha, "autodetect")) alpha <- unique(unlist(x))
  arity <- length(alpha)
  modes <- lapply(x, mode)
  if(!(all(modes == "integer") | all(modes == "numeric"))){
    fun <- function(y){
      ints <- 0:(arity - 1)
      res <- ints[match(y, alpha)]
      res[!is.na(res)]
    }
    x <- lapply(x, fun)
  }
  tuplecount <- function(y, k, arity){
    tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
    for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
    res <- apply(tuplemat, 2, decimal, from = arity) + 1
    tabulate(res, nbins = arity^k)
  }
  counts <- sapply(x, tuplecount, k, arity)
  denoms <- sapply(x, length) - k + 1
  freqs <- t(counts)/denoms
  return(dist(freqs, ... = ...))
}



