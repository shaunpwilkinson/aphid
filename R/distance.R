#' K-mer distance matrix calculation.
#'
#' Computes the matrix of k-tuple distances between all pairwise comparisons
#' of the input sequence set.
#'
#' @param x a list of sequences, possibly an object of class
#'   \code{"DNAbin"} or \code{"AAbin"}.
#' @param k an integer specifying the size of the k-tuples.
#' @param alpha the residue alphabet.
#' @param ... further arguments to be passed to \code{"dist"}.
#' @return a distance matrix of class \code{"dist"}
#' @name kdistance
#'
#'
kdistance <- function(x, k = 5, alpha = "Dayhoff6", ...){
  UseMethod("kdistance")
}


#' @rdname kdistance
#'
kdistance.AAbin <- function(x, k = 5, alpha = "Dayhoff6", ...){
  #x <- lapply(x, compress.AA, alpha = compress)
  if(min(sapply(x, length)) < k) stop("minimum sequence length is shorter than k")
  x <- encode.AA(x, arity = 6, na.rm = TRUE)
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


#' @rdname kdistance
#'
kdistance.DNAbin <- function(x, k = 5, alpha = NULL, ...){
  if(min(sapply(x, length)) < k) stop("minimum sequence length is shorter than k")
  x <- lapply(x, function(y) y[y != as.raw(2)])
  counts <- kcount.DNA(x, k = k)
  denoms <- sapply(x, length) - k + 1
  freqs <- counts/denoms
  return(dist(freqs, ... = ...))
}


#' @rdname kdistance
#'
kdistance.default <- function(x, k = 5, alpha = "autodetect", ...){
  if(is.DNA(x)){
    return(kdistance.DNAbin(x, k = k, alpha = NULL, ... = ...))
  }else if(is.AA(x)){
    return(kdistance.DNAbin(x, k = k, alpha = "Dayhoff6", ... = ...))
  }
  if(min(sapply(x, length)) < k) stop("minimum sequence length is shorter than k")
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


#' K-mer counting for DNA.
#'
#' \code{kcountDNA} takes a list of DNA sequences and returns a matrix
#' of k-mer counts with one row for each sequence and 4^k columns.
#'
#' @param x a list of DNA sequences in the "DNAbin" format.
#' @param k integer representing the length of the k-mers to be counted.
#'
#' @details
#' This function deals with ambiguities by assigning counts proportionally.
#' For example the motif ACRTG would assign the 5-mers ACATG and ACGTG counts of
#' 0.5 each
#'
#' This algorithm is of the order n * 4^k in memory and time so can be very
#' slow and memory hungry for large values of k (> 8).
#'
#' @author Shaun P. Wilkinson
#'
#'
kcount.DNA <- function(x, k = 5){
  if(!is.DNA(x)) stop("sequence list must be a 'DNAbin' object")
  if(!(is.integer(k) | is.numeric(k))) stop("k must be an integer")
  k <- as.integer(k)
  x <- lapply(x, function(y) y[y != as.raw(2)])
  return(.kcountDNA(x, k = k))
}
