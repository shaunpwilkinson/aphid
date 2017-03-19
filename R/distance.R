#' K-mer distance matrix computation.
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
################################################################################
kdistance <- function(x, k = 5, alpha = "Dayhoff6", ...){
  UseMethod("kdistance")
}

#' @rdname kdistance
################################################################################
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
################################################################################
kdistance.DNAbin <- function(x, k = 5, alpha = NULL, ...){
  if(min(sapply(x, length)) < k) stop("minimum sequence length is shorter than k")
  x <- lapply(x, function(y) y[y != as.raw(2)])
  counts <- aphid:::.kcountDNA(x, k = k)
  denoms <- sapply(x, length) - k + 1
  freqs <- counts/denoms
  return(dist(freqs, ... = ...))
}

#' @rdname kdistance
################################################################################
kdistance.default <- function(x, k = 5, ...){
  if(is.DNA(x)){
    return(kdistance.DNAbin(x, k = k, alpha = NULL, ... = ...))
  }else if(is.AA(x)){
    return(kdistance.AAbin(x, k = k, alpha = "Dayhoff6", ... = ...))
  }
  if(min(sapply(x, length)) < k) stop("minimum sequence length is shorter than k")
  if(length(alpha) < 2) alpha <- unique(unlist(x))
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
    res <- tabulate(res, nbins = arity^k)
    return(res)
  }
  counts <- sapply(x, tuplecount, k, arity)
  denoms <- sapply(x, length) - k + 1
  freqs <- t(counts)/denoms
  return(dist(freqs, ... = ...))
}
################################################################################

#' Convert sequences to vectors of distances to \emph{n} seed sequences.
#'
#' \code{embed.DNA} takes a list of sequences and returns a matrix of
#'   distances to a subset of seed sequences using the method outlined
#'   in Blacksheilds et al. (2010).
#'
#' @param x a list of sequences, either character vectors or raw vectors in
#'   "DNAbin" or "AAbin" format.
#' @param k integer representing the length of the k-mers to be counted.
#' @param seeds integer vector indicating which sequences should be used as
#'   the seed sequences.
#' @param alpha the compression alphabet to use if the input sequence list
#'   is an "AAbin" object. Currently the only supported options are "Dayhoff6"
#'   and "none". Note that if \code{alpha = "none"} or the input sequences have
#'   20 or more unique residues (e.g. protein sequences),
#'   the value of k should be kept low (k < 4)
#'   otherwise the number of k-mers to count will become excessively high.
#' @return returns a N x  matrix of class "embed". The returned object also has
#'   attributes including a matrix of k-mer counts with one row for each
#'   sequence and 4^k columns.
#' @details This function deals with ambiguities by assigning counts proportionally.
#'   For example the motif ACRTG would assign the 5-mers ACATG and ACGTG counts of
#'   0.5 each.
#'   This algorithm is order n * 4^k in memory and time complexity so can be very
#'   slow and memory hungry for large values of k (> 8).
#' @author Shaun Wilkinson
#'
################################################################################
mbed <- function(x, seeds = NULL, k = 5, alpha = "Dayhoff6"){
  if(!is.list(x)) stop("x must be a list")
  DNA <- is.DNA(x)
  AA <- is.AA(x)
  nseq <- length(x)
  seqalongx <- seq_along(x)
  seqlengths <- sapply(x, length)
  if(min(seqlengths) < k) stop("minimum sequence length is shorter than k")
  if(DNA){
    x <- lapply(x, function(s) s[s != as.raw(2)])
    kcounts <- aphid:::.kcountDNA(x, k = k)
  }else{
    tuplecount <- function(y, k, arity){
      tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
      for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
      res <- apply(tuplemat, 2, decimal, from = arity) + 1
      res <- tabulate(res, nbins = arity^k)
      return(res)
    }
    if(AA){
      arity <- switch(alpha, "Dayhoff6" = 6, "none" = 20, stop("invalid 'alpha'"))
      #arity <- if(alpha == "Dayhoff6") 6 else if(alpha == "none") 20 else stop("invalid alpha")
      x <- encode.AA(x, arity = arity, na.rm = TRUE)
    }else{
      if(length(alpha) < 2) alpha <- unique(unlist(x))
      arity <- length(alpha)
      modes <- lapply(x, mode)
      if(!(all(modes == "integer") | all(modes == "numeric"))){
        encode.char <- function(y){
          ints <- 0:(arity - 1)
          res <- ints[match(y, alpha)]
          res[!is.na(res)]
        }
        x <- lapply(x, encode.char)
      }
    }
    kcounts <- t(sapply(x, tuplecount, k, arity))
  }
  if(is.null(seeds)){
    duplicates <- duplicated(x)
    nseeds <- min(sum(!duplicates), ceiling(log(nseq, 2)^2))
    # LLR algorithm see Blacksheilds et al. 2010
    seeds <- sort(sample(seqalongx[!duplicates], size = nseeds))
    ## TODO could expand on this using pivot objects etc
  }
  res <- aphid:::.kdist(kcounts, from = seqalongx - 1, to = seeds - 1,
                seqlengths = seqlengths, k = k)
  attr(res, "seeds") <- seeds
  attr(res, "kcounts") <- kcounts
  attr(res, "duplicates") <- duplicates
  class(res) <- "mbed"
  return(res)
}
################################################################################
