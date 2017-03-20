#' K-mer distance matrix computation.
#'
#' Computes the matrix of k-tuple distances between all pairwise comparisons
#'   of a set of sequences.
#'
#' @param x a matrix of aligned sequences or a list of unaligned sequences.
#'   Accepted modes are "character" and "raw" (for "DNAbin" and "AAbin" objects).
#' @param k integer representing the k-mer size to be used for calculating
#'   the distance matrix. Defaults to 5. Note that high values of k
#'   may be slow to compute and use a lot of memory due to the large numbers
#'   of calculations required, particularly when the residue alphabet is
#'   also large.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gapchar the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param ... further arguments to be passed to \code{"dist"}.
#' @return a distance matrix of class \code{"dist"}
#' @details this function uses the k-mer distance measure outlined in
#'   Edgar (2004) to compute a N x N distance matrix (where N is the
#'   number of sequences).
#' @author Shaun Wilkinson
#' @references
#'   Edgar RC (2004) Local homology recognition and distance measures in
#'   linear time using compressed amino acid alphabets.
#'   \emph{Nucleic Acids Research}, \strong{32}, 380-385.
#' @seealso \code{\link{mbed}} for leaner distance matrices
#' @examples
#'   ## compute the k-mer distance for the woodmouse dataset in the ape package
#'   ## with a k-mer size of 4
#'   library(ape)
#'   data(woodmouse)
#'   woodmouse <- woodmouse[, 47:962] ## trim gappy ends for true global alignment
#'   woodmouse.dist <- kdistance(woodmouse, k = 4)
#'   ## compute and plot the UPGMA tree
#'   woodmouse.tree <- as.dendrogram(hclust(woodmouse.dist, "average"))
#'   plot(woodmouse.tree)
#' @name kdistance
################################################################################
kdistance <- function(x, k = 5, residues = NULL, gapchar = "-", ...){
  UseMethod("kdistance")
}

#' @rdname kdistance
################################################################################
kdistance.AAbin <- function(x, k = 5, ...){
  if(is.matrix(x)) x <- unalign(x)
  #x <- lapply(x, compress.AA, alpha = compress)
  if(min(sapply(x, length)) < k) stop("minimum sequence length is less than k")
  arity <- if(k > 2) 6 else 20 # compress AA alphabet for high k values
  x <- .encodeAA(x, arity = arity, na.rm = TRUE)
  tuplecount <- function(y, k, arity){
    tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
    for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
    res <- apply(tuplemat, 2, .decimal, from = arity) + 1
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
kdistance.DNAbin <- function(x, k = 5, ...){
  if(is.matrix(x)) x <- unalign(x)
  if(min(sapply(x, length)) < k) stop("minimum sequence length is shorter than k")
  #x <- lapply(x, function(y) y[y != as.raw(2)]) # N's had too much impact
  x <- lapply(x, function(y) y[!(y %in% as.raw(c(2, 4, 240)))])
  counts <- .kcountDNA(x, k = k)
  denoms <- sapply(x, length) - k + 1
  freqs <- counts/denoms
  return(dist(freqs, ... = ...))
}

#' @rdname kdistance
################################################################################
kdistance.default <- function(x, k = 5, residues = NULL, gapchar = "-", ...){

  if(is.matrix(x)) x <- unalign(x, gapchar = gapchar)
  if(.isDNA(x)){
    class(x) <- "DNAbin"
    return(kdistance.DNAbin(x, k = k, ... = ...))
  }else if(.isAA(x)){
    class(x) <- "AAbin"
    return(kdistance.AAbin(x, k = k, ... = ...))
  }
  if(min(sapply(x, length)) < k) stop("minimum sequence length is less than k")
  residues <- .alphadetect(x, residues = residues, gapchar = gapchar)
  arity <- length(residues)
  if(k > 2 & arity >= 20) stop("Unable to calculate distance matrix for
                               large k and large alphabet size. If residues
                               are amino acids consider converting to AAbin
                               object for compression")
  modes <- lapply(x, mode)
  if(!(all(modes == "integer") | all(modes == "numeric"))){
    x <- .encodeCH(x, residues = residues, na.rm = TRUE)
  }
  tuplecount <- function(y, k, arity){
    tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
    for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
    res <- apply(tuplemat, 2, .decimal, from = arity) + 1
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
#' The \code{mbed} function takes a list of sequences and returns a matrix of
#'   distances to a subset of seed sequences using the method outlined
#'   in Blacksheilds et al. (2010).
#'
#' @param x a matrix of aligned sequences or a list of unaligned sequences.
#'   Accepted modes are "character" and "raw" (for "DNAbin" and "AAbin" objects).
#' @param k integer representing the k-mer size to be used for calculating
#'   the distance matrix. Defaults to 5. Note that high values of k
#'   may be slow to compute and use a lot of memory due to the large numbers
#'   of calculations required, particularly when the residue alphabet is
#'   also large.
#' @param seeds optional integer vector indicating which sequences should
#'   be used as the seed sequences. If \code{seeds = NULL} a set of
#'   log(N, 2)^2 non-identical sequences is randomly selected from the
#'   sequence set (where N is the number of sequences; see Blacksheilds et al.
#'   2010).
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gapchar the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @return returns a N x log(N, 2)^2 matrix of class "mbed" (where N is the
#'   number of sequences). The returned
#'   object has additional attributes including a the 'seeds' vector and
#'   a matrix of k-mer counts with one row for each sequence and 4^k columns.
#' @details This function deals with ambiguities by assigning counts proportionally.
#'   For example the motif ACRTG would assign the 5-mers ACATG and ACGTG counts of
#'   0.5 each. This algorithm is O N * 4^k in memory and time complexity so can be very
#'   slow and memory hungry for larger values of k. Further details of the embedding
#'   process can be found in Blacksheilds et al. (2010).
#' @author Shaun Wilkinson
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#' @seealso \code{\link{kdistance}} for full N x N distance matrix computation
#' @examples
#'   ## embed the woodmouse dataset from the ape package
#'   library(ape)
#'   data(woodmouse)
#'   woodmouse.mbed <- mbed(woodmouse)
#'   woodmouse.mbed
################################################################################
mbed <- function(x, seeds = NULL, k = 5, residues = NULL, gapchar = "-"){
  if(is.matrix(x)) x <- unalign(x, gapchar = gapchar)
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  nseq <- length(x)
  seqalongx <- seq_along(x)
  seqlengths <- sapply(x, length)
  if(min(seqlengths) < k) stop("minimum sequence length is less than k")
  if(DNA){
    x <- lapply(x, function(s) s[s != as.raw(2)])
    kcounts <- .kcountDNA(x, k = k)
  }else{
    tuplecount <- function(y, k, arity){
      tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
      for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
      res <- apply(tuplemat, 2, .decimal, from = arity) + 1
      res <- tabulate(res, nbins = arity^k)
      return(res)
    }
    if(AA){
      arity <- if(k > 2) 6 else 20 # compress AA alphabet for high k values
      x <- .encodeAA(x, arity = arity, na.rm = TRUE)
    }else{
      residues <- .alphadetect(x, residues = residues, gapchar = gapchar)
      arity <- length(residues)
      if(k > 2 & arity >= 20) stop("Unable to calculate distance matrix for
                               large k and large alphabet size. If residues
                               are amino acids consider converting to AAbin
                               object for compression")
      modes <- lapply(x, mode)
      if(!(all(modes == "integer") | all(modes == "numeric"))){
        x <- .encodeCH(x, residues = residues, na.rm = TRUE)
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
  res <- .kdist(kcounts, from = seqalongx - 1, to = seeds - 1,
                seqlengths = seqlengths, k = k)
  attr(res, "seeds") <- seeds
  attr(res, "kcounts") <- kcounts
  attr(res, "duplicates") <- duplicates
  class(res) <- "mbed"
  return(res)
}
################################################################################
