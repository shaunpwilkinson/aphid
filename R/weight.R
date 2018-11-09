#' Sequence weighting.
#'
#' Weighting schemes for DNA and amino acid sequences.
#'
#' @param x a list or matrix of sequences
#'   (usually a "DNAbin" or "AAbin" object).
#'   Alternatively x can be an object of class \code{"dendrogram"}
#'   for tree-base weighting.
#' @param method a character string indicating the weighting method to be used.
#'   Currently the only methods available are a modified version of the
#'   maximum entropy weighting scheme proposed by
#'   Henikoff and Henikoff (1994) (\code{method = "Henikoff"})
#'   and the tree-based weighting scheme of Gerstein et al (1994)
#'   (\code{method = "Gerstein"}).
#' @param k integer representing the k-mer size to be used.
#'   Defaults to 5. Note that higher
#'   values of k may be slow to compute and use excessive memory due to
#'   the large numbers of calculations required.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Furthermore, the default setting
#'   \code{residues = NULL} will not detect rare residues that are not present
#'   in the sequences, and thus will not assign them emission probabilities
#'   in the model. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param ... additional arguments to be passed between methods.
#' @return a named vector of weights, the sum of which is equal to
#'    the total number of sequences (average weight = 1).
#' @details
#'   This is a generic function.
#'   If \code{method = "Henikoff"} the sequences are weighted
#'   using a modified version of the maximum entropy method proposed by
#'   Henikoff and Henikoff (1994). In this case the
#'   maximum entropy weights are calculated from a k-mer presence absence
#'   matrix instead of an alignment as originally described by
#'   Henikoff and Henikoff (1994).
#'   If \code{method = "Gerstein"} the agglomerative method of
#'   Gerstein et al (1994) is used to weight sequences based
#'   on their relatedness as derived from a phylogenetic tree.
#'   In this case a dendrogram is first derived using the
#'   \code{\link[kmer]{cluster}} function in the
#'   \code{\link[kmer]{kmer}} package.
#'   Methods are available for
#'   \code{"dendrogram"} objects, \code{"DNAbin"} and \code{"AAbin"}
#'   sequence objects (as lists or matrices) and sequences in standard
#'   character format provided either as lists or matrices.
#'
#'   For further details on sequence weighting schemes see Durbin et al
#'   (1998) chapter 5.8.
#'
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Gerstein M, Sonnhammer ELL, Chothia C (1994) Volume changes in protein evolution.
#'   \emph{Journal of Molecular Biology}, \strong{236}, 1067-1078.
#'
#'   Henikoff S, Henikoff JG (1994) Position-based sequence weights.
#'   \emph{Journal of Molecular Biology}, \strong{243}, 574-578.
#'
#' @examples
#'   ## weight the sequences in the woodmouse dataset from the ape package
#'   library(ape)
#'   data(woodmouse)
#'   woodmouse.weights <- weight(woodmouse)
#'   woodmouse.weights
#' @name weight
################################################################################
weight <- function(x, ...){
  UseMethod("weight")
}
################################################################################
#' @rdname weight
################################################################################
weight.DNAbin <- function(x, method = "Henikoff", k = 5, ...){
  if(is.list(x)){
    weight.list(x, method = method, k = k)
  }else{
    x <- unalign(x, gap = as.raw(4))
    weight.list(x, method = method, k = k)
  }
}
################################################################################
#' @rdname weight
################################################################################
weight.AAbin <- function(x, method = "Henikoff", k = 5, ...){
  if(is.list(x)){
    weight.list(x, method = method, k = k)
  }else{
    x <- unalign(x, gap = as.raw(45))
    weight.list(x, method = method, k = k)
  }
}
################################################################################
#' @rdname weight
################################################################################
weight.list <- function(x, method = "Henikoff", k = 5, residues = NULL,
                        gap = "-", ...){
  nsq <- length(x)
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  for(i in 1:nsq) x[[i]] <- x[[i]][x[[i]] != gap]
  if(nsq > 2){
    if(identical(method, "Henikoff")){
      hashes <- .digest(x)
      pointers <- .point(hashes)
      catchnames <- names(x)
      x <- x[!duplicated(pointers)]
      nsq <- length(x)
      suppressMessages(kmers <- round(kmer::kcount(x, k = k, residues = residues, gap = gap)))
      kmers <- kmers > 0
      ksums <- apply(kmers, 2, sum)
      uwfs <- ksums/nsq # unweighted freqs
      ftweights <- rbind(1/(2*(1-uwfs)), 1/(2*uwfs))
      f <- function(a, b) mean(b[rbind(!a, a)])
      res <- apply(kmers, 1, f, ftweights)
      rm(kmers)
      spl <- split(seq_along(hashes), f = factor(pointers))
      ## divide weight evenly among duplicates
      res <- res/vapply(spl, length, 0L)
      res <- res[pointers]
      names(res) <- catchnames
      res <- res * length(res)/sum(res)
    }else if(identical(method, "Gerstein")){
      tmpnames <- names(x)
      names(x) <- paste0("S", 1:nsq)
      guidetree <-suppressMessages(kmer::cluster(x, k = k, residues = residues, gap = gap))
      res <- weight.dendrogram(guidetree, method = "Gerstein")[names(x)]
      names(res) <- tmpnames
    }else stop("Accepted methods are 'Gerstein' and 'Henikoff'\n")
  }else if(nsq == 2){
    res <- structure(c(1, 1), names = names(x))
  }else if(nsq == 1){
    res <- structure(1, names = names(x))
  }else{
    res <- numeric(0)
  }
  return(res)
}
################################################################################
#' @rdname weight
################################################################################
weight.dendrogram <- function(x, method = "Gerstein", ...){
  if(!identical(method, "Gerstein")) stop("Only Gerstein method supported")
  if(is.leaf(x)) return(structure(1, names = attr(x, "label")))
  acal <- function(d) !any(sapply(d, is.list)) # all children are leaves?
  md <- function(d) all(sapply(d, acal)) & !acal(d) # mergable dendro?
  ch <- function(d) sapply(d, attr, "height") # child heights
  if(acal(x)){
    res <- attr(x, "height") - ch(x) + 0.0000001
    res <- length(res) * (res/sum(res))
    names(res) <- sapply(x, attr, "label")
    return(res)
  }
  Gerstein <- function(x){ # x is a dendrogram
    ngrandchildren <- sapply(x, length)
    childisdendro <- ngrandchildren > 1
    if(md(x)){
      childheights <- ch(x[childisdendro])
      childedges <- attr(x, "height") - childheights
      # ch only works on dendro lists
      grandchildheights <- lapply(x[childisdendro], ch)
      #  list same length as childedges
      grandchildedges <- mapply("-", childheights, grandchildheights,
                                SIMPLIFY = FALSE)
      grandchildedges <- lapply(grandchildedges, function(e) e + 0.0000001)
      # this just safeguards against 0 denominators (but is a bit of a hack)
      ratios <- lapply(grandchildedges, function(v) v/sum(v))
      inheritances <- mapply("*", childedges, ratios, SIMPLIFY = FALSE)
      newgrandchildedges <- mapply("+", grandchildedges, inheritances,
                                   SIMPLIFY = FALSE)
      lcounter <- 1 #leaf counter
      dcounter <- 1 #dendro counter
      tmp <- x
      for(i in seq_along(ngrandchildren)){
        if(ngrandchildren[i] > 1){
          for(j in 1:ngrandchildren[i]){
            leafj <- tmp[[i]][[j]]
            attr(leafj, "height") <- attr(tmp, "height") -
              newgrandchildedges[[dcounter]][j]
            x[[lcounter]] <- leafj
            lcounter <- lcounter + 1
          }
          dcounter <- dcounter + 1
        }else{
          x[[lcounter]] <- tmp[[i]]
          lcounter <- lcounter + 1
        }
      }
    }else{
      x[childisdendro] <- lapply(x[childisdendro], Gerstein)
    }
    x
  }
  while(any(sapply(x, is.list))) x <- Gerstein(x)
  res <- sapply(x, function(d) attr(x, "height") - attr(d, "height"))
  res <- res * length(res)/sum(res)
  names(res) <- sapply(x, function(d) attr(d, "label"))
  res
}
################################################################################
#' @rdname weight
################################################################################
weight.default <- function(x, method = "Henikoff", k = 5, residues = NULL,
                           gap = "-", ...){
  x <- unalign(x, gap = gap)
  weight.list(x, method = method, k = k, residues = residues, gap = gap)
}
################################################################################
