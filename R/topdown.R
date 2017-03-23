#' Top down phylogenetic trees.
#'
#' Create phylogenetic trees by successively splitting the sequence dataset
#'   into smaller and smaller subsets.
#'
#' @param x a list of sequences, possibly an object of class
#'   \code{"DNAbin"} or \code{"AAbin"}.
#' @param seeds an integer vector indicating which sequences should be used as
#'   seeds for the clustering procedure (see details section).
#'   If NULL the seeds are randomly selected from the input sequence list.
#' @param k integer representing the the k-mer size to be used for calculating
#'   the distance matrix.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless the sequence list is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param weighted logical indicating whether the edge lengths in the
#'   tree should correspond to the degree of relatedness
#'   between the sequences. Defaults to TRUE (slower).
#' @return an object of class \code{"dendrogram"}.
#' @details This function creates phylogenetic trees by successively splitting
#'   the sequence dataset into smaller and smaller subsets (a.k.a. recursive
#'   partitioning). This is a "top down" approach to tree-building, as
#'   opposed to agglomerative "bottom up" methods such as neighbour joining
#'   and UPGMA. It is useful for large sequence sets (> 20, 000) since the
#'   usual N x N distance matrix is not required (where N is the number of
#'   sequences). Instead, an N x log(N, 2)^2
#'   distance matrix is first derived from the input sequences using the
#'   \code{\link{mbed}} function of Blacksheilds et al. 2010 and the k-mer
#'   istance measure of (Edgar 2004).
#'   This requires substantially less memory and time. The 'embedded' sequences
#'   are then split recursively using the k-means (clusters = 2) technique.
#' @author Shaun Wilkinson
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#'
#'   Edgar RC (2004) Local homology recognition and distance measures in
#'   linear time using compressed amino acid alphabets.
#'   \emph{Nucleic Acids Research}, \strong{32}, 380-385.
#' @examples
#' \dontrun{
#' ## Cluster the woodmouse dataset from the ape package
#' library(ape)
#' data(woodmouse)
#' woodmice <- unalign(woodmouse[, 47:962]) ## trim gappy ends
#' set.seed(999)
#' woodmice.tree <- topdown(woodmice)
#' par(mar = c(5, 2, 4, 8) + 0.1)
#' plot(woodmice.tree, main = "Woodmouse tree", horiz = TRUE)
#' }
################################################################################
topdown <- function(x, seeds = NULL, k = 5, residues = NULL, gap = "-",
                    weighted = TRUE){
  if(is.matrix(x)) x <- unalign(x, gap = gap)
  # first embed the seqs in a N x log(N, 2)^2 distmat as in Blackshields 2010
  nseq <- length(x)
  seqlengths <- sapply(x, length)
  #generate long thin distance matrix with counts attributes etc
  M <- mbed(x, seeds = seeds, k = k, residues = residues, gap = gap)
  # throw some minor random variation into duplicate rows so kmeans doesn't crash
  Mduplicates <- duplicated.matrix(M, MARGIN = 1)
  M[Mduplicates, ] <- M[Mduplicates, ] + runif(length(M[Mduplicates, ]),
                                                min = -0.0001, max = 0.0001)
  seeds <- attr(M, "seeds")
  nseeds <- length(seeds)
  kcounts <- attr(M, "kcounts")
  # initialize the tree
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "sequences") <- 1:nseq
  attr(tree, "height") <- 0
  # define recursive splitting functions
  topdown1 <- function(tree, d){
    tree <- topdown2(tree, d = d)
    if(is.list(tree)) tree[] <- lapply(tree, topdown1, d = d)
    return(tree)
  }
  topdown2 <- function(node, d){
    if(!is.list(node)){ # fork leaves only
      seqs <- d[attr(node, "sequences"), , drop = FALSE]
      nseq <- nrow(seqs)
      if(nseq == 1) return(node)
      membership <- if(nseq == 2) 1:2 else kmeans(seqs, centers = 2)$cluster
      tmpattr <- attributes(node)
      node <- vector(mode = "list", length = 2)
      attributes(node) <- tmpattr
      attr(node, "leaf") <- NULL
      for(i in 1:2){
        node[[i]] <- 1
        attr(node[[i]], "height") <- attr(node, "height") - 1
        attr(node[[i]], "leaf") <- TRUE
        attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
      }
    }
    return(node)
  }
  # recursively build the tree
  tree <- topdown1(tree, d = M)
  tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  # model the branch lengths
  if(weighted){
    avdist <- function(node){
      if(is.list(node)){
        tof <- matrix(TRUE, nrow = nseq, ncol = nseeds)
        rownames(tof) <- names(x)
        colnames(tof) <- names(x)[seeds]
        node1seqs <- attr(node[[1]], "sequences")
        node2seqs <- attr(node[[2]], "sequences")
        n1sis <- node1seqs %in% seeds
        n2sis <- node2seqs %in% seeds
        if(any(n1sis) & any(n2sis)){
          tof[-node1seqs[n1sis], ] <- FALSE
          tof[, -(match(node2seqs[n2sis], seeds))] <- FALSE
          attr(node, "avdist") <- mean(M[tof])
        }else{
          dists <- .kdist(kcounts, from = node1seqs - 1,
                          to = node2seqs - 1, seqlengths = seqlengths, k = k)
          attr(node, "avdist") <- mean(dists)
        }
      }else{
        attr(node, "avdist") <- 0
      }
      return(node)
    }
    tree <- dendrapply(tree, avdist)
    reheight <- function(node){
      if(is.list(node)){
        node1edge <- max(0.001, (attr(node, "avdist") - attr(node[[1]], "avdist"))/2)
        node2edge <- max(0.001, (attr(node, "avdist") - attr(node[[2]], "avdist"))/2)
        attr(node[[1]], "height") <- attr(node, "height") - node1edge
        attr(node[[2]], "height") <- attr(node, "height") - node2edge
      }
      return(node)
    }
    reheight1 <- function(node){
      node <- reheight(node)
      if(is.list(node)) node[] <- lapply(node, reheight1)
      return(node)
    }
    tree <- reheight1(tree)
    class(tree) <- "dendrogram"
    tree <- phylogram::reposition(tree)
  }else{
    tree <- phylogram::reposition(tree)
    tree <- phylogram::ultrametricize(tree)
  }
  label <- function(node, x){
    if(is.leaf(node)) attr(node, "label") <- names(x)[attr(node, "sequences")]
    return(node)
  }
  tree <- dendrapply(tree, label, x = x)
  return(tree)
}
################################################################################

# # first temporarily remove duplicated sequences
# hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
# duplicates <- duplicated(hashes)
# hasduplicates <- any(duplicates)
# if(hasduplicates){
#   pointers <- integer(length(x))
#   dhashes <- hashes[duplicates]
#   ndhashes <- hashes[!duplicates]
#   pointers[!duplicates] <- seq_along(ndhashes)
#   pointers[duplicates] <- sapply(dhashes, match, ndhashes)
#   fullx <- x
#   x <- x[!duplicates]
# }
# # then after reheighting tree:
# if(hasduplicates){
#   add_duplicates <- function(node, pointers){
#     attr(node, "nunique") <- length(attr(node, "sequences"))
#     attr(node, "sequences") <- which(pointers %in% attr(node, "sequences"))
#     attr(node, "ntotal") <- length(attr(node, "sequences"))
#     return(node)
#   }
#   #if(!quiet) cat("Repatriating duplicate sequences with tree\n")
#   tree <- dendrapply(tree, add_duplicates, pointers)
#   forknode2 <- function(node){
#     if(!is.list(node)){# fork leaves only
#       lns <- length(attr(node, "sequences"))
#       if(lns == 1) return(node)
#       membership <- c(1, rep(2, lns - 1))
#       tmpattr <- attributes(node)
#       node <- vector(mode = "list", length = 2)
#       attributes(node) <- tmpattr
#       attr(node, "leaf") <- NULL
#       for(i in 1:2){
#         node[[i]] <- 1
#         attr(node[[i]], "height") <- attr(node, "height") - 0.0001
#         attr(node[[i]], "leaf") <- TRUE
#         attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
#       }
#     }
#     return(node)
#   }
#   forknode1 <- function(node){
#     node <- forknode2(node)
#     if(is.list(node)) node[] <- lapply(node, forknode1)
#     return(node)
#   }
#   tree <- forknode1(tree)
#   tree <- phylogram::remidpoint(tree)
#   x <- fullx
# }

#
# topdown2 <- function(node, d){
#   if(!is.list(node)){ # fork leaves only
#     seqs <- d[attr(node, "sequences"), , drop = FALSE]
#     nseq <- nrow(seqs)
#     allident <- FALSE
#     if(nseq == 1) return(node)
#     if(nseq == 2){
#       membership <- 1:2
#       if(identical(d[1,], d[2,])) allident <- TRUE
#     }else{
#       if(all(apply(seqs, 2, function(v) all(v[2:nseq] == v[1])))){
#         allident <- TRUE
#         membership <- c(1, rep(2, nseq - 1))
#       }else{
#         #test <<- seqs
#         membership <- kmeans(seqs, centers = 2)$cluster
#       }
#     }
#     tmpattr <- attributes(node)
#     node <- vector(mode = "list", length = 2)
#     attributes(node) <- tmpattr
#     attr(node, "leaf") <- NULL
#     for(i in 1:2){
#       node[[i]] <- 1
#       attr(node[[i]], "height") <- attr(node, "height") - if(allident) 0.0001 else 1
#       attr(node[[i]], "leaf") <- TRUE
#       attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
#     }
#   }
#   return(node)
# }
#
#
#
# #
# rdname topdown
# ################################################################################
# topdown.default <- function(x){
#   # x is a matrix with one row for each sequence
#   # each column is a dimension eg count for tuple k, dist to seed r, etc
#   # produces a dendogram via recursive k-means clustering
#   dupes <- duplicated.array(x, MARGIN = 1)
#   hasdupes <- any(dupes)
#   if(hasdupes){
#     whichdupes <- which(dupes)
#     lwd <- length(whichdupes)
#     tmp <- rbind(x[dupes, ], x)
#     dupes2 <- duplicated.array(tmp, MARGIN = 1)
#     whichdupes2 <- which(dupes2)
#     whichdupes2 <- whichdupes2[whichdupes2 > lwd] - lwd
#     m <- as.matrix(dist(x[whichdupes2, ]))
#     test <- as.dendrogram(hclust(dupdist, "average"))
#     test2 <- which(dupdist == 0)
#     zerodists <- which(m == 0 & lower.tri(m), arr.ind = T)
#   }
#   tree <- 1
#   attr(tree, "leaf") <- TRUE
#   attr(tree, "sequences") <- 1:nrow(x)
#   attr(tree, "height") <- 0
#   tree <- topdown1(tree, d = x)
#   tree <- phylogram::remidpoint(tree)
#   class(tree) <- "dendrogram"
#   tree <- phylogram::reposition(tree)
#   label <- function(node, d){
#     if(is.leaf(node)) attr(node, "label") <- rownames(d)[attr(node, "sequences")]
#     return(node)
#   }
#   tree <- dendrapply(tree, label, d = x)
#   return(tree)
# }
