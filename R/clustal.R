#' Multiple sequence alignment.
#'
#' Align multiple DNA or protein
#' sequences using the clustal omega algorithm.
#'
#' @param x a named list of character vectors representing the sequences to be aligned.
#' @param quick logical value indicating whether the Wilbur Lipman k-tuples heuristic
#' should be used for the calculation of the pairwise distance matrix.
#' @param pivot logical value indicating whether pivot objects should be used
#' in mbed to allocate the most dissimilar seeds.
#' @param type a character string specifying whether the alignment should be
#' 'global' (penalized end gaps) or 'semiglobal' (default; free end gaps).
#' @return a character matrix of aligned sequences.
#'
clustalo <- function(x, quick = TRUE, pivot = FALSE, type = 'semiglobal'){
  # x is a list of sequences
  x <- x[order(lengths(x))]
  N <- length(x)
  nseeds <- if(N < 20) N else (log(N, 2)^2)
  seeds <- round(seq(from = 1, to = N, by = N/nseeds))
  seeds[length(seeds)] <- N
  ## distfun takes 2 character vectors and returns distance value
  distfun <- if(quick) kmer else function(x, y)
    JC69(WilburLipman(x, y)$globalAlignment)$K
  if(pivot){
    outliers <- integer()
    for(i in 1:length(seeds)){
      seed <- x[[seeds[i]]]
      dists1 <- sapply(x, distfun, y = seed)
      l <- x[which.max(dists1)]
      dists2 <- sapply(x, distfun, y = l)
      candidate <- which.max(dists2)
      if(!candidate %in% seeds) outliers <- c(outliers, candidate)
    }
    seeds <- c(seeds, outliers)
  }
  ## mbed takes seq char vector and returns vector of dists from each seed
  mbed <- function(z) sapply(x[seeds], distfun, y = z)
  mbdsqs <- t(as.matrix(data.frame(lapply(x, mbed))))
  guidetree <- as.dendrogram(hclust(dist(mbdsqs), method = "average"))
  newick <- write.dendrogram(guidetree, strip.edges = TRUE)
  newick <- gsub(";", "", newick)
  newick <- gsub("\\(", "alignpair\\(", newick)
  if(type == 'semiglobal') newick <- gsub("\\)", ", type = 'semiglobal'\\)", newick)
  res <- with(x, eval(parse(text = newick)))
  res
}


