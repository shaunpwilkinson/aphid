################################################################################
topdown <- function(x){
  # x is a matrix with one row for each sequence
  # each column is a dimension eg count for tuple k, dist to seed r, etc
  # uses recursive k-means clustering to produce a dendogram
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "sequences") <- 1:nrow(x)
  attr(tree, "height") <- 0
  tree <- topdown1(tree, sequences = x)
  tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  tree <- phylogram::reposition(tree)
  tree <- phylogram::ultrametricize(tree)
  label <- function(node, sequences){
    if(is.leaf(node)){
      attr(node, "label") <- rownames(sequences)[attr(node, "sequences")]
    }
    return(node)
  }
  tree <- dendrapply(tree, label, sequences = x)
  return(tree)
}

topdown1 <- function(tree, sequences){
  tree <- topdown2(tree, sequences = sequences)
  if(is.list(tree)) tree[] <- lapply(tree, topdown1, sequences = sequences)
  return(tree)
}

topdown2 <- function(node, sequences){
  if(!is.list(node)){ # fork leaves only
    seqs <- sequences[attr(node, "sequences"), , drop = FALSE]
    nseq <- nrow(seqs)
    allident <- FALSE
    if(nseq == 1) return(node)
    if(nseq == 2){
      membership <- 1:2
      if(identical(sequences[1,], sequences[2,])) allident <- TRUE
    }else{
      # if(all(apply(sequences, 2, function(v) all(v[2:nseq] == v[1])))){
      #   allident <- TRUE
      #   membership <- c(1, rep(2, nseq))
      # }else{
        membership <- kmeans(seqs, centers = 2)$cluster
     # }
    }
    tmpattr <- attributes(node)
    node <- vector(mode = "list", length = 2)
    attributes(node) <- tmpattr
    attr(node, "leaf") <- NULL
    for(i in 1:2){
      node[[i]] <- 1
      attr(node[[i]], "height") <- attr(node, "height") - if(allident) 0.0001 else 1
      attr(node[[i]], "leaf") <- TRUE
      attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
    }
  }
  return(node)
}
