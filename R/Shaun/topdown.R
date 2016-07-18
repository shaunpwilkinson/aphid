#' Top-down reconstruction of phylogenetic trees.
#' 
#' \code{topdown} recursively partitions a multiple sequence alignment to create
#' a dichotomous tree object of class \code{"dendrogram"}. Nodes are evaluated with
#' Akaike weights. 
#' 
#' @param x a character matrix or object of class \code{"DNAbin"} representing
#' a multiple sequence alignment.
#' @param root an integer indicating which row of the input matrix is the root sequence, 
#' or "last" (the default setting, indicating that the last row of x is the 
#' root sequence. 
#' 
#' @return an object of class \code{"dendrogram"}. 
#' 
#' @examples 
#' set.seed(99999)
#' myseq <- randomSequence(100, alphabet = c(0,1))
#' mynewick <- "((((A:3,(B:1,C:1):2):2,D:5):5,(E:1,F:1):9):11,G:20);"
#' mydendrogram <- read.dendrogram(text = mynewick) 
#' mysim <- treesim(myseq, K = 0.01, n = 25, d = mydendrogram, alphabet = c(0,1))
#' mysimx <- treesim.as.matrix(mysim)
#' reconstructed.dendro <- topdown(mysimx, root.index = "last", exhaustive = TRUE)
#' par(mfrow = c(1,2))
#' plot(mydendrogram)
#' plot(reconstructed.dendro)
#' 
enumerateSplits <- function(n){
  if(n > 15) stop("n is too large, can't compute matrix")
  res <- matrix(0, nrow = n, ncol = 2^(n-1) - 1)
  for(i in 1:ncol(res)){
    col.i <- dec2bin(i)
    col.i <- append(col.i, rep(0, nrow(res) - length(col.i)))
    res[,i] <- col.i
  }
  return(res)
}





approximateSplits <- function(x, switch.to.exhaustive = 7){
  n <- nrow(x)
  if(n <= switch.to.exhaustive){
    return(enumerateSplits(n))
  } else {
    x2 <- x[, apply(x, 2, function(v) length(unique(v)) == 2)]    
    bindices <- t(apply(x2, 1, function(y) y != x2[n,]))
    bindices.short <- unique(bindices, MARGIN = 2)
    dindices <- apply(bindices.short, 2, bin2dec)
    mergeOpp <- function(x, y, poss.splits){
      if((x <= poss.splits/2) & (y > poss.splits/2) & (poss.splits - y == x)) TRUE
      else FALSE
    }
    tof <- outer(dindices, dindices, Vectorize(mergeOpp), 2^(n - 1) - 1)
    tof <- apply(tof, 1, any)
    if(any(tof)) {
      return(bindices.short[, tof, drop = FALSE])
    } else {
      print("using all binary indices")
      return(bindices.short)
    } 
  }
}


  
nodeWeight <- function(x, M, p = 0.01, model = "JC69"){ 
  # x is a > three row alignment matrix, with root(s) at bottom
  # M is a binary matrix of candidate splits (columns) and split membership (rows)
  # can be generated from enumerateSplits() for exhaustive search, or 
  # approximateSplits() for fast heuristic searching 
  # nrow(M) should be less than nrow(x) since it doesnt include the root sequence(s)
  # p is the intrinsic mutation rate of the system and is set arbitrarily
  if(inherits(x, "DNAbin")) x <- as.character(x)
  #if(!(inherits(x, "DNAbin"))) x <- ape::as.DNAbin(x)
  s <- ncol(x) 
  taxa <- nrow(x)
  #bases.to.keep <- .C(GlobalDeletionDNA, x, taxa, s, rep(1L, s))[[4]]
  bases.to.keep <- apply(as.character(x), 2, function(v) all(v != "n"))
  
  print(sum(bases.to.keep))
  
  x <- x[, as.logical(bases.to.keep)]
  n <- ncol(x)
  nonroot.taxa <-  nrow(M)
  root.taxa <- taxa - nonroot.taxa
  root.split <- c(rep(1, nonroot.taxa), rep(0, root.taxa))
  #dis <- dist.dna(x, model = model, as.matrix = TRUE, pairwise.deletion = TRUE)
  dis <- JC69dist(x)
  dis <- round(dis *  n)
  dis[dis == 0] <- 0.0001
  nonroot.dis <- dis[1:nonroot.taxa, 1:nonroot.taxa]
  aic <- function(v, X, p = 0.01){ 
    # v is a binary split vector
    # X is a distance matrix where nrow(X) = ncol(X) = length(v)
    # outputs AIC value
    dis <- X
    tof <- outer(v, v, Vectorize(function(x, y) x != y))
    xhat <- round(mean(dis[tof]))
    dis[tof] <- xhat
    dis[dis > xhat] <- xhat
    N <- dis/p
    fun <- function(x, size) log(dbinom(x = x, size = size, prob = p))
    llk <- sum(mapply(fun, as.dist(X), as.dist(N)))
    par <- 1 + sum(as.dist(!tof))
    aic <- (2 * par) - (2 * llk)
    attr(aic, 'log.likelihood') <- llk
    attr(aic, 'no.parameters') <- par
    attr(aic, 'xhat') <- xhat
    aic
  }
  aics <- apply(M, 2, aic, nonroot.dis, p = p)
  best.split <- which(aics == min(aics))[1]
  ind <- bin2dec(M[,best.split]) # best split decimal index  
  tof2 <- outer(root.split, root.split, Vectorize(function(x, y) x != y))
  dis2 <- dis
  xhat2 <- round(mean(dis[tof2]))
  dis2[tof2] <- xhat2
  dis2[dis2 > xhat2] <- xhat2
  fun <- function(x, size) log(dbinom(x = x, size = size, prob = p))
  tof <- outer(M[,best.split], M[,best.split], Vectorize(function(x, y) x != y))
  # null model (one node)
  tof0 <- tof2
  dis0 <- dis
  tof0[1:nonroot.taxa, 1:nonroot.taxa] <- tof
  xhat0 <- round(mean(dis0[tof0]))
  dis0[tof0] <- xhat0
  dis0[dis0 > xhat0] <- xhat0
  n0 <- dis0/p
  llk0 <- sum(mapply(fun, as.dist(dis), as.dist(n0)))
  par0 <- 1 + sum(as.dist(!tof0))
  aic0 <- (2 * par0) - (2 * llk0)
  # alternative model (2 nodes)
  tof1 <- matrix(FALSE, nrow = taxa, ncol = taxa)
  nrt <- tof1
  nrt[1:nonroot.taxa, 1:nonroot.taxa] <- TRUE
  tof1[1:nonroot.taxa, 1:nonroot.taxa] <- tof # will be new candidate node
  dis1 <- dis2
  xhat1 <- round(mean(dis1[tof1]))
  dis1[tof1] <- xhat1
  dis1[nrt & dis1 > xhat1] <- xhat1
  n1 <- dis1/p
  llk1 <- sum(mapply(fun, as.dist(dis), as.dist(n1)))
  par1 <- 2 + sum(as.dist(!tof1 & !tof2))
  aic1 <- (2 * par1) - (2 * llk1)
  aics2 <- c(aic0, aic1)
  AW.numerators <- exp(-0.5 * (aics2 - min(aics2)))
  AW.denominator <- sum(AW.numerators)
  wis <- (AW.numerators/AW.denominator) 
  res <- 100 * round(wis[2], 2)
  attr(res, 'split.binary') <- M[,best.split]
  attr(res, 'split.decimal') <- ind
  attr(res, 'no.generations.to.rootnode') <- ceiling((xhat2 - xhat1)/(2 * n * p))
  attr(res, 'no.generations.to.leaf') <- ceiling(xhat1/(2 * n * p))
  res
}


axe <- function(x,  switch.to.exhaustive = 7, model = "JC69", root.index = NULL){
  if(inherits(x, "DNAbin")) x <- as.character(x)
  if(inherits(x, "matrix"){
    n <- nrow(x)
    root.split <- NULL
    if(is.null(root.index)){      
      x2 <- x[, apply(x, 2, function(v) length(unique(v)) == 2)]    
      bindices <- t(apply(x2, 1, function(y) y != x2[n,]))
      bindices <- unique(bindices, MARGIN = 2)
      root.split <- nodeWeight(x, M = bindices, model = model)
      root.index <- attr(root.split, 'split.binary')
    }
    if(identical(root.index, "last")) root.index <- n
    if(!(is.logical(root.index))) root.index <- 1:n %in% root.index
    tmp <- list(rownames(x[!root.index, , drop = FALSE]), rownames(x[root.index,
                                                                     , drop = FALSE]))
    attr(tmp[[1]], 'alig') <- x[!root.index, , drop = FALSE]
    attr(tmp[[2]], 'alig') <- x[root.index, , drop = FALSE]
    if(sum(root.index) == 1){
      if(is.null(root.split)){ 
        root.split <- nodeWeight(x, M = as.matrix(!root.index, ncol = 1), model = model)
      }
      attr(tmp[[2]], 'edge') <- attr(root.split, 'no.generations.to.leaf')
    }
    x <- tmp
  }
  if(is.list(x)){
    if(length(x[[1]]) > length(x[[2]])) x[] <- x[2:1]
    is.splittable <- sapply(x, function(y) length(y) > 1)
    if(is.splittable[1]){
      is.nonroot.left <- c(rep(1, length(x[[1]])), rep(0, length(x[[2]])))
      no.nonroot.left <- sum(is.nonroot.left)
      combined.left <- rbind(attr(x[[1]], 'alig'), attr(x[[2]], 'alig'))
      M.left <- approximateSplits(combined.left[1:(no.nonroot.left + 1),],
                                  switch.to.exhaustive = switch.to.exhaustive
                                  )[1:no.nonroot.left, , drop = FALSE]
      fit.left <- nodeWeight(combined.left, M.left, model = model)
      alig.left.left <- attr(x[[1]], 'alig')[which(attr(fit.left, 'split.binary') == 1),
                                             , drop = FALSE]
      alig.left.right <- attr(x[[1]], 'alig')[which(attr(fit.left, 'split.binary') == 0),
                                              , drop = FALSE]
      left.replacement <- list(rownames(alig.left.left), rownames(alig.left.right))
      attr(left.replacement[[1]], 'alig') <- alig.left.left
      if(nrow(alig.left.left) == 1) {
        attr(left.replacement[[1]], 'edge') <- attr(fit.left, 'no.generations.to.leaf')
      }
      attr(left.replacement[[2]], 'alig') <- alig.left.right
      if(nrow(alig.left.right) == 1) {
        attr(left.replacement[[2]], 'edge') <- attr(fit.left, 'no.generations.to.leaf')
      }
      attr(left.replacement, 'edge') <- attr(fit.left, 'no.generations.to.rootnode')
      attr(left.replacement, 'edgetext') <- fit.left[1]
    }
    if(is.splittable[2]){
      is.nonroot.right <- c(rep(1, length(x[[2]])), rep(0, length(x[[1]])))
      no.nonroot.right <- sum(is.nonroot.right)
      combined.right <- rbind(attr(x[[2]], 'alig'), attr(x[[1]], 'alig'))
      M.right <- approximateSplits(combined.right[1:(no.nonroot.right + 1),],
                                  switch.to.exhaustive = switch.to.exhaustive
                                  )[1:no.nonroot.right, , drop = FALSE]
      fit.right <- nodeWeight(combined.right, M.right, model = model)
      alig.right.left <- attr(x[[2]], 'alig')[which(attr(fit.right, 'split.binary') == 1),
                                              , drop = FALSE]
      alig.right.right <- attr(x[[2]], 'alig')[which(attr(fit.right, 'split.binary') == 0),
                                               , drop = FALSE]
      right.replacement <- list(rownames(alig.right.left), rownames(alig.right.right))
      attr(right.replacement[[1]], 'alig') <- alig.right.left
      if(nrow(alig.right.left) == 1) {
        attr(right.replacement[[1]], 'edge') <- attr(fit.right, 'no.generations.to.leaf')
      }
      attr(right.replacement[[2]], 'alig') <- alig.right.right
      if(nrow(alig.right.right) == 1) {
        attr(right.replacement[[2]], 'edge') <- attr(fit.right, 'no.generations.to.leaf')
      }
      attr(right.replacement, 'edge') <- attr(fit.right, 'no.generations.to.rootnode')
      attr(right.replacement, 'edgetext') <- fit.right[1]
    }
    if(is.splittable[1]){
      x[[1]] <- left.replacement
    }
    if(is.splittable[2]){
      x[[2]] <- right.replacement
    }
    x[] <- lapply(x, axe)
  }
  x
}

topdown <- function(x, root.index = "last", switch.to.exhaustive = 7, model = "JC69"){
  #if(!(inherits(x, "DNAbin"))) x <- ape::as.DNAbin(x)
  if(inherits(x, "DNAbin")) x <- as.character(x)
  tree.as.list <- axe(x, root.index = root.index, 
                      switch.to.exhaustive = switch.to.exhaustive, model = model)
  tree.as.dendro <- list2dendrogram(tree.as.list)
  attr(tree.as.dendro, 'class') <- 'dendrogram'
  min.height <- min(unlist(dendrapply(tree.as.dendro, attr, 'height')))
  tree.as.dendro <- dendrapply(tree.as.dendro, function(y){
    attr(y, 'height') <- attr(y, 'height') - min.height
    y
  })
  return(tree.as.dendro)
}
