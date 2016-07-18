# distanceMatrix <- function(x){
#   # x is a seq alignment
#   n <- nrow(x)
#   fun <- function(row.i, row.j){
#     pairwiseDistance(x[row.i,], x[row.j,])
#   }
#   fun <- Vectorize(fun)
#   outer(1:n, 1:n, fun)
# }

alignemtLength <- function(x){# a DNAbin object. removes gaps and returns length
  cha <- as.character(x)
  tof <- !(cha == "n" | cha == "-")
  sum(apply(tof, 2, all))
}

# Recursively compute a Hadamard matrix 
hadamard<-function(m){
  if(m != round(m)) stop ("input argument must be an integer")
  if(m < 1) stop ("input argument must be an integer greater than or equal to 1")
  m1 <- matrix(1)
  if(m==1){
    return(m1)
  } else{
    for(i in 1:m){
      firstquad <- m1
      secondquad <- firstquad
      firsthalf <- cbind(firstquad, secondquad)
      secondhalf <- cbind(firstquad, (-1* secondquad))
      m1 <- rbind(firsthalf, secondhalf)
    }
  }  
  return(m1)
}

JC69 <- function(x, n){
  d <- x/n
  if(!(d < 0.75)) stop("fraction of differing sites should not exceed 3/4")
  K <- (-3/4) * log(1 - ((4 * d)/ 3))
  #VarK <- d * (1 - d)/(n * ((1 - ((4/3) * d))^2))
  #list(K = K, VarK = VarK, stdevK = sqrt(VarK))
  K
}

JC69dist <- function(x){
  s <- ncol(x)
  n <- nrow(x)
  outer(1:n, 1:n, Vectorize(function(y, z){
    (-3/4) * log(1 - ((4 * (sum(x[y,] != x[z,])/s))/ 3))
  }))
}


# JC69bin <- function(x, n){
#   d <- x/n
#   if(!(all(d < 0.5))) stop("fraction of differing sites should not exceed 1/2")
#   K <- (-1/2) * log(1 - (2 * d))
#   as.integer(round(n * K))
# }

JC69bin <- function(x, as.matrix = FALSE){ # x is a sequence alignment, returns dist mat
  n <- ncol(x)
  taxa <- nrow(x)
  tmp <- outer(1:taxa, 1:taxa, Vectorize(function(y, z) sum(x[y,] != x[z,])))
  fun <- function(y){
    d <- y/n
    K <- (-1/2) * log(1 - (2 * d))
    K
  }
  res <- apply(tmp, c(1,2), fun)
  dimnames(res) <- list(rownames(x), rownames(x))
  if(as.matrix){
    return(res)
  } else {
    return(as.dist(res))
  }  
}

joinTaxa <- function(ind, dimn){
  # ind is the lower tri index
  # dim is the dimension of a square matrix
  stopifnot (ind <= (dimn^2 - dimn)/2)
  tfmat <- lower.tri(matrix(NA, nrow = dimn, ncol = dimn))
  actual.ind <- which(tfmat)[ind]
  i <- actual.ind %% dimn
  j <- 1 + (actual.ind %/% dimn)
  if(i == 0){ 
    i <- dimn
    j <- j - 1
  }
  return(c(i,j))
}


# how much memory does each object in global environment use?
memoryUse <- function() sort(sapply(ls(),function(x){object.size(get(x))}))


nodeBuster <- function(x){ # x is a list
  res <- x
  isnode <- sapply(res, function(y) length(y) == 2)
  extendBranch <- function(z, extension){ 
    # z is a list of length 1 with branchlength attribute (leafnode)
    attr(z, "branchlength") <- attr(z, "branchlength") + extension
    z
  }
  while(any(isnode)){
    for (i in which(isnode)){
      appendage_i <- lapply(res[[i]], extendBranch, extension = attr(res[[i]], "branchlength"))
      res <- append(res, appendage_i)
    }
    res <- res[-which(isnode)]
    isnode <- sapply(res, function(y) length(y) == 2)
  }
  attr(res, "branchlength") <- 0
  res
}

# pairwiseDistance <- function(x, y){
#   #x and y are sequence vectors of the same length
#   observed.diffs <- sum(x != y)
#   sequence.length <- length(x)
#   JC69(observed.diffs, sequence.length)[[1]]
# }

purgeGaps <- function(x){ # a DNAbin object
  cha <- as.character(x)
  tof <- !(cha == "n" | cha == "-")
  tof <- apply(tof, 2, all)
  x[, tof]
}


qMatrix <- function(x){
  # x is a distance matrix (not as.dist)
  n <- nrow(x)
  row.totals <- apply(x, 1, sum)
  col.totals <- apply(x, 2, sum)
  fun <- function(i, j){
    (n - 2) * x[i, j] - row.totals[i] - col.totals[j]
  }
  fun <- Vectorize(fun)
  outer(1:n, 1:n, fun)
}

relationship <- function (x, y, format = "binary"){
  if (format == "decimal"){
    x <- dec2bin(x)
    y <- dec2bin(y)
  }
  # append zeros to shorter vector
  lengthdiff <- length(x) - length(y)
  addon <- rep(0, abs(lengthdiff))
  if (lengthdiff < 0){
    x <- append(x, addon)
  } else if (lengthdiff > 0){
    y <- append(y, addon)
  }
  ln <- length(x)
  if (!(x[ln] == 0 & y[ln] == 0)){
    x <- append(x, 0)
    y <- append(y, 0)
  } 
  if (all(x >= y)){
    res <- "p" # x is superset of y
  } else if (all(x <= y)){
    res <- "b" # x is subset of y    
  } else if (any(x & y)){
    res <- "i" # x is incompatible with y
  } else if (sum(x, y) == length(x) - 1){
    res <- "o" # x is opposite to y
  } else {
    res <- "c" # x is compatible with y
  }
  res
}
 relationship(c(1,0,0,1), c(0,1,1,0))






