#' Pairwise distance calculation.
#'
ktup <- function(x, y, k = 4, correct = TRUE){
  N1 <- length(x)
  N2 <- length(y)
  n2q <- function(v) switch(v, 'a' = 0, '88' = 0, 'c' = 1, '28' = 1,
                            'g'= 2, '48' = 2, 't' = 3, '18' = 3,
                            sample(0:3, 1)) ### this is a hack
  xqu <- unname(sapply(x, n2q)) #x expressed as a quaternary vector
  yqu <- unname(sapply(y, n2q)) #y expressed as a quaternary vector
  xkm <- matrix(nrow = k, ncol = N1 - k + 1) # kmers for x
  ykm <- matrix(nrow = k, ncol = N2 - k + 1) # kmers for y
  for(i in 1:k){
    xkm[i, ] <- xqu[i:(N1 - (k - i))]
    ykm[i, ] <- yqu[i:(N2 - (k - i))]
  }
  xdc <- apply(xkm, 2, decimal, 4) + 1
  xtf <- 1:k^4 %in% xdc
  ydc <- apply(ykm, 2, decimal, 4) + 1
  ytf <- 1:k^4 %in% ydc
  dis <- 1 - sum(xtf & ytf)/min(sum(xtf), sum(ytf))
  if(correct){
    nmin <- min(N1, N2)
    nmax <- max(N1, N2)
    fxy <- nmin/nmax * 0.1 + 10000 /(nmax + 10000) + 0.01 # see MAFFT docs:
    # http://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html
    dis <- dis/fxy
  }
  return(dis)
}



JC69 <- function(x, fivecharstates = TRUE){
  if(!fivecharstates) x <- x[, x[1,] != "-" & x[2,] != "-"]
  d <- sum(x[1,] != x[2,])/ncol(x)
  if(fivecharstates){
    if(!(d < 0.8)) stop("fraction of differing sites should not exceed 4/5")
    K <- (-4/5) * log(1 - ((5 * d)/ 4))
    VarK <- d * (1 - d)/(n * ((1 - ((5/4) * d))^2))
  }else{
    if(!(d < 0.75)) stop("fraction of differing sites should not exceed 3/4")
    K <- (-3/4) * log(1 - ((4 * d)/ 3))
    VarK <- d * (1 - d)/(n * ((1 - ((4/3) * d))^2))
  }
  list(K = K, VarK = VarK)
}
