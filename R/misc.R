#-------------------------------------------------------------------------------

### other functions

longisland <- function(x) max(which(c(F,x) & !c(x,F)) - which(!c(F,x) & c(x,F)))

logsumR <- function(x){
  n <- length(x)
  if(n == 1) return(x)
  res <- x[1]
  for(i in 2:n){
    if(res == -Inf) res <- x[i]
    else if(x[i] == -Inf) res <- res
    else if(res > x[i]) res <- res + log1p(exp(x[i] - res))
    else res <- x[i] + log1p(exp(res - x[i]))
  }
  return(unname(res))
}


addmats <- function(x){ # list of matrices all of same dimension
  dmn <- dimnames(x[[1]])
  n <- nrow(x[[1]])
  m <- ncol(x[[1]])
  tmp <- unlist(x, use.names = FALSE)
  tmp <- matrix(tmp, nrow = n * m)
  tmp <- apply(tmp, 1, sum)
  res <- matrix(tmp, nrow = n, dimnames = dmn)
  res
}


pathfinder <- function(v){
  tmp <- v[-(length(v))] - v[-1]
  transtype <- function(x){
    if(x == 0) "I" else if(x == -1) "M" else if(
      x < -1) rep("D", abs(x + 1)) else stop("expected ascending vector")
  }
  unlist(lapply(tmp, transtype))
}

insertcolumn <- function(x, what, where){
  if(nrow(x) != length(what)) stop("new column has different length to destination")
  tmp <- rep(TRUE, ncol(x) + length(where))
  tmp[where] <- FALSE
  index <- as.numeric(tmp)
  index[tmp] <- 1:ncol(x)
  index <- index + 1
  cbind(what, x, deparse.level = 0)[, index]
}

insertgaps <- function(x, positions, lengths, gapchar = "-"){
  if(length(lengths) != length(positions)){
    stop("arguments 'lengths' and 'positions' should of be equal length")
  }
  #positions: vector, which matrix columns should the gaps be inserted after?
  # (can include zero)
  #lengths: vector, how long is each gap?
  xisvec <- is.vector(x)
  if(xisvec) x <- matrix(x, nrow = 1)
  n <- nrow(x)
  m <- ncol(x)
  tmp <- rep(TRUE, ncol(x) + sum(lengths))
  tab <- rep(0, m + 1)
  names(tab) <- 0:m
  tab[positions + 1] <- lengths
  fun <- function(e) if(e == 0) TRUE else c(TRUE, rep(FALSE, e))
  notgap <- unlist(lapply(tab, fun))[-1]   #remember to delete true #1
  indices <- as.numeric(notgap)
  indices[notgap] <- 1:m
  indices <- indices + 1
  res <- cbind(rep(gapchar, n), x, deparse.level = 0)[, indices]
  if(xisvec) res <- as.vector(res)
  return(res)
}
# x <- rbind(c("A","B","C","D","E","F","G","H"), c("A","B","C","D","E","F","G","H"))
# insertgaps(x, positions = c(0, 8), lengths = c(3, 5))
# x <- c("A","B","C","D","E","F","G","H")
# insertgaps(x, positions = c(2, 6), lengths = c(3, 5))
#x <- 1:10
#insertgaps(x, positions = c(2, 6), lengths = c(3, 5), gapchar = NA)

# given logical vector, how many falses are after each true?
# note - also outputs a zero position
insertlengths <- function(x){
  tuples <- rbind(c(T, x), c(x, T))
  decs <- apply(tuples, 2, decimal, 2)
  startsl <- decs == 2
  starts <- which(startsl) #insert start positions
  ends <- which(decs == 1)
  lengths <- ends - starts
  tmp <- as.numeric(c(T, x))
  tmp[startsl] <- lengths
  tmp[!startsl] <- 0
  res <- tmp[c(T, x)]
  names(res) <- 0:(length(res) - 1)
  res
}
#x <- c(F,F,T,T,T,F,F,T,T,T,F,T,T,F,F,F,F,F,T,F,F,F)
#insertlengths(x)

insert <- function(X, v, index, margin = 2){
  #x is a matrix, v is a vector, index is the row/col no. to insert
  if(is.null(index)) return(X)
  if(is.na(index)) return(X)
  xisvec <- is.vector(X)
  if(xisvec) X <- matrix(X, nrow = 1)
  n <- nrow(X)
  m <- ncol(X)
  q <- length(index)
  if(margin == 1){
    #if(length(v) != m) stop("new row has wrong length")
    v <- recycle(v, m)
    res <- rbind(X, matrix(rep(v, q), nrow = q, byrow = T), deparse.level = 0)
    guide <- c((1:(n + q))[-index], index)
    if(any(guide > n + q)) stop("index out of bounds")
    res <- res[order(guide), ]
  }else if(margin == 2){
    #if(length(v) != n) stop("new column has wrong length")
    v <- recycle(v, n)
    res <- cbind(X, matrix(rep(v, q), ncol = q), deparse.level = 0)
    guide <- c((1:(m + q))[-index], index)
    if(any(guide > m + q)) stop("index out of bounds")
    res <- res[, order(guide)]
  }
  if(xisvec & margin == 2) res <- as.vector(res)
  res
}

recycle <- function(v, len){
  if(length(v) < len){
    quotient <- floor(len/length(v))
    remainder <- len %% length(v)
    v <- c(rep(v, quotient), if(remainder > 0) v[1:remainder] else NULL)
  }
  return(v)
}


