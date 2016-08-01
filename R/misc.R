#-------------------------------------------------------------------------------

### geometric functions
arc <- function(x, y, radx, rady = radx, from = 0, to = 2 * pi,
                no.points = 100, fill = NULL,
                arrow = NULL, arrowsize = 0.08, code = 2, ...){
  piseq <- seq(from, to, by = (to - from)/no.points)
  coords <- matrix(nrow = no.points + 1, ncol = 2)
  coords[, 1] <- x + radx * sin(piseq)
  coords[, 2] <- y + rady * cos(piseq)
  if(from == 0 & to == 2 * pi){
    polygon(coords, col = fill, ... = ...)
  }else{
    lines(coords, ... = ...)
  }
  if(!(is.null(arrow))){
    if(arrow < 0 | arrow > 1) stop ("arrow argument must be between 0 and 1")
    if(arrow == 0){
      arrows(x0 = coords[1, 1], y0 = coords[1, 2], x1 = coords[2, 1],
             y1 = coords[2, 2], code = 1, length = arrowsize, ... = ...)
    }else{
      l <- ceiling(arrow * no.points)
      arrows(x0 = coords[l, 1], y0 = coords[l, 2], x1 = coords[l + 1, 1],
             y1 = coords[l + 1, 2], code = code, length = arrowsize, ... = ...)
    }
  }
}

chord <- function(x, y, rad, type = 'outer', no.points = 100,
                  arrow = TRUE, arrowlength = 0.08, reversearrow = FALSE,
                  ...){#x and y are vectors of from, to
  distance <- sqrt(diff(x)^2 + diff(y)^2)
  tmp <- sqrt(rad^2 - (distance/2)^2) #distance from midpoint to centre
  midpoint <- c(mean(x), mean(y)) #coordinates
  angle <- atan2(diff(x), diff(y))
  opp <- sin(angle) * tmp
  adj <- cos(angle) * tmp
  centcoords <- midpoint + (c(adj, -opp))
  segangle <- 2 * acos(tmp/rad)
  startang <- atan2(x[1]- centcoords[1], y[1] - centcoords[2])
  if(startang < 0) startang <- -startang + (2 * (pi + startang))
  if(type == 'outer'){
    piseq <- seq(startang, startang + segangle, by = segangle/no.points)
  } else if (type == 'inner'){
    piseq <- seq(startang, startang - (2 * pi - segangle),
                 by = -(2 * pi - segangle)/no.points)
  } else{
    stop("type argument must be set to either 'outer' or 'inner'")
  }
  coords <- matrix(nrow = no.points + 1, ncol = 2)
  coords[, 1] <- centcoords[1] + rad * sin(piseq)
  coords[, 2] <- centcoords[2] + rad * cos(piseq)
  lines(coords, ... = ...)
  if(arrow){
    arrowdir <- if(reversearrow) -1 else 1
    arrows(x0 = coords[no.points/2 , 1],
           y0 = coords[no.points/2, 2],
           x1 = coords[no.points/2 + arrowdir, 1],
           y1 = coords[no.points/2 + arrowdir, 2],
           length = arrowlength,
           ... = ...)
  }
}

ellipse <- function(x, y, radx, rady = radx,
                    no.points = 100, col = NULL, lwd = 1){
  piseq <- seq(0, 2 * pi, by = (2 * pi)/no.points)
  coords <- matrix(nrow = no.points + 1, ncol = 2)
  coords[, 1] <- x + radx * sin(piseq)
  coords[, 2] <- y - rady * cos(piseq)
  polygon(coords, col = col, lwd = lwd)
}

diamond <- function(x, y, radx, rady = radx, col = NULL, lwd = 1){
  xcoords <- c(x + c(0, radx, 0, -radx))
  ycoords <- c(y + c(rady, 0, -rady, 0))
  polygon(xcoords, ycoords, col = col, lwd = lwd)
}

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


convert <- function(x, from = 10, to = 2, collapse = FALSE){
  if(to %% 1 > 0) stop("Non-integers are not supported yet")
  if(length(x) == 1) x <- as.integer(strsplit(paste(x), split = "")[[1]])
  x <- sum(x * from^rev(seq_along(x) - 1))
  if(to == 10){
    if(collapse) x else as.integer(strsplit(paste(x), split = "")[[1]])
  }
  dividend <- as.integer(to)
  quotient <- floor(x/dividend)
  result <- x %% dividend
  while(quotient > 0){
    remainder <- quotient %% dividend
    result = c(remainder, result)
    quotient <- floor(quotient/dividend)
  }
  if(collapse) result <- as.integer(paste(result, collapse = ""))
  return(result)
}

# x is an integer vector in 'from' numbering system (eg for binary, from = 2)
# decimal(x, from) is the same as convert(x, from, to = 10, collapse = FALSE)
decimal <- function(x, from) sum(x * from^rev(seq_along(x) - 1))



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

quickdist <- function(x, y, k = 4){
  N1 <- length(x)
  N2 <- length(y)
  n2t <- function(x) switch(x, 'a' = 0, '88' = 0, 'c' = 1, '28' = 1,
                            'g'= 2, '48' = 2, 't' = 3, '18' = 3,
                            sample(0:3, 1)) ### this is a quick hack
  xtr <- unname(sapply(x, n2t)) #x expressed as a ternary vector
  ytr <- unname(sapply(y, n2t)) #y expressed as a ternary vector
  xnc <- N1 - k + 1 # number of columns for the x matrix
  ync <- N2 - k + 1 # number of columns for the y matrix
  xsq <- 1:xnc # seq along x
  ysq <- 1:ync # seq along y
  xdf <- xtr[xsq]
  ydf <- ytr[ysq]
  for(i in 2:k){
    xsq <- xsq + 1
    ysq <- ysq + 1
    xdf <- rbind(xdf, xtr[xsq])
    ydf <- rbind(ydf, ytr[ysq])
  }
  tmp <- apply(ydf, 2, convert, from = 4, to = 10) + 1
  pointer <- apply(matrix(1:k^4, nrow = 1), 2, function(x) which(x == tmp))
  S1 <- apply(xdf, 2, function(x) pointer[[convert(x, from = 4, to = 10) + 1]])
  diags <- table(factor(unlist(mapply("-", S1, 1:xnc)), levels = -xnc:ync))
  ndgs <- length(diags)
  # need a vector of expected freqs of 'chance tuples' same length as diags vec
  # this is to correct for white noise
  # first generate vector of diag lengths
  vdglns <- k:N1 #vertical diags including 0
  hdglns <- (N2 - 1):k #vertical diags including 0
  dif <- N1 - N2
  if(dif > 0) vdglns[xnc:(xnc - dif)] <- N2
  if(dif < 0) hdglns[1:abs(dif)] <- N1
  dglns <- c(vdglns, hdglns)
  tuplns <- c(0, dglns - k + 1, 0)
  names(tuplns) <- names(diags)
  cutoffs <- sapply(tuplns, qbinom, p = 1 - 0.01/ndgs, prob = 1/length(pointer))
  expecteds <- sapply(tuplns, function(z) round(z/length(pointer)))
  signals <- diags - cutoffs >= 0
  # now subtract expected noise value from the diags with signif signal
  diags[!signals] <- 0
  diags[signals] <- diags[signals] - expecteds[signals]
  mjdgln <- tuplns[names(which.max(diags))] #major diag length
  diss <- 1 - sum(diags)/mjdgln
  diss
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
  decs <- apply(tuples, 2, convert, 2, 10)
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

whichismax <- function(v){
  ind <- which(v == max(v))
  if(length(ind) > 1) ind <- sample(ind, 1)
  ind
}
