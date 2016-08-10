#' Pairwise alignment of sequences and/or sequence alignments.
#'
#' \code{align} uses the Viterbi algorithm to find the optimal alignment
#' between two sequences, a sequence and an alignment, or two alignments.
#'
#' @param x,y character vectors or character matrices of aligned sequences.
#' @param d,e gap opening and gap extension penalties for pairwise sequence
#' alignment.
#' @param S an optional substitution matrix with \code{dimnames} attributes
#' corresponding to the residue alphabet. If NULL matches are
#' scored as 1 and mismatches as -1.
#' @param type a character string specifying whether the alignment should be
#' 'global' (penalized end gaps), 'semiglobal' (default; free end gaps) or
#' local (highest scoring subalignment).
#' @return a character matrix of aligned sequences.
#' @examples
#' x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#' y <- c("P", "A", "W", "H", "E", "A", "E")
#' z <- align(x, y)
#' align(x, z)
#'
align <- function(x, y, d = 8, e = 2, S = NULL, qe = NULL,
                  logspace = FALSE, type = 'global',
                  offset = -0.1, itertab = NULL, method = 'background',
                  residues = 'auto', gapchar = "-"){
  if(is.vector(x) & is.vector(y)){
    alig <- Viterbi(x, y, d = d, e = e, S = S, type = type,
                    itertab = itertab)#, offset = offset) ###not necessary for vec vs vec
    xind <- yind <- alig$path
    xind[alig$path != 3] <- 1:length(x)
    xind[alig$path == 3] <- 0
    newx <- c(gapchar, x)[xind + 1]
    yind[alig$path != 1] <- 1:length(y)
    yind[alig$path == 1] <- 0
    newy <- c(gapchar, y)[yind + 1]
    res <- rbind(newx, newy)
    rownames(res) <- c(deparse(substitute(x)), deparse(substitute(y)))
    return(res)
  }else if((is.matrix(x) & is.vector(y)) | (is.vector(x) & is.matrix(y))){ ###also need option to flip
    vm <- is.vector(x) & is.matrix(y)
    if(vm){
      tmp1 <- deparse(substitute(x))
      tmp2 <- x
      x <- y
      y <- tmp2
    }else{
      tmp1 <- deparse(substitute(y))
    }
    if(identical(residues, "auto")) residues <- sort(unique(c(as.vector(x), y)))
    n <- nrow(x)
    z <- derivePHMM(x, method = method, residues = residues)
    l <- z$size
    alig <- Viterbi(z, y, qe = qe, logspace = logspace, type = type,
                    offset = offset, itertab = itertab)
    yind <- alig$path
    yind[alig$path != 1] <- 1:length(y)
    yind[alig$path == 1] <- 0
    newy <- c(gapchar, y)[yind + 1]
    #also need to account for inserts in x
    ynotinsert <- alig$path != 3 #logical vector
    yinsertlengths <- insertlengths(ynotinsert) #tabulate insert lengths
    xinsertlengths <- z$insertlengths
    #reconcile x and y insert lengths
    # gls = gap lengths, gps = gap positions
    diffsx <- diffsy <- yinsertlengths - xinsertlengths
    newx <- x
    if(any(diffsx > 0)){
      diffsx[diffsx < 0] <- 0
      xgls <- diffsx[diffsx > 0]
      xgps <- c(z$alignment, ncol(x) + 1)[which(diffsx > 0)] - 1
      newx <- insertgaps(newx, xgps, xgls, gapchar = gapchar)
    }
    if(any(diffsy < 0)){
      diffsy[diffsy > 0] <- 0
      #align y to model
      yprog <- alig$progression
      ygls <- -1 * diffsy[diffsy < 0]
      ygps <- sapply(which(diffsy < 0), match, c(yprog, l + 1)) - 1
      newy <- insertgaps(newy, ygps, ygls, gapchar = gapchar)
    }
    res <- rbind(newx, newy)
    rownames(res)[n + 1] <- tmp1
    return(res)
  }else if(is.matrix(x) & is.matrix(y)){
    if(identical(residues, "auto")) residues <- sort(unique(c(as.vector(x), as.vector(y))))
    nx <- nrow(x)
    ny <- nrow(y)
    zx <- derivePHMM(x, method = method, residues = residues)
    zy <- derivePHMM(y, method = method, residues = residues)
    lx <- zx$size
    ly <- zy$size
    alig <- Viterbi(zx, zy, qe = qe, logspace = logspace, type = type,
            offset = offset, itertab = itertab)
    #vectors same length as model (+ 1 for begin state) with counts of
    #no of gaps to insert after each position
    xinsertlengths <- insertlengths(alig$path < 4) - zx$insertlengths
    yinsertlengths <- insertlengths(alig$path > 2) - zy$insertlengths
    #these vecs correspond to orig alig cols,
    #keeps track of how many gaps to insert after each
    resx <- rep(0, ncol(x) + 1)
    names(resx) <- 0:ncol(x)
    resy <- rep(0, ncol(y) + 1)
    names(resy) <- 0:ncol(y)
    zxali <- c(0, zx$alignment, ncol(x) + 1)
    zyali <- c(0, zy$alignment, ncol(y) + 1)
    names(zxali) <- 0:(lx + 1)
    names(zyali) <- 0:(ly + 1) # includes begin and end states
    prog <- alig$progression
    if(!(all(prog[, 1] == 0))) prog <- cbind(c(0, 0), prog)
    ypve <- which(yinsertlengths > 0)
    ynve <- which(yinsertlengths < 0)
    xpve <- which(xinsertlengths > 0)
    xnve <- which(xinsertlengths < 0)
    #ptiga = positions to insert gaps after
    #ctiga = columns (of original alignment) to insert gaps after
    if(length(ypve) > 0){
      yctiga <- zyali[ypve + 1] - 1
      ygls <- yinsertlengths[ypve]
      resy[yctiga + 1] <- ygls
    }
    if(length(ynve) > 0){
      xptiga <- prog[1, sapply(ynve - 1, match, prog[2,])] ## #may need checking
      xctiga <- zxali[xptiga + 2] - 1
      xgls <- -1 * yinsertlengths[ynve]
      resx[xctiga + 1] <- xgls
    }
    if(length(xpve) > 0){
      xctiga <- zxali[xpve + 1] - 1
      xgls <- xinsertlengths[xpve]
      resx[xctiga + 1] <- resx[xctiga + 1] + xgls
    }
    if(length(xnve) > 0){
      yptiga <- prog[2, sapply(xnve - 1, match, prog[1,])]
      yctiga <- zyali[yptiga + 2] - 1
      ygls <- -1 * xinsertlengths[xnve]
      resy[yctiga + 1] <- resy[yctiga + 1] + ygls
    }
    newx <- insertgaps(x, which(resx > 0) - 1, resx[resx > 0], gapchar = gapchar)
    newy <- insertgaps(y, which(resy > 0) - 1, resy[resy > 0], gapchar = gapchar)
    res <- rbind(newx, newy)
    res <- res[, apply(res, 2, function(v) !all(v == gapchar))]
    return(res)
  }else{
    stop("invalid arguments provided for x and or y")
  }
}
