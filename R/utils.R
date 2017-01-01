#' Detect residue alphabet.
#'
#' \code{"alphadetect"} performs checks on the format of the "residues" argument
#' to be passed to a variety of other functions. Single-element string arguments
#' such as "DNA" and "AA" are
#' converted to their respective alphabet in a character vector format.
#' @param sequences a character matrix or vector, or a list of character matrices
#' and/or vectors
#'
alphadetect <- function(sequences, residues = NULL, gapchar = "-", endchar = "?"){
  if(identical(toupper(residues), "RNA")){
    residues <- c("A", "C", "G", "U")
  }else if(is.DNA(sequences) | identical(toupper(residues), "DNA")){
    residues <- c("A", "C", "G", "T")
  }else if(is.AA(sequences) | identical(residues, "AA") | identical (toupper(residues), "AMINO")){
    residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  }
  else if(is.null(residues)){
    residues <- sort(unique(as.vector(unlist(sequences))))
    if(!is.null(gapchar)) residues <- residues[residues != gapchar]
    if(!is.null(endchar)) residues <- residues[residues != endchar]
  }else{
    if(!is.null(gapchar)) residues <- residues[residues != gapchar]
    if(!is.null(endchar)) residues <- residues[residues != endchar]
  }
  if(!(length(residues) > 0 & mode(residues) == "character")){
    stop("invalid residues argument")
  }
  return(residues)
}

#' Detect if model parameters are in log space.
logdetect <- function(x){
  if(inherits(x, "HMM")){
    if(all(x$A <= 0) & all(x$E <= 0)){
      return(TRUE)
    } else if(all(x$A >= 0) & all(x$A <= 1) & all(x$E >= 0) & all(x$E <= 1)){
      return(FALSE)
    } else stop("unable to detect if model probabilities are in log space")
  } else if(inherits(x, "PHMM")){
    if(all(x$A <= 0) & all(x$E <= 0) & all(x$qa <= 0) & all(x$qe <= 0)){
      return(TRUE)
    } else if(all(x$A >= 0) & all(x$A <= 1) & all(x$E >= 0) & all(x$E <= 1) &
              all(x$qa >= 0) & all(x$qa <= 1) & all(x$qe >= 0) & all(x$qe <= 1)){
      return(FALSE)
    } else stop("unable to detect if model probabilities are in log space")
  } else stop("x must be an object of class 'HMM' or 'PHMM'")
}

decimal <- function(x, from) sum(x * from^rev(seq_along(x) - 1))

trim <- function(x, gapchar = "-", endchar = "?", DNA = FALSE, AA = FALSE){
  #X is a raw or character matrix
  # gapchar can have length > 1
  gapchar <- if(DNA) as.raw(c(4, 240)) else if(AA) as.raw(c(45, 88)) else gapchar
  endchar <- if(DNA) as.raw(2) else if (AA) as.raw(63) else endchar
  L <- ncol(x)
  n <- nrow(x)
  if(!any(x[, 1] %in% gapchar) & !any(x[, L] %in% gapchar)) return(x)
  for(i in 1:n){
    if(x[i, 1] %in% gapchar){
      counter <- 1
      advance = TRUE
      while(advance & counter <= L){
        x[i, counter] <- endchar
        counter <- counter + 1
        advance <- if(counter <= L) x[i, counter] %in% gapchar else FALSE
      }
    }
    if(x[i, L] %in% gapchar){
      counter <- L
      advance = TRUE
      while(advance & counter >= 1){
        x[i, counter] <- endchar
        counter <- counter - 1
        advance <- if(counter >= 1) x[i, counter] %in% gapchar else FALSE
      }
    }
  }
  return(x)
}

is.DNA <- function(x){
  if(inherits(x, "DNAbin")){
    return(TRUE)
  }else if(inherits(x, "AAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                                       224, 176, 208, 112, 240, 4, 2))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x) %in% as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                             224, 176, 208, 112, 240, 4, 2))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

is.AA <- function(x){
  if(inherits(x, "AAbin")){
    return(TRUE)
  }else if(inherits(x, "DNAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(65:90, 42, 45, 63))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x) %in% as.raw(c(65:90, 42, 45, 63))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

tabulate.char <- function(x, residues, seqweights = NULL){
  if(is.null(seqweights)) seqweights <- rep(1, length(x))
  #if(identical(seqweights, 1)) seqweights <- rep(1, length(v))
  #stopifnot(length(seqweights) == length(v) & sum(seqweights) == length(v))
  res <- structure(integer(length(residues)), names = residues)
  for(i in residues) res[i] <- sum(seqweights[x == i], na.rm = TRUE)
  return(res)
}

tabulate.DNA <- function(x, ambiguities = FALSE, seqweights = NULL){
  # x is a DNAbin vector
  if(is.null(seqweights)) seqweights <- rep(1, length(x))
  #stopifnot(length(seqweights) == length(x) & sum(seqweights) == length(x))
  res <- structure(numeric(4), names = c("A", "C", "G", "T"))
  res["A"] <- sum(seqweights[x == 136])
  res["C"] <- sum(seqweights[x == 40])
  res["G"] <- sum(seqweights[x == 72])
  res["T"] <- sum(seqweights[x == 24])
  if(ambiguities){
    xambigs <- x != 4 & (x & as.raw(8)) != 8
    if(any(xambigs)){
      truncx <- x[xambigs]
      truncweights <- seqweights[xambigs]
      R <- sum(truncweights[truncx == 192])/2 # A or G
      M <- sum(truncweights[truncx == 160])/2 # A or C
      W <- sum(truncweights[truncx == 144])/2 # A or T
      S <- sum(truncweights[truncx == 96])/2 # G or C
      K <- sum(truncweights[truncx == 80])/2 # G or T
      Y <- sum(truncweights[truncx == 48])/2 # C or T
      V <- sum(truncweights[truncx == 224])/3 # A, G or C
      H <- sum(truncweights[truncx == 176])/3 # A, C or T
      d <- sum(truncweights[truncx == 208])/3 # A, G or T
      B <- sum(truncweights[truncx == 112])/3 # G, C or T
      N <- sum(truncweights[truncx == 240])/4 # A, G, C or T
      res["A"] <- res["A"] + R + M + W + V + H + d + N
      res["C"] <- res["C"] + M + S + Y + V + H + B + N
      res["G"] <- res["G"] + R + S + K + V + d + B + N
      res["T"] <- res["T"] + W + K + Y + H + d + B + N
    }
  }
  return(res)
}

tabulate.AA <- function(x, ambiguities = FALSE, seqweights = NULL){
  # x is an AAbin vector
  if(is.null(seqweights)) seqweights <- rep(1, length(x))
  tmp <- structure(numeric(26), names = LETTERS)
  guide <- as.raw(65:90)
  for(i in seq_along(tmp)) tmp[i] <- sum(seqweights[x == guide[i]])
  #tmp <- tabulate(AA2hexavigesimal(x) + 1, nbins = 26)
  res <- tmp[-c(2, 10, 15, 21, 24, 26)]
  #names(res) <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  if(ambiguities){
    #names(xambigs) <- LETTERS[c(2, 10, 15, 21, 24, 26)]
    if(sum(tmp[c(2, 10, 15, 21, 24, 26)]) > 0){
      xambigs <- x %in% guide[c(2, 10, 15, 21, 24, 26)]
      truncx <- x[xambigs]
      truncweights <- seqweights[xambigs]
      B <- sum(truncweights[truncx == guide[2]])/2 # D or N
      J <- sum(truncweights[truncx == guide[10]])/2 # L or I
      O <- sum(truncweights[truncx == guide[15]]) # Pyrrolysine
      U <- sum(truncweights[truncx == guide[21]]) # Selenocysteine
      X <- sum(truncweights[truncx == guide[24]])/20 #X
      Z <- sum(truncweights[truncx == guide[26]])/2 # E or Q
      res <- res + X
      res["D"] <- res["D"] + B
      res["N"] <- res["N"] + B
      res["L"] <- res["L"] + J
      res["I"] <- res["I"] + J
      res["K"] <- res["K"] + O
      res["C"] <- res["C"] + U
      res["E"] <- res["E"] + Z
      res["Q"] <- res["Q"] + Z
    }
  }
  return(res)
}

disambiguate.DNA <- function(a, probs = rep(0.25, 4), random = TRUE){
  # a is a raw byte in Paradis (2007) format
  # probs is a 4-element numeric vector of background probabilities for the set {a,c,g,t}
  # returns a sampled base
  if((a & as.raw(55)) == as.raw(0)){ # is purine?
    if(a != 136 & a != 72){
      if(random){
        sample(as.raw(c(136, 72)), size = 1, prob = probs[c(1, 3)]) # unknown A or G
      }else{
        as.raw(c(136, 72))[which.max(probs[c(1, 3)])]
      }
    }else{
      return(a) #known base A or G
    }
  }else if((a & as.raw(199)) == as.raw(0)){ # is pyrimidine
    if(a != 40 & a != 24){
      if(random){
        sample(as.raw(c(40, 24)), size = 1, prob = probs[c(2, 4)]) # unknown base C or T
      }else{
        as.raw(c(40, 24))[which.max(probs[c(2, 4)])]
      }
    }else{
      return(a) # known base C or T
    }
    # a,c,g,t = 136, 40, 72, 24
  }else if(a == 160){ # M (A or C)
    if(random){
      sample(as.raw(c(136, 40)), size = 1, prob = probs[1:2])
    }else{
      as.raw(c(136, 40))[which.max(probs[c(1, 2)])]
    }
  }else if(a == 144){ # W (A or T)
    if(random){
      sample(as.raw(c(136, 24)), size = 1, prob = probs[c(1, 4)])
    }else{
      as.raw(c(136, 24))[which.max(probs[c(1, 4)])]
    }
  }else if(a == 96){ # S (G or C)
    if(random){
      sample(as.raw(c(40, 72)), size = 1, prob = probs[c(2, 3)])
    }else{
      as.raw(c(40, 72))[which.max(probs[c(2, 3)])]
    }
  }else if(a == 80){ # K (G or T)
    if(random){
      sample(as.raw(c(72, 24)), size = 1, prob = probs[c(3, 4)])
    }else{
      as.raw(c(72, 24))[which.max(probs[c(3, 4)])]
    }
  }else if(a == 224){ # V (A or C or G)
    if(random){
      sample(as.raw(c(136, 40, 72)), size = 1, prob = probs[-4])
    }else{
      as.raw(c(136, 40, 72))[which.max(probs[-4])]
    }
  }else if(a == 176){
    if(random){
      sample(as.raw(c(136, 40, 24)), size = 1, prob = probs[-3]) # H (A or C or T)
    }else{
      as.raw(c(136, 40, 24))[which.max(probs[-3])]
    }
  }else if(a == 208){
    if(random){
      sample(as.raw(c(136, 72, 24)), size = 1, prob = probs[-2]) # D (A or G or T)
    }else{
      as.raw(c(136, 72, 24))[which.max(probs[-2])]
    }
  }else if(a == 112){
    if(random){
      sample(as.raw(c(40, 72, 24)), size = 1, prob = probs[-1]) # B (C or G or T)
    }else{
      as.raw(c(40, 72, 24))[which.max(probs[-1])]
    }
  }else if(a == 240){
    if(random){
      sample(as.raw(c(136, 40, 72, 24)), size = 1, prob = probs) #N
    }else{
      as.raw(c(136, 40, 72, 24))[which.max(probs)]
    }
  }else if(a == 2 | a == 4){
    return(a)
  }else stop("invalid byte for class 'DNAbin'")
}

disambiguate.AA <- function(a, probs = rep(0.05, 20), random = TRUE){
  # a is a raw byte in AAbin format
  guide <- as.raw(c(65:90, 42, 45)) #length = 28
  nonambigs <- guide[1:25][-c(2, 10, 15, 21, 24)]
   #structure(guide, class = "AAbin")
  if(a == guide[24]){
    if(random){
      sample(nonambigs, size = 1, prob = probs)
    }else{
      return(nonambigs[which.max(probs)])
    }
  }else if(a == guide[2]){# B
    if(random){
      sample(nonambigs[c(3, 12)], size = 1, prob = probs[c(3, 12)]) # D or N
    }else{
      return(nonambigs[which.max(probs[c(3, 12)])])
    }
  }else if(a == guide[10]){# J
    if(random){
      sample(nonambigs[c(8, 10)], size = 1, prob = probs[c(8, 10)]) # I or L
    }else{
      return(nonambigs[which.max(probs[c(8, 10)])])
    }
  }else if(a == guide[26]){ # Z
    if(random){
      sample(nonambigs[c(4, 14)], size = 1, prob = probs[c(4, 14)]) # E or Q
    }else{
      return(nonambigs[which.max(probs[c(4, 14)])])
    }
  }else if(a == guide[15]){# O (Pyrrolysine)
    return(nonambigs[9]) # K (Lysine)
  }else if(a == guide[21]){# U (Selenocysteine)
    return(nonambigs[2]) #C )(Cysteine)
  }else if(a == guide[27] | a == guide[28]) {
    return(NULL)
  }else stop("invalid byte for class 'AAbin'")
}

encode.DNA <- function(x, arity = 4, probs = NULL, random = TRUE, na.rm = FALSE){
  # x is a raw vector in Paradis (2007) coding scheme, possibly containing ambiguities
  # arity is an integer either 4 or 15
  if(arity == 4){
    #converts A to 0, T to 1, G to 2, and C to 3
    if(is.null(probs)) probs <- rep(0.25, 4)
    fun <- function(v){
      v <- unclass(v)
      ambigs <- (v & as.raw(8)) != 8
      if(any(ambigs)) v[ambigs] <- sapply(v[ambigs], disambiguate.DNA, probs, random)
      bits <- c(136, 24, 72, 40)
      ints <- 0:3
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(v)
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
  }else if(arity == 15){
    fun <- function(v){
      # return order A, T, G, C, S, W, R, Y, K, M, B, V, H, D, N (same as NUC4.4)
      bits <- c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240)
      ints <- 0:14
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
  }else stop("invalid 'arity' argument")
  if(is.list(x)) lapply(x, fun) else fun(x)
}

encode.AA <- function(x, arity = 20, probs = NULL, random = TRUE, na.rm = FALSE){
  # x is a vector in AAbin format, possibly containing ambiguties
  # arity is an integer, either 20, 22, 24, 26, 27, or 6 (Dayhoff6 compression)
  if(is.null(probs)) probs <- rep(0.05, 20)
  if(arity == 20){
    fun <- function(v){
      ambigs <- !(v %in% as.raw((65:89)[-c(2, 10, 15, 21, 24)]))
      if(any(ambigs)) v[ambigs] <- sapply(v[ambigs], disambiguate.AA, probs)
      bits <- (65:89)[-c(2, 10, 15, 21, 24)]
      ints <- 0:19
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
  }else if(arity == 22){
    # for use with Gonnet matrix
    # return order "C" "S" "T" "P" "A" "G" "N" "D" "E" "Q" "H" "R"
    # "K" "M" "I" "L" "V" "F" "Y" "W" "X" "*"
    # Ambig codes B, J and Z, special codes O and Z, are returned as 20 (X),
    # and gaps are returned as NA
    fun <- function(v){
      bits <- c(67, 83, 84, 80, 65, 71, 78, 68, 69, 81, 72, 82, 75, 77,
                73, 76, 86, 70, 89, 87, 88, 42, 66, 74, 79, 85, 90)
      ints <- c(0:21, rep(20, 5))
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      isnares <- is.na(res) #placeholder for ambig treatment
      if(na.rm) if(any(isnares)) res <- res[!isnares]
      return(res)
    }
  }else if(arity == 24){
    # for use with PAM and BLOSUM matrices
    # return order "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M"
    # "F" "P" "S" "T" "W" "Y" "V" "B" "Z" "X" "*"
    # Ambig code J, and special codes O and U are returned as 22 (X).
    # Gaps are returned as NA or removed if na.rm = T
    fun <- function(v){
      bits <- c(65, 82, 78, 68, 67, 81, 69, 71, 72, 73, 76, 75,
                77, 70, 80, 83, 84, 87, 89, 86, 66, 90, 88, 42, 74, 79, 85)
      ints <- c(0:23, 22, 22, 22)
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      isnares <- is.na(res)
      if(na.rm) if(any(isnares)) res <- res[!isnares]
      return(res)
    }
  }else if(arity == 26){
    ### return order = LETTERS
    fun <- function(v){
      res <- as.integer(v) - 65
      res[res < 0 | res > 25] <- NA
      attributes(res) <- attributes(unclass(v))
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
  }else if(arity == 27){
    ### return order ACDEFGHIKLMNPQRSTVWY, X, BJZ, OU, *
    ### for input into AAprobC2
    fun <- function(v){
      bits <- c(65, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 80,
                81, 82, 83, 84, 86, 87, 89, 88, 66, 74, 90,  79,85, 42)
      ints <- c(0:26)
      res <- ints[match(as.numeric(v), bits)]
      attributes(res) <- attributes(unclass(v))
      isnares <- is.na(res)
      if(na.rm) if(any(isnares)) res <- res[!isnares]
      return(res)
    }
  }else if(arity == 6){
    fun <- function(y){
      y <- unclass(y)
      bits <- 65:90
      ints <- c(0, 2, 1, 2, 2, 3, 0, 4, 5, 5, 4, 5, 5, 2, 4, 0, 2, 4, 0, 0, 1, 5, 3, NA, 3, 2)
      res <- ints[match(as.numeric(y), bits)]
      attributes(res) <- attributes(y)
      if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
      return(res)
    }
  }else stop("invalid 'arity' argument")
  if(is.list(x)) lapply(x, fun) else fun(x)
}

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

insertgaps <- function(x, positions, lengths, gapchar = "-"){
  if(length(lengths) != length(positions)){
    stop("arguments 'lengths' and 'positions' should of be equal length")
  }
  if(length(lengths) == 0) return(x)
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
  gapcol <- matrix(rep(gapchar, n), ncol = 1, dimnames = list(rownames(x), NULL))
  res <- cbind(gapcol, x, deparse.level = 0)[, indices]
  #res <- cbind(gapcol, x)[, indices]
  if(xisvec) res <- as.vector(res)
  return(res)
}
# x <- rbind(c("A","B","C","D","E","F","G","H"), c("A","B","C","D","E","F","G","H"))
# insertgaps(x, positions = c(0, 8), lengths = c(3, 5))
# x <- c("A","B","C","D","E","F","G","H")
# insertgaps(x, positions = c(2, 6), lengths = c(3, 5))
#x <- 1:10
#insertgaps(x, positions = c(2, 6), lengths = c(3, 5), gapchar = NA)


# this function is used in alignpair for matrix x matrix alignment
insert <- function(x, into, at){
  if(ncol(x) == 0) return(into)
  into[, at:(at + ncol(x) - 1)] <- x
  return(into)
}
