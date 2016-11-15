#' Utilities.
#'

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
    return(all(x %in% as.raw(c(65:90, 42, 45))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x) %in% as.raw(c(65:90, 42, 45))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

#-----------------------------------------------------------------------------

decimal <- function(x, from) sum(x * from^rev(seq_along(x) - 1))


# whichismax <- function(v){
#   ind <- which(v == max(v, na.rm = TRUE))
#   if(length(ind) > 1) ind <- sample(ind, 1)
#   ind
# }

#
# progression <- function(x){ # an object of class "Viterbi"
#   res <- matrix(nrow = 2, ncol = length(x$path))
#   xpos <- res[1, 1] <- x$start[1]
#   ypos <- res[2,1 ] <- x$start[2]
#   for(i in 1:(length(x$path) - 1)){
#     if(x$path[i] == 0) {
#       xpos <- xpos + 1
#     } else if(x$path[i] == 1){
#       xpos <- xpos + 1
#       ypos <- ypos + 1
#     } else if(x$path[i] == 2){
#       ypos <- ypos + 1
#     } else stop("path contains unknown elements")
#     res [1, i + 1] <- xpos
#     res [2, i + 1] <- ypos
#   }
#   res
# }

#' Detect residue alphabet.
#'
#' \code{"alphadetect"} performs checks on the format of the "residues" argument
#' to be passed to a variety of other functions. Single-element string arguments
#' such as "DNA" and "AA" are
#' converted to their respective alphabet in a character vector format.
#' @param sequences a character matrix or vector, or a list of character matrices
#' and/or vectors
#'
alphadetect <- function(sequences, residues = NULL, gapchar = "-"){
  if(is.DNA(sequences) | identical(residues, "DNA")){
    residues <- c("A", "C", "G", "T")
  } else if(is.AA(sequences) | identical(residues, "AA")){
    residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  }
  else if(is.null(residues)){
    residues <- sort(unique(as.vector(unlist(sequences))))
    if(!is.null(gapchar)) residues <- residues[residues != gapchar]
  }else{
    if(!is.null(gapchar)) residues <- residues[residues != gapchar]
  }
  if(!(length(residues) >= 1 & mode(residues) == "character")){
    stop("invalid residues argument")
  }
  return(residues)
}


tabulate.char <- function(x, residues, seqweights = NULL){
  if(is.null(seqweights)) seqweights <- rep(1, length(x))
  #if(identical(seqweights, 1)) seqweights <- rep(1, length(v))
  #stopifnot(length(seqweights) == length(v) & sum(seqweights) == length(v))
  res <- structure(integer(length(residues)), names = residues)
  for(i in residues) res[i] <- sum(seqweights[x == i], na.rm = TRUE)
  return(res)
}

#-------------------------------------------------------------------------------------------
# x is a DNAbin vector
tabulate.DNA <- function(x, ambiguities = FALSE, seqweights = NULL){
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


# compress.AA2 <- function(x){
#   indices <- c(rep(0, 5), rep(1, 2), rep(2, 6), rep(3, 3), rep(4, 4), rep(5, 5), 6, 6)
#   names(indices) <- unlist(strsplit("AGPSTCUDENQBZFWYHKROILMVJX-", split = ""))
#   res <- indices[as.character.AAbin(x)]
#   res <- res[res < 6]
#   return(res)
# }
#
# compress.AA3 <- function(x, alphabet = "Dayhoff6"){
#   if(!identical(alphabet, "Dayhoff6")) stop("only Dayhoff6 alphabet supported")
#   res <- integer(length(x))
#   res[x %in% as.raw(c(65, 71, 80, 83, 84))] <- 1
#   res[x %in% as.raw(c(67, 85))] <- 2
#   res[x %in% as.raw(c(68, 69, 78, 81, 66, 90))] <- 3
#   res[x %in% as.raw(c(70, 87, 89))] <- 4
#   res[x %in% as.raw(c(72, 75, 82, 79))] <- 5
#   res[x %in% as.raw(c(73, 76, 77, 86, 74))] <- 6
#   res <- res[res > 0]
#   return(res - 1)
# }

compress.AA <- function(x, alpha = "Dayhoff6", na.rm = FALSE){
  if(!identical(alpha, "Dayhoff6")) stop("only Dayhoff6 alphabet supported in this version")
  fun <- function(y){
    y <- unclass(y)
    bits <- 65:90
    ints <- c(0, 2, 1, 2, 2, 3, 0, 4, 5, 5, 4, 5, 5, 2, 4, 0, 2, 4, 0, 0, 1, 5, 3, NA, 3, 2)
    res <- ints[match(as.numeric(y), bits)]
    attributes(res) <- attributes(y)
    if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
    return(res)
  }
  if(is.list(x)) lapply(x, fun) else fun(x)
}


#' Diagnostic model checks.
#'
validate <- function(x){
  if(inherits(x, "HMM")){
    states <- rownames(x$A)
    if(is.null(colnames(x$A)) | is.null(colnames(x$E))){
      message('both A and E must have dimnames attributes')
    } else if(!(identical(states, colnames(x$A)) & identical(states[-1], rownames(x$E)))){
      message("rownames and colnames of the transitions matrix (A), and
              rownames of the emissions matrix (E) must all be non-NULL
              and identical")
    } else if(rownames(x$A)[1] != "BeginEnd") {
      message("transition probability matrix (A) should include a 'BeginEnd'
              state at row 1 and col 1")
    } else {
      return(TRUE)
    }
    return(FALSE)
  } else if(inherits(x, "PHMM")){
    states <- c("M", "I", "D")
    if(!(all(dimnames(x$A)[[1]] == states) & all(dimnames(x$A)[[3]] == states))){
      message("names for dimensions 1 and 3 of transitions
              array (A) must be c('D', 'M', 'I')")
    } else if(!(ncol(x$A) == ncol(x$E) + 1)){
      message("transitions array (A) should include a begin state")
    } else if(!all(names(x$qe) == rownames(x$E))){
      message("residue names for background emissions vector (qe) should be identical
              to rownames of emissions matrix (E)")
    } else if(!(all(colnames(x$qa) == states) & all(rownames(x$qa) == states))){
      message("rownames and colnames of background transitions matrix
              must be c('D', 'M', 'I')")
    } else {
      return(TRUE)
    }
    return(FALSE)
  } else {
    stop("x must be an object of class 'HMM' or 'PHMM'")
  }
}

#' Detect if model parameters are in log space.
#'
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

probDNA <- function(a, probs = rep(0.25, 4)){
  # a is a raw byte in Paradis (2007) format
  # probs is a 4-element numeric vector of probabilities for the set {a,c,g,t}
  # returns the weighted average probability
  if(a <= 4){
    stop("Input sequence contains gaps or unknown characters")
  }else if(length(probs) != 4){
    stop("Vector of probabilities should be of length 4")
  }else{
    if((a & as.raw(55)) == as.raw(0)){ # is purine?
      if(a == 136){
        probs[1] # A
      } else if(a == 72){
        probs[3] # G
      } else{
        mean(probs[c(1, 3)]) # A or G  ### what if probs are logged??
      }
    } else if((a & as.raw(199)) == as.raw(0)){
      if(a == 40){
        probs[2] # C
      } else if(a == 24){
        probs[4] # T
      } else{
        mean(probs[c(2, 4)]) # A or G  ### what if probs are logged??
      }
    } else if(a == 160){
      mean(probs[1:2])# M (A or C)
    }else if(a == 144){
      mean(probs[c(1, 4)])# W (A or T)
    }else if(a == 96){
      mean(probs[2:3])# S (G or C)
    }else if(a == 80){
      mean(probs[3:4])# K (G or T)
    }else if(a == 224){
      mean(probs[-4]) # V (A or C or G)
    } else if(a == 176){
      mean(probs[-3]) # H (A or C or T)
    } else if(a == 208){
      mean(probs[-2]) # D (A or G or T)
    } else if(a == 112){
      mean(probs[-1]) # B (C or G or T)
    } else if(a == 240){
      mean(probs) #N
    } else stop("invalid byte")
  }
}

is.ambiguous.DNA <- function(a) a != 4 & (a & as.raw(8)) != 8

disambiguate.DNA <- function(a, probs = rep(0.25, 4), random = TRUE){
  # a is a raw byte in Paradis (2007) format
  # probs is a 4-element numeric vector of background probabilities for the set {a,c,g,t}
  # returns a sampled base
  if((a & as.raw(55)) == as.raw(0)){ # is purine?
    if(a != 136 & a != 72){
      sample(as.raw(c(136, 72)), size = 1, prob = probs[c(1, 3)]) # unknown A or G
    }else{
      return(a) #known base A or G
    }
  }else if((a & as.raw(199)) == as.raw(0)){ # is pyrimidine
    if(a != 40 & a != 24){
      sample(as.raw(c(40, 24)), size = 1, prob = probs[c(2, 4)]) # unknown base C or T
    }else{
      return(a) # known base C or T
    }
    # a,c,g,t = 136, 40, 72, 24
  }else if(a == 160){
    sample(as.raw(c(136, 40)), size = 1, prob = probs[1:2])
    # M (A or C)
  }else if(a == 144){
    sample(as.raw(c(136, 24)), size = 1, prob = probs[c(1, 4)])
    # W (A or T)
  }else if(a == 96){
    sample(as.raw(c(40, 72)), size = 1, prob = probs[2:3])
    # S (G or C)
  }else if(a == 80){
    sample(as.raw(c(72, 24)), size = 1, prob = probs[3:4])
    # K (G or T)
  }else if(a == 224){ # V (A or C or G)
    sample(as.raw(c(136, 40, 72)), size = 1, prob = probs[-4])
  } else if(a == 176){
    sample(as.raw(c(136, 40, 24)), size = 1, prob = probs[-3]) # H (A or C or T)
  } else if(a == 208){
    sample(as.raw(c(136, 72, 24)), size = 1, prob = probs[-2]) # D (A or G or T)
  } else if(a == 112){
    sample(as.raw(c(40, 72, 24)), size = 1, prob = probs[-1]) # B (C or G or T)
  } else if(a == 240){
    sample(as.raw(c(136, 40, 72, 24)), size = 1, prob = probs) #N
  } else if(a == 2 | a == 4){
    return(a)
  }else stop("invalid byte for class 'DNAbin'")
}

disambiguate.AA <- function(a, probs = rep(0.05, 20)){
  # a is a raw byte in AAbin format
  guide <- as.raw(c(65:90, 42, 45)) #length = 28
  nonambigs <- guide[1:25][-c(2, 10, 15, 21, 24)]
   #structure(guide, class = "AAbin")
  if(a == guide[24]){
    sample(nonambigs, size = 1, prob = probs)
  }else if(a == guide[2]){# B
    sample(nonambigs[c(3, 12)], size = 1, prob = probs[c(3, 12)]) # D or N
  }else if(a == guide[10]){# J
    sample(nonambigs[c(8, 10)], size = 1, prob = probs[c(8, 10)]) # I or L
  }else if(a == guide[26]){ # Z
    sample(nonambigs[c(4, 14)], size = 1, prob = probs[c(4, 14)]) # E or Q
  }else if(a == guide[15]){# O (Pyrrolysine)
    return(nonambigs[9]) # K (Lysine)
  }else if(a == guide[21]){# U (Selenocysteine)
    return(nonambigs[2]) #C )(Cysteine)
  }else if(a == guide[26]){ # Z
    sample(nonambigs[c(4, 14)], size = 1, prob = probs[c(4, 14)]) # E or Q
  }else if(a == guide[27] | a == guide[27]) {
    return(NULL)
  }else stop("invalid byte for class 'AAbin'")
}

DNA2quaternary <- function(x, probs = rep(0.25, 4), random = TRUE, na.rm = FALSE){
  fun <- function(v){
    # v is a raw vector in Paradis (2007) scheme, possibly containing ambiguities
    #converts A to 0, T to 1, G to 2, and C to 3
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
  if(is.list(x)) lapply(x, fun) else fun(x)
}

DNA2pentadecimal <- function(x, na.rm = FALSE){
  fun <- function(v){
    # v is a raw vector in Paradis (2007) scheme, possibly containing ambiguities
    #converts A to 0, T to 1, G to 2, and C to 3
    # return order A, T, G, C, S, W, R, Y, K, M, B, V, H, D, N (same as NUC4.4)
    bits <- c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240)
    ints <- 0:14
    res <- ints[match(as.numeric(v), bits)]
    attributes(res) <- attributes(unclass(v))
    if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
    return(res)
  }
  if(is.list(x)) lapply(x, fun) else fun(x)
}

AA2hexavigesimal <- function(x, na.rm = FALSE){
  ### return order = LETTERS
  fun <- function(v){
    res <- as.integer(v) - 65
    res[res < 0 | res > 25] <- NA
    attributes(res) <- attributes(unclass(v))
    if(na.rm) if(any(is.na(res))) res <- res[!is.na(res)]
    return(res)
  }
  if(is.list(x)) lapply(x, fun) else fun(x)
}

AA2heptovigesimal <- function(x, na.rm = FALSE){
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
  if(is.list(x)) lapply(x, fun) else fun(x)
}

AA2vigesimal <- function(x, probs = rep(0.05, 20), random = TRUE, na.rm = FALSE){
  # return order "A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y"
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
  if(is.list(x)) lapply(x, fun) else fun(x)
}

AA2quadrovigesimal <- function(x, na.rm = FALSE){
  # for use with PAM and BLOSUM matrices
  # return order "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y" "V" "B" "Z" "X" "*"
  # Ambig code J, and special codes O and U are returned as 22 (X). Gaps are returned as NA or removed if na.rm = T
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
  if(is.list(x)) lapply(x, fun) else fun(x)
}

AA2duovigesimal <- function(x, na.rm = FALSE){
  # for use with Gonnet matrix
  # return order "C" "S" "T" "P" "A" "G" "N" "D" "E" "Q" "H" "R" "K" "M" "I" "L" "V" "F" "Y" "W" "X" "*"
  # Ambig codes B, J and Z, special codes O and Z, are returned as 20 (X), and gaps are returned as NA
  fun <- function(v){
    bits <- c(67, 83, 84, 80, 65, 71, 78, 68, 69, 81, 72, 82, 75, 77, 73, 76, 86, 70, 89, 87, 88, 42, 66, 74, 79, 85, 90)
    ints <- c(0:21, rep(20, 5))
    res <- ints[match(as.numeric(v), bits)]
    attributes(res) <- attributes(unclass(v))
    isnares <- is.na(res) #placeholder for ambig treatment
    if(na.rm) if(any(isnares)) res <- res[!isnares]
    return(res)
  }
  if(is.list(x)) lapply(x, fun) else fun(x)
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
