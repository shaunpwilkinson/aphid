#' Derive profile HMMs from sequence alignments.
#' 
#' \code{derive} converts a multiple sequence alignment to a profile HMM.
#' 
#' @param x a character matrix of aligned sequences.
#' @param residues one of the character strings 'auto' (default), aminos' or 
#' 'bases', or a case sensitive character vector matching the residue 
#' alphabet. Specifying the symbol type ('aminos' or 'bases') can increase speed 
#' for larger alignments.
#' @param gapchar the character used to represent gaps in the alignment matrix.
#' @param method the method used to account for the possible absence of certain
#' character states in the alignment columns. If \code{method = background} and 
#' some transition/emission types are absent from the alignment, Laplacean 
#' pseudocounts are automatically added.
#' @param qe an optional named vector of background emission probabilities the 
#' same length as the number of residues (i.e. 4 for nucleotides and 20 for 
#' amino acids) and with names corresponding to the residue alphabet 
#' (i.e. A C D E F G H I K L M N P Q R S T V W Y for amino acids or 
#' a c g t for nucleotides). If NULL, background emission probabilities are 
#' derived from the alignment.
#' @param qa an optional named 3 x 3 matrix of background transition probabilities with 
#' \code{dimnames(qa) = list(from = c('D', 'M', 'I'), to = c('D', 'M', 'I'))}, 
#' where M, I and D represent matches, inserts and deletions, respectively. If
#' NULL background transition probabilities are derived from the alignment.
#' 
derive <- function(x, residues = 'auto', gapchar = "-",
                   method = 'background', qe = NULL, qa = NULL){
  if(!(method %in% c('background', 'Laplace'))) stop("invalid method")
  if(identical(residues, "auto")) residues <- sort(unique(as.vector(x)))
  if(identical(residues, "aminos")) residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  if(identical(residues, "bases")) residues <- c("a", "c", "g", "t")
  if(gapchar %in% residues) residues <- residues[residues != gapchar]
  if(!(length(residues) > 1)) stop("invalid residues argument")
  symblen <- length(residues)
  n <- nrow(x)
  m <- ncol(x)
  states <- c("D", "M", "I")
  gaps <- x == gapchar
  gapsn <- as.vector(t(gaps))
  inserts <- apply(gaps, 2, sum) > 0.5 * n
  l <- sum(!inserts) # PHMM length (excluding B & E positions)
  tab <- function(v){
    res <- rep(NA, symblen)
    for(i in 1:symblen) res[i] <- sum(v == residues[i], na.rm = TRUE)
    return(res)
  }
  #tab <- function(v) sapply(residues, function(x) sum(v == x, na.rm = TRUE))
  ecs <- apply(x[, !inserts], 2, tab) # emission counts
  dimnames(ecs) <- list(symbol = residues, position = 1:l)
  # background emission probabilities (qe)
  if(is.null(qe)){
    qe <- apply(ecs, 1, sum)
    if(any(qe == 0)) qe <- qe + 1
    qe <- qe/sum(qe)
  }else{
    if(!(is.vector(qe) & length(qe) == length(residues))) stop("qe invalid")
    if(is.null(names(qe))) stop("qe argument missing names attribute")
    if(!(identical(names(qe), residues))) stop("qe names must match residues")
    if(!(round(sum(qe), 2) == 1)) stop("qe vector must sum to 1")
    if(any(qe == 0)) warning("qe contains at least one zero probabiliy")
    qe <- qe/sum(qe) #account for rounding errors
  }
  #transitions
  xtr <- rep(NA, n * m) # code x to ternary numbering system 
  insertsn <- rep(inserts, n)
  xtr[gapsn & !insertsn] <- 0 # Delete
  xtr[!gapsn & !insertsn] <- 1 # Match
  xtr[!gapsn & insertsn] <- 2 # Insert
  xtr <- matrix(xtr, nrow = n, byrow = TRUE)
  xtr <- cbind(1, xtr, 1) # begin and end match states
  tupler <- function(v) apply(rbind(v[-(length(v))], v[-1]), 2, convert, 3, 10)
  xdc <- t(apply(xtr, 1, tupler)) + 1 # transition types coded to decimal
  tab <- function(v){
    res <- rep(NA, 9)
    for(i in 1:9) res[i] <- sum(v == i, na.rm = TRUE)
    res
  } # this is much faster than table(factor(v, levels = 1:9)))
  #its even faster than this:
  #tab <- function(v) sapply(1:9, function(x) sum(v == x, na.rm = TRUE))
  tcs <- apply(xdc, 2, tab)
  rownames(tcs) <- c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  if(sum(inserts) > 0){
    # merge insert states
    itp <- apply(rbind(c(F, inserts), c(inserts, F)), 2, convert, 2, 10)
    ist <- which(itp == 1) #insert start positions
    ien <- which(itp == 2) #insert end positions
    tmp <- mapply(":", ist, ien, SIMPLIFY = FALSE)
    tcs[, ist] <- sapply(tmp, function(e) apply(tcs[, e], 1, sum))
    tcs <- tcs[, itp < 2]
  }
  # background transition probs
  if(is.null(qa)){
    #qa <- addmats(tcs)
    qa <- apply(tcs, 1, sum)
    if(any(qa == 0)) qa <- qa + 1
    qa <- matrix(qa, nrow = 3, byrow = TRUE)
    dimnames(qa) <- list(from = c("D", "M", "I"), to = c("D", "M", "I"))
    qa <- qa/sum(qa) 
    #qa <- qa/apply(qa, 1, sum) #? check in book
  }else{
    if(!(is.matrix(qa) & all(dim(qa) == 3))) stop("qa must be a 3 x 3 matrix") 
    if(!(round(sum(qa), 2) == 1)) stop("elements of qa matrix should sum to 1")
    qa <- qa/sum(qa) #account for rounding errors
  }
  # transition and emission probability arrays A & E
  notfromD <- rep(c(0, 1, 1), times = 3)
  nottoD <- c(rep(0, 3), rep(1, 6))
  qavec <- as.vector(qa)
  if(method == 'background'){
    ecs <- ecs + qe * symblen
    tcs[, 2:l] <- tcs[, 2:l] + qavec * 9
    tcs[, 1] <- tcs[, 1] + qavec * notfromD * 6
    tcs[, l + 1] <- tcs[, l + 1] + qavec * nottoD * 6 
  }else if(method == 'Laplace'){
    ecs <- ecs + 1
    tcs[, 2:l] <- tcs[, 2:l]  + 1
    tcs[, 1] <- tcs[, 1] + notfromD
    tcs[, l + 1] <- tcs[, l + 1] + nottoD 
  }
  E <- t(t(ecs)/apply(ecs, 2, sum))
  toD <- tcs[1:3,]
  toM <- tcs[4:6,]
  toI <- tcs[7:9,]
  tot <- toD + toM + toI
  A <- cbind(toD/tot, toM/tot, toI/tot)
  A <- unlist(A, use.names = FALSE)
  A[is.nan(A)] <- 0 
  A <- array(A, dim = c(3, l + 1, 3))
  dimnames(A) <- list(from = states, 0:l, to = states)
  #inserts <- if(any(inserts)) which(inserts) else NULL
  inslens <- insertlengths(!inserts)
  #which alignment columns correspond to which model positions?
  alignment <- which(!inserts)
  names(alignment) <- 1:l
  res <- structure(list(E = E, A = A, qe = qe, qa = qa, size = l,
                        insertlengths = inslens, 
                        alignment = alignment), 
                   class = "PHMM")
  return(res)
}




