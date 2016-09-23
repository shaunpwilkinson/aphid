#' Derive a profile HMM from a multiple sequence alignment.
#'
#' \code{derivePHMM} generates a profile HMM from a given multiple sequence alignment.
#'
#' This function is the homologue of \code{hmmbuild} in HMMER3.1 (http://hmmer.org/) and
#' \code{modelfromalign} and \code{buildmodel} in the SAM package
#' (https://compbio.soe.ucsc.edu/sam.html).
#'
#' @param x a character matrix of aligned sequences.
#' @param seqweights integer 1 or numeric vector of sequence weights used to
#' derive the model. The sum of these weights should be equal to the number of
#' sequences in the alignment.
#' @param residues one of the character strings 'autodetect' (default),
#' 'AA' (amino acids) or
#' 'DNA' (A, C, G, and T), or a case sensitive character vector specifying the
#' residue
#' alphabet. Specifying the residue alphabet (i.e 'DNA' or 'AA') can increase speed
#' for larger alignments. Note that setting \code{residues = 'autodetect'} will not
#' detect rare residues that are not present in the alignment and thus will
#' not assign them emission probabilities.
#' @param gapchar the character used to represent gaps in the alignment matrix.
#' @param pseudocounts the method used to account for the possible absence of certain
#' character states in the alignment. Currently only \code{"Laplace"}
#' and \code{"background"} are supported. If \code{pseudocounts = "background"} and
#' some transition/emission types are absent from the alignment, Laplacean
#' pseudocounts are automatically added.
#' @param logspace logical indicating whether the emission and transition
#' probability on the returned model should be logged.
#' @param qe an optional named vector of background emission probabilities the
#' same length as the residue alphabet (i.e. 4 for nucleotides and 20 for
#' amino acids) and with  corresponding names
#' (i.e. A C D E F G H I K L M N P Q R S T V W Y for amino acids or
#' A C G T for nucleotides). If \code{NULL}, background emission probabilities are
#' generated from the alignment.
#' @param qa an optional named 3 x 3 matrix of background transition probabilities with
#' \code{dimnames(qa) = list(from = c('D', 'M', 'I'), to = c('D', 'M', 'I'))},
#' where M, I and D represent matches, inserts and deletions, respectively. If
#' \code{NULL}, background transition probabilities are generated from the alignment.
#' @param inserts model construction method. Either 'threshold' (default),'map',
#' or if insert columns are to be specified manually, a logical vector of length ncol(x),
#' specifying the columns of x that should be designated
#' as inserts (\code{TRUE} for insert columns, \code{FALSE} for match states).
#' If \code{inserts = "threshold"} columns are assigned as insert states if the proportion of
#' gap characters in the column is higher than a specified threshold (defaults to 0.5).
#' If set to 'map' insert columns are found using the maximum \emph{a posteriori}
#' dynamic programming method outlined in Durbin et al. (1998).
#' @param threshold the maximum proportion of gaps for an alignment column
#' to be considered a module in the PHMM (defaults to 0.5). Only applicable for
#' \code{inserts = "threshold"}.
#' @param lambda penalty parameter used to favour models with fewer match states. Equivalent
#' to the log of the prior probability of marking each column (Durbin et al. 1998, pg 124).
#' Only applicable for \code{inserts = "map"}.
#' @references Durbin..
#'
derivePHMM <- function(x, seqweights = 1, residues = "autodetect",
                       gapchar = "-", pseudocounts = "background",
                       logspace = TRUE, qa = NULL, qe = NULL,
                       inserts = "threshold", threshold = 0.5,
                       lambda = 0, DI = TRUE, ID = TRUE){
  if(!(pseudocounts %in% c("background", "Laplace", "none"))){
    stop("invalid pseudocounts argument")
  }
  ### include option to enter manually
  if(!(is.matrix(x))) stop("invalid object type, x must be a matrix")
  DNA <- inherits(x, "DNAbin")
  if(DNA) gapchar <- as.raw(4)
  residues <- alphabet(x, residues = residues, gapchar = gapchar)
  nres <- length(residues)
  n <- nrow(x)
  m <- ncol(x)
  states <- c("D", "M", "I")
  transitions <- c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  if(identical(seqweights, 1)) {
    seqweights <- rep(1, n)
  }else{
    if(round(sum(seqweights), 2) != n){
      if(round(sum(seqweights), 2) == 1){
        seqweights <- seqweights * n
      }else stop("invalid seqweights argument")
    }
  }
  if(length(seqweights) != n) stop("invalid seqweights argument")
  # background emission probabilities (qe)
  if(is.null(qe)){
    if(DNA){
      allecs <- apply(x, 2, tabDNA, ambiguities = TRUE, seqweights = seqweights)
    } else{
      allecs <- apply(x, 2, tab, residues = residues, seqweights = seqweights)
    }
    allecs <- apply(allecs, 1, sum)
    qe <- (allecs + 1)/sum(allecs + 1)
  }else{
    if(!(is.vector(qe) & length(qe) == length(residues))) stop("qe invalid")
    if(is.null(names(qe))) stop("qe argument is missing names attribute")
    if(!(identical(names(qe), residues))) stop("qe names must match residues")
    if(all(qe <= 0)) qe <- exp(qe)
    if(!(round(sum(qe), 2) == 1)) stop("background emissions (qe) must sum to 1")
    if(any(qe == 0)) warning("at least one background emission probability is zero")
    qe <- qe/sum(qe) #account for rounding errors
  }
  # designate insert-columns
  gaps <- x == gapchar
  gapweights <- gaps * seqweights
  if(identical(inserts, "threshold")){
    inserts <- apply(gapweights, 2, sum) > threshold * n
  } else if(identical(inserts, "map")){
    inserts <- !map(x, residues = residues, gapchar = gapchar, seqweights = seqweights,
                    pseudocounts = pseudocounts, qa = qa, qe = qe)
  } else if(!(mode(inserts) == "logical" & length(inserts) == ncol(x))){
    stop("invalid inserts argument")
  }
  l <- sum(!inserts) # PHMM length (excluding B & E positions)
  # emission counts
  if(DNA){
    ecs <- apply(x[, !inserts], 2, tabDNA, ambiguities = TRUE, seqweights = seqweights)
  }else{
    ecs <- apply(x[, !inserts], 2, tab, residues = residues, seqweights = seqweights)
  }
  dimnames(ecs) <- list(residue = residues, position = 1:l)

  #transitions
  xtr <- matrix(nrow = n, ncol = m)
  insertsn <- matrix(rep(inserts, n), nrow = n, byrow = T)
  xtr[gaps & !insertsn] <- 0L # Delete
  xtr[!gaps & !insertsn] <- 1L # Match
  xtr[!gaps & insertsn] <- 2L # Insert
  xtr <- cbind(1L, xtr, 1L) # append begin and end match states
  tcs <- tab9C(xtr, seqweights = seqweights)
  alltcs <- apply(tcs, 1, sum) # forced addition of Laplacian pseudos

  # background transition probs
  if(is.null(qa)){
    if(!DI) alltcs["DI"] <- 0
    if(!ID) alltcs["ID"] <- 0
    qa <- (alltcs + 1)/sum(alltcs + 1)
  }else{
    if(!is.vector(qa) | length(qa) != 9) stop("qa must be a numeric vector of length 9")
    if(all(qa <= 0)) qa <- exp(qa)
    if(round(sum(qa), 2) != 1) stop("qa vector must sum to 1")
    if(!DI & (qa[3] != 0)){
      stop("DI is set to FALSE but delete -> insert transitions are assigned non-zero
            probabilities in qa vector. Change DI to TRUE or change qa['DI'] to zero")
    }
    if(!ID & (qa[7] != 0)){
      stop("ID is set to FALSE but insert -> delete transitions are assigned non-zero
            probabilities in qa vector. Change ID to TRUE or change qa['ID'] to zero")
    }
    # qa <- qa/apply(qa, 1, sum) #account for rounding errors
    qa <- qa/sum(qa) #account for rounding errors
  }

  # transition and emission probability arrays A & E
  if(pseudocounts == 'background'){
    ecs <- ecs + qe * nres
    tcs <- tcs + qa * (7 + sum(c(DI, ID)))
  }else if(pseudocounts == 'Laplace'){
    ecs <- ecs + 1
    tcs <- tcs + 1
    # tcs[, 2:l] <- tcs[, 2:l]  + 1
    # tcs[, 1] <- tcs[, 1] + notfromD
    # tcs[, l + 1] <- tcs[, l + 1] + nottoD
  }
  tcs[1:3, 1] <- tcs[c(1, 4, 7), l + 1] <- 0
  if(!DI) tcs[3, ] <- 0
  if(!ID) tcs[7, ] <- 0
  E <- t(t(ecs)/apply(ecs, 2, sum))
  A <- t(tcs)
  for(i in c(1, 4, 7)) A[, i:(i + 2)] <- A[, i:(i + 2)]/apply(A[, i:(i + 2)], 1, sum)
  A[1, 1:3] <- 0 # gets rid of NaNs caused by dividing by zero
  A <- t(A)

  inslens <- insertlengths(!inserts)
  #which alignment columns correspond to which model positions?
  alignment <- which(!inserts)
  names(alignment) <- 1:l
  if(logspace){
    A <- log(A)
    E <- log(E)
    qa <- log(qa)
    qe <- log(qe)
  }
  res <- structure(list(A = A, E = E, qa = qa, qe = qe, size = l,
                        insertlengths = inslens,
                        alignment = alignment),
                   class = "PHMM")
  return(res)
}

#' Optimal profile HMM construction.
#'
#' Assigns columns to insert states using the maximum \emph{a posteriori} method
#' outlined in Durbin et al. (1998).
#'
#' @param x a character matrix of aligned sequences.
#' @param seqweights integer 1 or numeric vector of sequence weights used to
#' derive the emission and trasnision counts. The sum of these weights should
#' be equal to the number of
#' sequences in the alignment.
#' @param residues one of the character strings 'autodetect' (default), aminos' or
#' 'bases', or a case sensitive character vector specifying the residue
#' alphabet. Specifying the symbol type ('aminos' or 'bases') can increase speed
#' for larger alignments. Note that setting \code{residues = 'autodetect'} will not
#' detect rare residues that are not present in the alignment and thus will
#' not assign them emission probabilities.
#' @param gapchar the character used to represent gaps in the alignment matrix.
#' @param pseudocounts the method used to account for the possible absence of certain
#' character states in the alignment. Currently only \code{"Laplace"},
#'\code{"background"} and \code{"none"} are supported. If \code{pseudocounts = "background"} and
#' some transition/emission types are absent from the alignment, Laplacean
#' pseudocounts are automatically added.
#' @param lambda penalty parameter used to favour models with fewer match states. Equivalent
#' to the log of the prior probability of marking each column (Durbin et al. 1998, pg 124).
#' @return a logical vector of length = ncol(x) indicating the columns to be assigned
#' as match states (\code{TRUE}) and those assigned as inserts (\code{FALSE}).
#' @references Durbin...
#'
map <- function(x, seqweights = 1, residues = "autodetect",
                gapchar = "-", pseudocounts = "background",
                lambda = 0, qa = NULL, qe = NULL){
  if(!(pseudocounts %in% c("background", "Laplace", "none"))) {
    stop("invalid pseudocounts")
  }
  L <- ncol(x)
  n <- nrow(x)
  DNA <- inherits(x, "DNAbin")
  if(DNA) gapchar = as.raw(4)
  residues <- alphabet(x, residues = residues, gapchar = gapchar)
  nres <- length(residues)
  transitions = c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  S <- sigma <- c(0, rep(NA, L + 1))
  if(identical(seqweights, 1)) {
    seqweights <- rep(1, n)
  }else{
    if(round(sum(seqweights), 2) != n){
      if(round(sum(seqweights), 2) == 1){
        seqweights <- seqweights * n
      }else stop("invalid seqweights argument")
    }
  }
  if(length(seqweights) != n) stop("invalid seqweights argument")
  ecs <- if(DNA) apply(x, 2, tabDNA, TRUE, seqweights) else
    apply(x, 2, tab, residues, seqweights)
  # ecs <- t(t(ecs) * seqweights)
  allecs <- apply(ecs, 1, sum)
  if(is.null(qe)) {
    qe <- log((allecs + 1)/sum(allecs + 1))
  }else if(all(qe >= 0 & qe <= 1) & round(sum(qe), 2) == 1){
    qe <- log(qe)
  } else if(all(qe <= 0) & round(sum(exp(qe)), 2) == 1){
    # qe <- qe
  } else stop("invalid qe argument")
  gaps <- x == gapchar
  #if(!DI) alltcs[c(3, 7)] <- 0 ### need to work out DI strategy
  if(is.null(qa)) {
    gapweights <- gaps * seqweights
    inserts <- apply(gapweights, 2, sum) > 0.5 * nrow(x)
    xtr <- matrix(nrow = nrow(x), ncol = ncol(x))
    insertsn <- matrix(rep(inserts, n), nrow = n, byrow = T)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- tab9C(xtr, seqweights) # modules = sum(!inserts) + 2)
    alltcs <- apply(tcs, 1, sum)
    qa <- log((alltcs + 1)/sum(alltcs + 1)) # force addition of Laplace pseudos
  }else if(all(qa >= 0 & qa <= 1) & round(sum(qa), 2) == 1){
    qa <- log(qa)
  } else if(all(qa <= 0) & round(sum(exp(qa)), 2) == 1){
    # qa <- qa
  } else stop("invalid qa argument")
  # calculate Mj for j = 1 ... L + 1
  if(pseudocounts == "background"){
    ecs2 <- ecs + exp(qe) * nres
  } else if(pseudocounts == "Laplace"){
    ecs2 <- ecs + 1
  } else if(pseudocounts == "none"){
    ecs2 <- ecs
  } else stop("invalid pseudocounts argument")
  term2 <- t(t(ecs2)/apply(ecs2, 2, sum))
  term2[ecs != 0] <- log(term2[ecs != 0]) # increase speed for conserved alignments
  M <- apply(ecs * term2, 2, sum)
  M <- c(0, M, 0)
  # calculate Tij for j = 1 ... L + 1, -1 < i < j
  Tij <- Iij <- matrix(0, nrow = L + 2, ncol = L + 2)
  notgaps <- cbind(TRUE, !gaps, TRUE) # logical matrix, ncol = L + 2, nrow = n
  ivec <- 1:(L + 1)
  jvec <- (L + 2):2 ###
  no.insertsi <- apply(!gaps, 1, sum) * seqweights
  allecsi <- allecs # sum(allecsi) == sum(no.insertsi)
  alltcsi <- structure(numeric(9), names = transitions) # empty container
  alltcsi["MM"] <- sum(seqweights[no.insertsi == 0]) # number of blank sequences
  alltcsi[c("MI", "IM")] <- n - alltcsi["MM"]
  alltcsi["II"] <- sum(allecs) - alltcsi["MI"] # equiv to sum(no.insertsi) - alltcsi["IM"], etc
  alphaxy <- if(pseudocounts == "background") exp(qa) * 9 else rep(1, 9)
  for(i in ivec){
    no.insertsj <- no.insertsi # numeric, length n
    allecsj <- allecsi # numeric, length nres (would be integer but for ambigs)
    alltcsj <- alltcsi # integer, length 9
    for(j in jvec){
      #print(i); print(j); print(alltcsj)
      cxy <- alltcsj + alphaxy
      axy <- cxy/c(rep(sum(cxy[1:3]), 3), rep(sum(cxy[4:6]), 3), rep(sum(cxy[7:9]), 3))
      Tij[i, j] <- sum(cxy * log(axy))
      if(j - i > 1){
        Iij[i, j] <- sum(allecsj * qe)
        #now update
        allecsj <- allecsj - ecs[, j - 2]
        no.insertsj <- no.insertsj - notgaps[, j - 1] * seqweights
        zeroinserts <- round(no.insertsj, 5) == 0 # logical, length n
        if(any(zeroinserts)){
          alltcsj["DD"] <- sum(seqweights[!notgaps[, i] & !notgaps[, j - 1] & zeroinserts])
          alltcsj["DM"] <- sum(seqweights[!notgaps[, i] & notgaps[, j - 1] & zeroinserts])
          alltcsj["MD"] <- sum(seqweights[notgaps[, i] & !notgaps[, j - 1] & zeroinserts])
          alltcsj["MM"] <- sum(seqweights[notgaps[, i] & notgaps[, j - 1] & zeroinserts])
          alltcsj["DI"] <- sum(seqweights[!notgaps[, i] & !zeroinserts])
          alltcsj["MI"] <- sum(seqweights[notgaps[, i] & !zeroinserts])
          alltcsj["ID"] <- sum(seqweights[!notgaps[, j - 1] & !zeroinserts])
          alltcsj["IM"] <- sum(seqweights[notgaps[, j - 1] & !zeroinserts])
          # equiv to (alltcsj["ID"] + alltcsj["IM"])
        }else{
          alltcsj[c("DD", "DM", "MD", "MM")] <- 0
          alltcsj["DI"] <- sum(seqweights[!notgaps[, i]])
          alltcsj["MI"] <- sum(seqweights[notgaps[, i]])
          alltcsj["ID"] <- sum(seqweights[!notgaps[, j - 1]])
          alltcsj["IM"] <- sum(seqweights[notgaps[, j - 1]])
        }
        alltcsj["II"] <- sum(no.insertsj) - (alltcsj["DI"] + alltcsj["MI"])
      }else{
        alltcsj["DD"] <- sum(seqweights[!notgaps[, i] & !notgaps[, j - 1]])
        alltcsj["DM"] <- sum(seqweights[!notgaps[, i] & notgaps[, j - 1]])
        alltcsj["MD"] <- sum(seqweights[notgaps[, i] & !notgaps[, j - 1]])
        alltcsj["MM"] <- sum(seqweights[notgaps[, i] & notgaps[, j - 1]])
        alltcsj[c("DI", "MI", "ID", "IM", "II")] <- 0
      }
    }
    if(i < L + 1){
      allecsi <- allecsi - ecs[, i]
      no.insertsi <- no.insertsi - notgaps[, i + 1] * seqweights
      zeroinserts <- round(no.insertsi, 5) == 0
      alltcsi[c("DD", "MD", "ID")] <- 0
      if(any(zeroinserts)){
        alltcsi["DM"] <- sum(seqweights[!notgaps[, i + 1] & zeroinserts])
        alltcsi["MM"] <- sum(seqweights[notgaps[, i + 1] & zeroinserts])
        alltcsi["DI"] <- sum(seqweights[!notgaps[, i + 1] & !zeroinserts])
        alltcsi["MI"] <- sum(seqweights[notgaps[, i + 1] & !zeroinserts])
        alltcsi["IM"] <- sum(seqweights[!zeroinserts])
      }else{
        alltcsi[c("DM", "MM")] <- 0
        alltcsi["DI"] <- sum(seqweights[!notgaps[, i + 1]])
        alltcsi["MI"] <- n - alltcsi["DI"] # sum(seqweights[notgaps[, i + 1]])
        alltcsi["IM"] <- n
      }
      alltcsi["II"] <- sum(no.insertsi) - alltcsi["IM"] ### because alltcsi["ID"] = 0
    }else{
      alltcsi[c("DD", "DI", "MD", "MI", "ID", "IM", "II")] <- 0
      alltcsi["DM"] <- sum(seqweights[!notgaps[, i + 1]])
      alltcsi["MM"] <- n - alltcsi["DM"]
    }
    jvec <- jvec[-length(jvec)] # gives upper.tri
  }
  #tmp <- t(t(Tij + Iij) + M + lambda)
  #tmp[lower.tri(tmp, diag = TRUE)] <- NA
  for(j in 2:(L + 2)){
    tmp <- S[1:(j - 1)] + Tij[1:(j - 1), j] + Iij[1:(j - 1), j] + M[j] + lambda
    sigma[j] <- whichismax(tmp)#tmp[1:(j - 1), j])
    S[j] <- tmp[sigma[j]] #tmp[1:(j - 1), j])[sigma[j]]
  }
  res <- structure(logical(L + 2), names = 0:(L + 1))
  res[L + 2] <- TRUE
  j <- sigma[L + 2]
  while(j > 0){
    res[j] <- TRUE
    j <- sigma[j]
  }
  res <- res[-(c(1, L + 2))]
  return(res)
}

