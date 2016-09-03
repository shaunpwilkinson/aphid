#' Derive a profile HMM from a multiple sequence alignment.
#'
#' \code{derivePHMM} generates a profile HMM from a given multiple sequence alignment.
#'
#' This function is the homologue of \code{hmmbuild} in HMMER3.1 (http://hmmer.org/) and
#' \code{modelfromalign} and \code{buildmodel} in the SAM package
#' (https://compbio.soe.ucsc.edu/sam.html).
#'
#' @param x a character matrix of aligned sequences.
#' @param residues one of the character strings 'autodetect' (default), aminos' or
#' 'bases', or a case sensitive character vector specifying the residue
#' alphabet. Specifying the symbol type ('aminos' or 'bases') can increase speed
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
derivePHMM <- function(x, residues = "autodetect", gapchar = "-",
                   pseudocounts = "background", logspace = TRUE,
                   qa = NULL, qe = NULL,
                   inserts = "threshold", threshold = 0.5,
                   lambda = 0, DI = TRUE, ID = TRUE){
  if(!(pseudocounts %in% c("background", "Laplace", "none"))){
    stop("invalid pseudocounts argument")
  }
  ### include option to enter manually

  if(!(is.matrix(x))) stop("invalid inout object type, x must be a matrix")
  residues <- alphabet(x, residues = residues, gapchar = gapchar)
  nres <- length(residues)
  n <- nrow(x)
  m <- ncol(x)
  states <- c("D", "M", "I")
  transitions <- c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  # background emission probabilities (qe)
  if(is.null(qe)){
    allecs <- tab(x, residues = residues) + 1
    #allecs <- apply(x, 2, tab, residues)
    #allecs <- apply(allecs, 1, sum) + 1
    qe <- allecs/sum(allecs)
  }else{
    if(!(is.vector(qe) & length(qe) == length(residues))) stop("qe invalid")
    if(is.null(names(qe))) stop("qe argument is missing names attribute")
    if(!(identical(names(qe), residues))) stop("qe names must match residues")
    if(all(qe <= 0)) qe <- exp(qe)
    if(!(round(sum(qe), 2) == 1)) stop("background emission probs (qe) must sum to 1")
    if(any(qe == 0)) warning("qe contains at least one zero-probabiliy")
    qe <- qe/sum(qe) #account for rounding errors
  }

  # designate insert-columns
  gaps <- x == gapchar
  if(identical(inserts, "threshold")){
    inserts <- apply(gaps, 2, sum) > threshold * n
  } else if(identical(inserts, "map")){
    inserts <- !map(x, residues = residues, gapchar = gapchar,
                    pseudocounts = pseudocounts, qa = qa, qe = qe)
  } else if(!(mode(inserts) == "logical" & length(inserts) == ncol(x))){
    stop("invalid inserts argument")
  }
  l <- sum(!inserts) # PHMM length (excluding B & E positions)

  # emission counts
  ecs <- apply(x[, !inserts], 2, tab, residues)
  dimnames(ecs) <- list(symbol = residues, position = 1:l)

  #transitions
  xtr <- matrix(nrow = n, ncol = m)
  insertsn <- matrix(rep(inserts, n), nrow = n, byrow = T)
  xtr[gaps & !insertsn] <- 0L # Delete
  xtr[!gaps & !insertsn] <- 1L # Match
  xtr[!gaps & insertsn] <- 2L # Insert
  xtr <- cbind(1L, xtr, 1L) # append begin and end match states
  tcs <- tab9C(xtr, modules = l + 2)
  transtotals <- apply(tcs, 1, sum) + 1 # forced addition of Laplacian pseudos

  # background transition probs
  if(is.null(qa)){
    if(!DI) transtotals["DI"] <- 0
    if(!ID) transtotals["ID"] <- 0
    #qa <- matrix(transtotals, nrow = 3, byrow = TRUE)
    #dimnames(qa) <- list(from = c("D", "M", "I"), to = c("D", "M", "I"))
    # qa <- qa/apply(qa, 1, sum)
    qa <- transtotals/sum(transtotals)
  }else{
    # if(!(is.matrix(qa) & all(dim(qa) == 3))) stop("qa must be a 3 x 3 matrix")
    if(!is.vector(qa) | length(qa) != 9) stop("qa must be a numeric vector of length 9")
    if(all(qa <= 0)) qa <- exp(qa)
    # if(!all(round(apply(qa, 1, sum), 2) == 1)) stop("rows of qa matrix must sum to 1")
    if(round(sum(qa), 2) != 1) stop("qa vector must sum to 1")
    if(!DI & (qa[3] != 0)){
      stop("DI is set to FALSE but delete -> insert transitions are assigned non-zero
            probabilities in qa vector. Change DI to TRUE or change qa[3] to zero")
    }
    if(!ID & (qa[7] != 0)){
      stop("ID is set to FALSE but insert -> delete transitions are assigned non-zero
            probabilities in qa vector. Change ID to TRUE or change qa[7] to zero")
    }
    # qa <- qa/apply(qa, 1, sum) #account for rounding errors
    qa <- qa/sum(qa) #account for rounding errors
  }

  # transition and emission probability arrays A & E

  # transprops <- transtotals/sum(transtotals)
  #transprops <- as.vector(qa)
  #notfromD <- c(rep(0, 3), rep(1, 6))
  #nottoD <- rep(c(0, 1, 1), times = 3)
  if(pseudocounts == 'background'){
    ecs <- ecs + qe * nres
    tcs <- tcs + qa * if(DI & ID) 9 else if (DI | ID) 8 else 7
    # tcs[, 2:l] <- tcs[, 2:l] + transprops * if(DI) 9 else 7
    # tcs[, 1] <- tcs[, 1] + transprops * notfromD * if(DI) 6 else 5
    ### doesn't sum to 6 or 5
    # tcs[, l + 1] <- tcs[, l + 1] + transprops * nottoD * if(DI) 6 else 5
    ### doesn't sum to 6 or 5
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
map <- function(x, residues = "autodetect", gapchar = "-",
                pseudocounts = "background",
                lambda = 0, qa = NULL, qe = NULL){
  if(!(pseudocounts %in% c("background", "Laplace", "none"))) {
    stop("invalid pseudocounts")
  }
  L <- ncol(x)
  n <- nrow(x)
  residues <- alphabet(x, residues = residues, gapchar = gapchar)
  nres <- length(residues)
  S <- sigma <- c(0, rep(NA, L + 1))
  ecs <- apply(x, 2, tab, residues)
  ### qa & qe checks?
  if(is.null(qe)){
    allecs <- apply(ecs, 1, sum) + 1
    qe <- allecs/sum(allecs)
  }
  if(is.null(qa)){
    gaps <- x == gapchar
    inserts <- apply(gaps, 2, sum) > 0.5 * nrow(x)
    xtr <- matrix(nrow = nrow(x), ncol = ncol(x))
    insertsn <- matrix(rep(inserts, n), nrow = n, byrow = T)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- tab9C(xtr, modules = sum(!inserts) + 2)
    transtotals <- apply(tcs, 1, sum) + 1 # forced addition of Laplacian pseudos
    #if(!DI) transtotals[c(3, 7)] <- 0 ### need to work out DI strategy
    # qa <- matrix(transtotals, nrow = 3, byrow = TRUE)
    # dimnames(qa) <- list(from = c("D", "M", "I"), to = c("D", "M", "I"))
    # qa <- qa/apply(qa, 1, sum)
    qa <- transtotals/sum(transtotals)
  }
  # calculate Mj for j = 1 ... L + 1
  if(pseudocounts == "background"){
    ecs2 <- ecs + qe * nres
  } else if(pseudocounts == "Laplace"){
    ecs2 <- ecs + 1
  } else if(pseudocounts == "none"){
    ecs2 <- ecs
  } else stop("invalid pseudocounts argument")
  emtots <- apply(ecs2, 2, sum)
  term2 <- log(t(t(ecs2)/emtots))
  M <- apply(ecs * term2, 2, sum)
  ### could speed this bit up using sparse matrices
  M <- c(0, M, 0)
  # calculate Tij for j = 1 ... L + 1, -1 < i < j
  Tij <- Iij <- matrix(0, nrow = L + 2, ncol = L + 2)
  x2 <- cbind(residues[1], x, residues[1]) # dummy match states at beginning and end
  for(j in 1:ncol(Tij)){
    for(i in 1:nrow(Tij)){
      if(i < j){
        salen <- j - i + 1 # sub-alignment length
        subalig <- x2[, i:j, drop = FALSE]
        if(salen > 2){
          #insres = vector of all residues emitted while in the insert state
          insres <- subalig[, 2:(salen - 1)][subalig[, 2:(salen - 1)] != gapchar]
          # inscounts = insres tabulated
          inscounts <- tab(insres, residues = residues)
          Iij[i, j] <- sum(inscounts * log(qe))
        }
        satr <- matrix(nrow = n, ncol = salen) # subalignment as ternary
        rownames(satr) <- rownames(subalig) # not really necessary
        satr[, 1][subalig[, 1] != gapchar] <- 1L
        satr[, 1][subalig[, 1] == gapchar] <- 0L
        satr[, salen][subalig[, salen] != gapchar] <- 1L
        satr[, salen][subalig[, salen] == gapchar] <- 0L
        if(salen > 2) satr[, 2:(salen - 1)][subalig[, 2:(salen - 1)] != gapchar] <- 2L
        cxy <- tab9C(satr, modules = 2)
        # replace D -> in begin state w/ 0
        cxy <- apply(cxy, 1, sum) # transition totals
        # cxy <- matrix(cxy, nrow = 3, byrow = TRUE)
        # dimnames(cxy) <- list(from = c("D", "M", "I"), to = c("D", "M", "I"))
        if(pseudocounts == "background"){
          cxy2 <- cxy + qa * 9
        } else if(pseudocounts == "Laplace"){
          cxy2 <- cxy + 1
          ### also need DI option
        } else{
          cxy2 <- cxy
        }
        cxy2 <- matrix(cxy2, 3, byrow = TRUE, dimnames = list(from = c("D", "M", "I"), to = c("D", "M", "I")))
        axy <- cxy2/apply(cxy2, 1, sum)
        Tij[i, j] <- sum(cxy * log(axy))
        ### note only 6 poss states in col 0 etc
      }
    }
  }
  tmp <- t(t(Tij + Iij) + M + lambda)
  tmp[lower.tri(tmp, diag = TRUE)] <- NA
  for(j in 2:(L + 2)){
    sigma[j] <- whichismax(S[1:(j - 1)] + tmp[1:(j - 1), j])
    S[j] <- (S[1:(j - 1)] + tmp[1:(j - 1), j])[sigma[j]]
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

