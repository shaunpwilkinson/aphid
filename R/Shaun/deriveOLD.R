#' Derive profile HMMs from sequence alignments.
#' 
#' \code{derivePHMM} converts a multiple sequence alignment to a profile HMM.
#' 
#' @param x a character matrix of aligned sequences.
#' @param residues a character string specifying whether the residues in x
#' represent 'aminos' or 'nucleotides'.
#' @param gapChar the character that represents gaps in the alignment matrix.
#' @param method the method used to account for the possible absence of certain
#' character states in the alignment columns. Only Laplacean pseudocounts
#' ("method = "Laplace") are supported in the current version.
#' @param qe a named vector of background emission probabilities the same length as
#' the number of residues (i.e. 4 for nucleotides and 20 for amino acids) and with 
#' names corresponding to the residue alphabet (i.e. A C D E F G H I K L M N P Q 
#' R S T V W Y for amino acids or a c g t for nucleotides). 
#' If NULL, background emission probabilities are derived from the alignment.
#' @param qa a 3 x 3 matrix of background transition probabilities with 
#' \code{dimnames(qa) = list(from = c('M', 'I', 'D'), to = c('M', 'I', 'D'))}, 
#' where M, I and D represent matches, inserts and deletions, respectively. If
#' NULL background transition probabilities are derived from the alignment.
#' @param quiet logical argument indicating whether warnings should be printed 
#' to the console.
#'    
derivePHMM <- function(x, residues = "aminos", gapChar = "-", method = "Laplace", 
                       qe = NULL, qa = NULL, quiet = FALSE){
  n <- nrow(x)
  gaps <- x == gapChar
  inserts <- apply(gaps, 2, sum) > 0.5 * n
  m <- sum(!inserts)
  states <- c("M", "I", "D")
  if(identical(residues, "aminos")) residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  if(identical(residues, "nucleotides")) residues <- c("a", "c", "g", "t")
  symblen <- length(residues)
  tpc <- cbind(rep(states, each = 3), rep(states, 3)) # trans pseudocounts
  epc <- matrix(rep(residues, m), ncol = m) # emission pseudocounts
  counts <- apply(rbind(x[, !inserts], epc), 2, table, exclude = gapChar) - 1
  if(!(is.null(qe))){
    if(!(is.vector(qe) & length(qe) == length(residues))){
      stop("background emissions argument (qe) must be a vector of the same length 
           as the residue alphabet (i.e. 20 for amino acids or 4 for nucleotides)")
    }
    if(is.null(names(qe))) stop("background emissions argument (qe) must be named")
    ###maybe change to warning and add names alphabetically?###
    qe <- qe[order(names(qe))]
    if(!(all(names(qe) == residues))){
      stop("background emissions argument (qe) must have names attribute matching
           the residue alphabet (i.e. A C D E F G H I K L M N P Q R S T V W Y for
           amino acids or a c g t for nucleotides)")
    }
    if(!(sum(qe) == 1)) stop("background emissions vector (qe) must sum to 1")
    if(any(qe == 0)) warning("background emissions (qe) contain zero probabilities")
  }else{
    qe <- apply(counts, 1, sum)
    if(any(qe == 0) | method == "Laplace"){
      if(!quiet & method == 'background'){
        warning("alignment gave zero emission probabilities for at least 
                  one residue, Laplacian pseudocounts added")
      } 
      qe <- qe + 1
    }
    qe <- qe/sum(qe)
  }
  if(method == "Laplace"){
    counts <- counts + 1
    totals <- apply(counts, 2, sum)
    E <- t(t(counts)/totals)
  }else if (method == "background"){
    totals <- apply(counts, 2, sum)
    E <- t(t(counts + matrix(rep(qe * symblen, m), ncol = m))/(totals + symblen))
  }
  dimnames(E) <- list(residues = residues, position = 1:m)
  pathfinder <- function(v, inserts){
    v[v != gapChar & !inserts] <- "M"
    v[v == gapChar & !inserts] <- "D"
    v[v != gapChar & inserts] <- "I"
    return(c("M", v[v != gapChar], "M"))
  }
  paths <- apply(x, 1, pathfinder, inserts = inserts)
  if(is.matrix(paths)) paths <- as.list(data.frame(paths, stringsAsFactors = FALSE))
  A <- array(0, dim = c(3, m + 1, 3))
  dimnames(A) <- list(from = states, 0:m, to = states)
  tmp <- list()
  for(i in 1:(m + 1)){
    tmp[[i]] <- matrix(ncol = 2, nrow = 0)
    for(j in 1:n){
      proceed <- TRUE
      while(proceed){
        trans <- paths[[j]][1:2]
        paths[[j]] <- paths[[j]][-1]
        tmp[[i]] <- rbind(tmp[[i]], trans)
        proceed <- trans[2] == "I"
      }
    }
  }
  if(!(is.null(qa))){
    if(!(is.matrix(qa) & all(dim(qa) == 3) & all(apply(qa, 1, sum) == 1))){
      stop("background transitions argument (qa) must be a 3 x 3 matrix
           with all rows summing to 1") 
    } 
    dimnames(qa) <- list(states, states)
  }else{
    qa <- matrix(NA, 3, 3, dimnames = list(states, states))
    mergedtmp <- tmp[[1]]
    for(i in 2:(m + 1)) mergedtmp <- rbind(mergedtmp, tmp[[i]])
    vectmp <- unname(apply(mergedtmp, 1, paste0, collapse = ""))
    vectpc <- apply(tpc, 1, paste0, collapse = "")
    all.trans.types.in.alig <- all(vectpc %in% vectmp)
    if(!all.trans.types.in.alig | method == "Laplace"){
      if(!quiet & method == 'background'){
        warning("alignment gave zero transition probabilities for at 
                  least one transition type, Laplacian pseudocounts added")
      } 
      mergedtmp <- rbind(mergedtmp, tpc)
    }
    for(k in states){
      tfk <- mergedtmp[, 1] == k
      nk <- sum(tfk)
      for(l in states){
        qa[k, l] <- sum(tfk & mergedtmp[, 2] == l)/nk
      }
    }
  }
  if(method == "Laplace"){
    for(i in 1:(m + 1)){
      pc <- if(i == 1) tpc[1:6,] else if(i == m + 1) tpc[-c(3, 6, 9),] else tpc
      tmp[[i]] <- rbind(tmp[[i]], pc)
      for(k in states){
        tfk <- tmp[[i]][, 1] == k
        nk <- sum(tfk)
        for(l in states){
          transprob <- sum(tfk & tmp[[i]][, 2] == l)/nk
          if(!is.nan(transprob)) A[k, i, l] <- transprob
        }
      }
    }
  }else if(method == 'background'){
    for(i in 1:(m + 1)){
      symblen <- if(i == 1 | i == m + 1) 6 else 9
      for(k in states){
        tfk <- tmp[[i]][, 1] == k
        nk <- sum(tfk)
        for(l in states){
          if((i == 1 & k == "D") | (i == m + 1 & l == "D")){
            transprob <- 0
          }else{
            transprob <- (sum(tfk & tmp[[i]][, 2] == l) + 
                            (qa[k, l] * symblen))/(nk + symblen)
          }
          A[k, i, l] <- if(is.nan(transprob)) 0 else transprob 
        }
      }
    }
  }
  res <- structure(list(E = E, A = A, qe = qe, qa = qa, inserts = inserts), 
                   class = "PHMM")
  res
}
