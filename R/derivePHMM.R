#' Derive a profile HMM from a multiple sequence alignment.
#'
#' \code{derive.PHMM} generates a profile HMM from a given multiple sequence alignment.
#'
#' This function is the homologue of \code{hmmbuild} in HMMER3.1 (http://hmmer.org/) and
#' \code{modelfromalign} and \code{buildmodel} in the SAM package
#' (https://compbio.soe.ucsc.edu/sam.html).
#'
#' @param x a character matrix of aligned sequences.
#' @param seqweights either NULL (default; all sequences are given a weight of 1)
#' or a numeric vector of sequence weights used to
#' derive the model. The sum of these weights should be equal to the number of
#' sequences in the alignment, so that mean(seqweights) = 1.
#' @param residues either NULL (default; emitted residues are automatically
#' detected from the alignment), or a case sensitive character vector specifying the
#' residue alphabet (e.g. \code{c("A", "T", "G", "C")} for DNA).
#' Note that the former option can be slow for large character matrices;
#' therefore specifying the residue alphabet can increase speed in these cases.
#' Also note that the default setting \code{residues = NULL} will not
#' detect rare residues that are not present in the alignment, and thus will
#' not assign them emission probabilities.
#' @param gapchar the character used to represent gaps in the alignment matrix.
#' Not required for \code{"DNAbin"} or \code{"AAbin"} objects.
#' @param pseudocounts the method used to account for the possible absence of certain
#' character states in the alignment. Currently only \code{"Laplace"}
#' and \code{"background"} or \code{"none"} are supported.
#' If \code{pseudocounts = "background"} and
#' some transition/emission types are absent from the alignment, Laplacean
#' pseudocounts are automatically added.
#' @param logspace logical indicating whether the emission and transition
#' probability on the returned model should be logged.
#' @param qe an optional named vector of background emission probabilities the
#' same length as the residue alphabet (i.e. 4 for nucleotides and 20 for
#' amino acids) and with corresponding names
#' (e.g. \code{c("A", "T", "G", "C")} for DNA).
#' If \code{NULL}, background emission probabilities are generated from the alignment.
#' @param qa an optional named 9-element vector of background transition probabilities with
#' \code{dimnames(qa) = c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")},
#' where M, I and D represent match, insert and delete states, respectively. If
#' \code{NULL}, background transition probabilities are generated from the alignment.
#' @param inserts model construction method. Either 'threshold' (default),'map', 'inherited',
#' or if insert columns are to be specified manually, a logical vector of length ncol(x),
#' specifying the columns of x that should be designated
#' as inserts (\code{TRUE} for insert columns, \code{FALSE} for match states).
#' If \code{inserts = "threshold"} columns are assigned as insert states if the proportion of
#' gap characters in the column is higher than a specified threshold (defaults to 0.5).
#' If set to 'map' insert columns are found using the maximum \emph{a posteriori}
#' dynamic programming method outlined in Durbin et al. (1998). If set to 'inherited' the
#' function looks for an 'inserts' attribute from the alignment.
#' @param threshold the maximum proportion of gaps for an alignment column
#' to be considered a module in the PHMM (defaults to 0.5). Only applicable for
#' \code{inserts = "threshold"}. Note that the maximum \emph{a posteriori} method works
#' poorly for small alignments so the 'threshold' method is automatically used when
#' the number of sequences is fewer than 5.
#' @param lambda penalty parameter used to favour models with fewer match states. Equivalent
#' to the log of the prior probability of marking each column (Durbin et al. 1998, pg 124).
#' Only applicable for \code{inserts = "map"}.
#' @return an object of class \code{"PHMM"}
#' @references Durbin et al. 1998.
#' @name derive.PHMM
#' @export
#'
derive.PHMM <- function(x, seqweights = "Gerstein", wfactor = 1, k = 5, residues = NULL,
                        gapchar = "-", endchar = "?", pseudocounts = "background",
                        logspace = TRUE, qa = NULL, qe = NULL,
                        inserts = "map", threshold = 0.5, lambda = 0, deltaLL = 1E-07,
                        DI = FALSE, ID = FALSE, omit.endgaps = TRUE,
                        name = NULL, description = NULL, compo = FALSE, consensus = FALSE,
                        seeds = "random", refine = "Viterbi",
                        maxiter = if(refine == "Viterbi") 10 else 100,
                        cpp = TRUE, quiet = FALSE, ...){
  UseMethod("derive.PHMM")
}

#' @rdname derive.PHMM
#' @export
#'
derive.PHMM.DNAbin <- function(x, seqweights = "Gerstein", wfactor = 1, k = 5, residues = NULL,
                               gapchar = "-", endchar = "?", pseudocounts = "background",
                               logspace = TRUE, qa = NULL, qe = NULL,
                               inserts = "map", threshold = 0.5, lambda = 0, deltaLL = 1E-07,
                               DI = FALSE, ID = FALSE, omit.endgaps = TRUE,
                               name = NULL, description = NULL,
                               compo = FALSE, consensus = FALSE, seeds = "random", refine = "Viterbi",
                               maxiter = if(refine == "Viterbi") 10 else 100,
                               cpp = cpp, quiet = FALSE, ...){ ##TODO don't need gapchar, residues etc?
  if(is.list(x)){
    derive.PHMM.list(x, seeds = seeds, refine = refine, maxiter = maxiter,
                     seqweights = seqweights, wfactor = wfactor,
                     k = k, residues = residues, gapchar = gapchar, pseudocounts = pseudocounts,
                     logspace = logspace, qa = qa, qe = qe,
                     inserts = inserts, lambda = lambda, DI = DI, ID = ID,
                     threshold = threshold, deltaLL = deltaLL, omit.endgaps = omit.endgaps,
                     name = name, description = description, compo = compo, consensus = consensus,
                     cpp = cpp, quiet = quiet, ... = ...)
  }else{
    derive.PHMM.default(x, seqweights = seqweights, wfactor = wfactor, k = k, residues = residues,
                        gapchar = gapchar, endchar = endchar, pseudocounts = pseudocounts,
                        logspace = logspace, qa = qa, qe = qe,
                        inserts = inserts, threshold = threshold,
                        lambda = lambda, DI = DI, ID = ID, omit.endgaps = omit.endgaps,
                        name = name, description = description,
                        compo = compo, consensus = consensus, cpp = cpp, quiet = quiet)
  }
}


#' @rdname derive.PHMM
#' @export
#'
derive.PHMM.AAbin <- function(x, seqweights = "Gerstein", wfactor = 1, k = 5, residues = NULL,
                              gapchar = "-", endchar = "?", pseudocounts = "background",
                              logspace = TRUE, qa = NULL, qe = NULL,
                              inserts = "map", threshold = 0.5, lambda = 0, deltaLL = 1E-07,
                              DI = FALSE, ID = FALSE, omit.endgaps = TRUE,
                              name = NULL, description = NULL,
                              compo = FALSE, consensus = FALSE, seeds = "random", refine = "Viterbi",
                              maxiter = if(refine == "Viterbi") 10 else 100,
                              cpp = TRUE, quiet = FALSE, ...){
  if(is.list(x)){
    derive.PHMM.list(x, seeds = seeds, refine = refine, maxiter = maxiter,
                     seqweights = seqweights, wfactor = wfactor,
                     k = k, residues = residues, gapchar = gapchar, pseudocounts = pseudocounts,
                     logspace = logspace, qa = qa, qe = qe,
                     inserts = inserts, lambda = lambda, DI = DI, ID = ID,
                     threshold = threshold, deltaLL = deltaLL, omit.endgaps = omit.endgaps,
                     name = name, description = description, compo = compo, consensus = consensus,
                     cpp = cpp, quiet = quiet, ... = ...)
  }else{
    derive.PHMM.default(x, seqweights = seqweights, wfactor = wfactor, k = k, residues = residues,
                        gapchar = gapchar, endchar = endchar, pseudocounts = pseudocounts,
                        logspace = logspace, qa = qa, qe = qe,
                        inserts = inserts, threshold = threshold,
                        lambda = lambda, DI = DI, ID = ID, omit.endgaps = omit.endgaps,
                        name = name, description = description,
                        compo = compo, consensus = consensus, cpp = cpp, quiet = quiet)
  }
}

#' @rdname derive.PHMM
#' @export
#'
derive.PHMM.list <- function(x, seeds = "random", refine = "Viterbi",
                             maxiter = if(refine == "Viterbi") 10 else 100,
                             seqweights = "Gerstein", wfactor = 1,
                             k = 5, residues = NULL, gapchar = "-", pseudocounts = "background",
                             logspace = TRUE, qa = NULL, qe = NULL,
                             inserts = "map", lambda = 0, DI = FALSE, ID = FALSE,
                             threshold = 0.5, deltaLL = 1E-07, omit.endgaps = TRUE,
                             name = NULL, description = NULL, compo = FALSE, consensus = FALSE,
                             cpp = TRUE, quiet = FALSE, ...){
  nsq <- length(x)
  DNA <- is.DNA(x)
  AA <- is.AA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- alphadetect(x, residues = residues, gapchar = gapchar)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  for(i in 1:nsq) x[[i]] <- x[[i]][x[[i]] != gapchar]
  if(nsq > 2){
    if(!quiet) cat("Calculating pairwise distances\n")
    names(x) <- paste0("S", 1:nsq)
    if(identical(seeds, "random")){
      if(nsq > 20){
        nseeds <- min(nsq, 20 + ceiling(log(nsq, 2)))
        seeds <- sample(1:length(x), size = nseeds)
      }else seeds <- seq_along(x)
    }else if(identical(seeds, "all")){
      seeds <- seq_along(x)
    }
    qds <- kdistance(x[seeds], k = k, alpha = if(AA) "Dayhoff6" else if(DNA) NULL else residues)
    if(!quiet) cat("Building guide tree\n")
    guidetree <- as.dendrogram(hclust(qds, method = "average"))
    seedweights <- weight.dendrogram(guidetree, method = "Gerstein")[names(x)[seeds]]
    attachseqs <- function(tree, sequences){
      if(!is.list(tree)) attr(tree, "sequences") <- sequences[attr(tree, "label")]
      return(tree)
    }
    guidetree <- dendrapply(guidetree, attachseqs, sequences = x)
    progressive <- function(tree, ...){
      if(is.list(tree)){
        if(!is.null(attr(tree[[1]], "sequences")) & !is.null(attr(tree[[2]], "sequences"))){
          attr(tree, "sequences") <- align.default(attr(tree[[1]], "sequences"), attr(tree[[2]], "sequences"), ... = ...)
          attr(tree[[1]], "sequences") <- attr(tree[[2]], "sequences") <- NULL
        }
      }
      return(tree)
    }
    progressive1 <- function(tree, ...){
      tree <- progressive(tree, ... = ...)
      if(is.list(tree)) tree[] <- lapply(tree, progressive1, ... = ...)
      return(tree)
    }
    if(!quiet) cat("Aligning seed sequences\n")
    while(is.null(attr(guidetree, "sequences"))) guidetree <- progressive1(guidetree, ... = ...)
    msa1 <- attr(guidetree, "sequences")
    # ### next few lines need to be rewritten TODO
    # newick <- write.dendrogram(guidetree, edges = FALSE)
    # newick <- gsub(";", "", newick)
    # newick <- gsub("\\(", "align\\(", newick)
    # newick <- gsub("\\)", ", ... = ...\\)", newick)
    # #newick <- gsub("\\)", ", type='global'\\)", newick)
    # msa1 <- with(x, eval(parse(text = newick)))
  }else if(nsq == 2){
    if(!quiet) cat("Aligning seed sequences\n")
    msa1 <- align(x[[1]], x[[2]], residues = residues, gapchar = gapchar,  ... = ...)
    seedweights <- c(1, 1)
  }else if(nsq == 1){
    msa1 <- matrix(x[[1]], nrow = 1)
    seedweights <- 1
  }else stop("Empty list")
  if(!quiet) cat("Deriving profile hidden Markov model\n")
  omniphmm <- derive.PHMM.default(msa1, seqweights = seedweights, wfactor = wfactor, k = k, residues = residues,
                          gapchar = gapchar, pseudocounts = pseudocounts, logspace = logspace,
                          qa = qa, qe = qe, DI = DI, ID = ID, omit.endgaps = omit.endgaps,
                          inserts = inserts, lambda = lambda, threshold = threshold,
                          name = name, description = description, compo = compo, consensus = consensus,
                          cpp = cpp, quiet = quiet)
  if(nsq < 3){
    if(!quiet) cat("Done\n")
    return(omniphmm)
  }
  if(is.null(refine)) refine <- "none"
  #refine <- toupper(refine)
  if(refine %in% c("Viterbi", "BaumWelch")){
    if(identical(seqweights, "Gerstein")){
      if(!quiet) cat("Calculating sequence weights\n")
      if(identical(sort(seeds), seq_along(x))){
        seqweights <- seedweights * wfactor
      }else{
        qds <- kdistance(x, k = k, alpha = if(AA) "Dayhoff6" else if(DNA) NULL else residues)
        guidetree <- as.dendrogram(hclust(qds, method = "average"))
        seqweights <- weight.dendrogram(guidetree, method = "Gerstein")[names(x)] * wfactor
      }
    }else if(is.null(seqweights)){
      seqweights <- rep(wfactor, nsq)
    }else{
      if(length(seqweights) != nsq) stop("invalid seqweights argument")
      seqweights <- seqweights * wfactor
    }
    if(length(seqweights) != nsq) stop("invalid seqweights argument")
    if(!quiet) cat("Refining model\n")
    finalphmm <- train(omniphmm, x, method = refine, seqweights = seqweights,
                       maxiter = maxiter, pseudocounts = pseudocounts,
                       inserts = inserts, lambda = lambda, threshold = threshold, deltaLL = deltaLL,
                       quiet = quiet, ... = ...)
  }else if (refine == "none"){
    finalphmm <- omniphmm
  }else stop("Argument 'refine' must be set to either 'Viterbi', 'BaumWelch' or 'none'.")
  if(!quiet) cat("Done\n")
  return(finalphmm)
}


#'@rdname derive.PHMM
#'@export
#'
derive.PHMM.default <- function(x, seqweights = "Gerstein", wfactor = 1, k = 5, residues = NULL,
                       gapchar = "-", endchar = "?", pseudocounts = "background",
                       logspace = TRUE, qa = NULL, qe = NULL,
                       inserts = "map", lambda = 0, threshold = 0.5,
                       DI = FALSE, ID = FALSE, omit.endgaps = TRUE,
                       name = NULL, description = NULL,
                       compo = FALSE, consensus = FALSE, cpp = TRUE, quiet = FALSE){
  if(!(is.matrix(x))) stop("invalid object type, x must be a matrix")
  DNA <- is.DNA(x) # raw DNA bytes
  AA <- is.AA(x) # raw AA bytes
  gapchar <- if(DNA) as.raw(4) else if(AA) as.raw(45) else gapchar
  endchar <- if(DNA) as.raw(2) else if(AA) as.raw(63) else endchar
  if(omit.endgaps) x <- trim(x, gapchar = gapchar, endchar = endchar, DNA = DNA, AA = AA)
  residues <- alphadetect(x, residues = residues, gapchar = gapchar, endchar = endchar)
  if(is.null(name)) name <- unname(sapply(match.call()[2], deparse))
  nres <- length(residues)
  n <- nrow(x)
  m <- ncol(x)
  states <- c("D", "M", "I")
  transitions <- c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  if(is.null(seqweights)){
    seqweights <- rep(wfactor, n)
  }else if(identical(seqweights, "Gerstein")){
    if(n > 2){
      xlist <- unalign(x, gapchar = gapchar)
      qds <- kdistance(xlist, k = k, alpha = if(AA) "Dayhoff6" else if(DNA) NULL else residues)
      guidetree <- as.dendrogram(hclust(qds, method = "average"))
      seqweights <- weight(guidetree, method = "Gerstein")[names(xlist)] * wfactor
    }else{
      seqweights <- rep(wfactor, n)
    }
  }else{
    if(length(seqweights) != n) stop("invalid seqweights argument")
    seqweights <- seqweights * wfactor
  }
  # background emission probabilities (qe)
  if(is.null(qe)){
    allecs <- if(AA){
      apply(x, 2, tabulate.AA, ambiguities = TRUE, seqweights = seqweights)
    }else if(DNA){
      apply(x, 2, tabulate.DNA, ambiguities = TRUE, seqweights = seqweights)
    }else{
      apply(x, 2, tabulate.char, residues = residues, seqweights = seqweights)
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
  #ends <- x == endchar
  gapweights <- gaps * seqweights
  if(identical(inserts, "none")){
    inserts <- rep(FALSE, m)
  }else if(identical(inserts, "inherited")){
    inserts <- attr(x, "inserts")
    if(is.null(inserts)) stop("No inserts available to inherit from alignment")
  }else if(identical(inserts, "threshold")){
    inserts <- apply(gapweights, 2, sum) > threshold * n
  }else if(identical(inserts, "map")){
    if(n < 5){
      # if(!quiet) cat("Maximum a priori insert assignment unsuitable for fewer than five sequences.
      #                Switching to threshold method\n")
      inserts <- apply(gapweights, 2, sum) > threshold * n
    }else{
      inserts <- !map(x, seqweights = seqweights, residues = residues,
                      gapchar = gapchar, endchar = endchar, pseudocounts = pseudocounts,
                      qa = qa, qe = qe)
      if(sum(!inserts) < 3) {
        if(!quiet) cat("Maximum a priori insert assignment produced model with fewer than three modules.
                       Switching to threshold method\n")
        inserts <- apply(gapweights, 2, sum) > threshold * n
      }
    }
  }else if(!(mode(inserts) == "logical" & length(inserts) == ncol(x))){
    stop("invalid inserts argument")
  }
  l <- sum(!inserts) # PHMM length (excluding B & E positions)
  # emission counts
  ecs <- if(AA){
    apply(x[, !inserts, drop = F], 2, tabulate.AA,
          ambiguities = TRUE, seqweights = seqweights)
  }else if(DNA){
    apply(x[, !inserts, drop = F], 2, tabulate.DNA,
          ambiguities = TRUE, seqweights = seqweights)
  }else{
    apply(x[, !inserts, drop = F], 2, tabulate.char,
          residues = residues, seqweights = seqweights)
  }
  if(length(ecs) > 0){
    dimnames(ecs) <- list(residue = residues, position = 1:l)
  }else ecs = NULL

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
  # convert pseudocounts arg to list of vectors if necessary
  if(is.list(pseudocounts)){
    if(!all(c("A", "E") %in% names(pseudocounts))){
      stop("pseudocounts list is missing
           named vectors 'A' (transition pseudocounts) and/or 'E'
           (emission pseudocounts)")
    }
  }else if(identical(pseudocounts, "background")){
    pseudocounts <- list(A = qa * (7 + sum(c(DI, ID))), E = qe * nres)
  }else if(identical(pseudocounts, "Laplace")){
    pseudocounts <- list(A = c(1,1,DI,1,1,1,ID,1,1), E = rep(1, nres))
  }else if(identical(pseudocounts, "none")){
    pseudocounts <- list(A = rep(0, 9), E = rep(0, nres))
  }else stop("invalid pseudocounts argument")
  tcs <- tcs + pseudocounts$A
  tcs[1:3, 1] <- tcs[c(1, 4, 7), l + 1] <- 0
  if(!DI) tcs[3, ] <- 0
  if(!ID) tcs[7, ] <- 0
  if(is.null(ecs)){
    E <- matrix(nrow = nres, ncol = 0)
    rownames(E) <- residues
  }else{
    ecs <- ecs + pseudocounts$E
    E <- t(t(ecs)/apply(ecs, 2, sum))
  }
  A <- t(tcs)
  for(i in c(1, 4, 7)) {
    A[, i:(i + 2)] <- A[, i:(i + 2)]/apply(A[, i:(i + 2), drop = F], 1, sum)
  }
  A[1, 1:3] <- 0 # gets rid of NaNs caused by division by zero
  A <- t(A)
  inslens <- insertlengths(!inserts)
  #which alignment columns correspond to which model positions?
  alignment <- which(!inserts)
  if(length(alignment) > 0) names(alignment) <- 1:l
  if(logspace){
    A <- log(A)
    E <- log(E)
    qa <- log(qa)
    qe <- log(qe)
  }
  reference <- character(m)
  reference[inserts] <- "."
  reference[!inserts] <- "x"
  mask <- character(m)
  mask[inserts] <- "m"
  mask[!inserts] <- "."
  alphabet <- if(identical(toupper(sort(residues)), c("A", "C", "G", "T"))){
    "dna"
  }else if(identical(toupper(sort(residues)), c("A", "C", "G", "U"))){
    "rna"
  }else if(identical(toupper(sort(residues)), LETTERS[-c(2, 10, 15, 21, 24, 26)])){
    "amino"
  } else "custom"
  res <- structure(list(name = name, description = description,
                        size = l, alphabet = alphabet,
                        A = A, E = E, qa = qa, qe = qe, inserts = inserts,
                        insertlengths = inslens, alignment = alignment,
                        date = date(), nseq = n, weights = seqweights,
                        reference = reference, mask = mask),
                   class = "PHMM")
  # if(consensus){
  #   res$consensus <- generate(res, size = l * 100, random = FALSE, gapchar = ".")
  # }
  if(compo) res$compo <- log(apply(exp(E), 1, mean))
  return(res)
}

#' Optimized profile HMM construction.
#'
#' Assigns columns to insert states using the maximum \emph{a posteriori} method
#' outlined in Durbin et al. (1998).
#'
#' @param x a matrix of aligned sequences.
#' @param seqweights either NULL (default; all sequences are given a weight of 1)
#' or a numeric vector of sequence weights used to
#' derive the model. The sum of these weights should be equal to the number of
#' sequences in the alignment, so that mean(seqweights) = 1.
#' @param residues either NULL (default; emitted residues are automatically
#' detected from the alignment), or a case sensitive character vector specifying the
#' residue alphabet (e.g. A, C, G, T for DNA).
#' Note that the former option can be slow for large character matrices;
#' therefore specifying the residue alphabet can increase speed in these cases.
#' Also note that the default setting \code{residues = NULL} will not
#' detect rare residues that are not present in the alignment, and thus will
#' not assign them emission probabilities.
#' @param gapchar the character used to represent gaps in the alignment matrix.
#' Not required for \code{"DNAbin"} or \code{"AAbin"} objects.
#' @param pseudocounts the method used to account for the possible absence of certain
#' character states in the alignment. Currently only \code{"Laplace"},
#'\code{"background"} and \code{"none"} are supported. If \code{pseudocounts = "background"} and
#' some transition/emission types are absent from the alignment, Laplacean
#' pseudocounts are automatically added.
#' @param lambda penalty parameter used to favour models with fewer match states. Equivalent
#' to the log of the prior probability of marking each column (Durbin et al. 1998, pg 124).
#' @return a logical vector of length = ncol(x) indicating the columns to be assigned
#' as match states (\code{TRUE}) and those assigned as inserts (\code{FALSE}).
#' @return a logical vector indicating the columns of the input matrix that were allocated as
#' insert states
#' @references Durbin...
#' @export
#'
map <- function(x, seqweights = NULL, residues = NULL,
                gapchar = "-", endchar = "?", pseudocounts = "background",
                lambda = 0, qa = NULL, qe = NULL, cpp = TRUE){
  if(!is.matrix(x)) stop("x must be a matrix")
  L <- ncol(x)
  n <- nrow(x)
  if(n < 4) return(setNames(rep(TRUE, L), 1:L))
  AA <- is.AA(x)
  DNA <- is.DNA(x)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  endchar <- if(DNA) as.raw(2) else if(AA) as.raw(63) else endchar
  residues <- alphadetect(x, residues = residues, gapchar = gapchar, endchar = endchar)
  nres <- length(residues)
  transitions = c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  S <- sigma <- c(0, rep(NA, L + 1))
  if(is.null(seqweights)){
    seqweights <- rep(1, n)
  }
  if(length(seqweights) != n) stop("invalid seqweights argument")
  if(AA){
    ecs <- apply(x, 2, tabulate.AA, ambiguities = TRUE, seqweights = seqweights)
  }else if(DNA){
    ecs <- apply(x, 2, tabulate.DNA, ambiguities = TRUE, seqweights = seqweights)
  }else{
    ecs <- apply(x, 2, tabulate.char, residues = residues, seqweights = seqweights)
  }
  # ecs <- t(t(ecs) * seqweights)
  allecs <- apply(ecs, 1, sum)
  if(is.null(qe)) {
    qe <- log((allecs + 1)/sum(allecs + 1))
  }else if(all(qe >= 0 & qe <= 1) & round(sum(qe), 2) == 1){
    qe <- log(qe)
  }else if(any(qe > 0) | round(sum(exp(qe)), 2) != 1) stop("invalid qe")
  gaps <- x == gapchar
  #ends <- x == endchar
  notgaps <- cbind(TRUE, !gaps, TRUE)
  #if(!DI) alltcs[c(3, 7)] <- 0 ### need to work out DI strategy
  if(is.null(qa)){
    gapweights <- gaps * seqweights
    inserts <- apply(gapweights, 2, sum) > 0.5 * nrow(x)
    xtr <- matrix(nrow = nrow(x), ncol = ncol(x))
    insertsn <- matrix(rep(inserts, n), nrow = n, byrow = TRUE)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- tab9C(xtr, seqweights) # modules = sum(!inserts) + 2)
    alltcs <- apply(tcs, 1, sum)
    qa <- log((alltcs + 1)/sum(alltcs + 1)) # force addition of Laplace pseudos
  }else if(all(qa >= 0 & qa <= 1) & round(sum(qa), 2) == 1){
    qa <- log(qa)
  }else if(any(qa > 0) | round(sum(exp(qa)), 2) != 1) stop("invalid qa")
  # calculate Mj for j = 1 ... L + 1
  # convert pseudocounts arg to list of vectors if necessary
  if(is.list(pseudocounts)){
    if(!names(pseudocounts)[1] == "A"){
      stop("first element of pseudocounts list should be named 'A' and
           be a vector of length 9 (corresponding to pseudocount values
           for DD, DM, DI, MD, MM, MI, ID, IM, II transitions, respectively)")
    }
    if(!names(pseudocounts)[2] == "E"){
      stop("second element of pseudocounts list should be named 'E' and
           be a vector the same length as the size of the residue alphabet
          (e.g. 4 for DNA, 20 for amino acids)")
    }
  }else if(identical(pseudocounts, "background")){
    pseudocounts <- list(A = exp(qa) * 9, E = exp(qe) * nres)
  }else if(identical(pseudocounts, "Laplace")){
    pseudocounts <- list(A = rep(1, 9), E = rep(1, nres))
  }else if(identical(pseudocounts, "none")){
    pseudocounts <- list(A = rep(0, 9), E = rep(0, nres))
  }else stop("invalid pseudocounts argument")
  if(cpp){
    res <- mapC(ecs, notgaps, pseudocounts, seqweights, qe, lambda)
  }else{
    ecs2 <- ecs + pseudocounts$E
    term2 <- t(t(ecs2)/apply(ecs2, 2, sum))
    term2[ecs != 0] <- log(term2[ecs != 0]) # increase speed for conserved alignments
    M <- apply(ecs * term2, 2, sum)
    M <- c(0, M, 0)
    ecsj <- structure(numeric(nres), names = residues)
    tcsij <- structure(numeric(9), names = transitions)
    icsj <- structure(numeric(n), names = rownames(x)) #insert counts
    #alphaxy <- if(pseudocounts == "background") exp(qa) * 9 else rep(1, 9)
    for(j in 2:(L + 2)){
      tau <- iota <- numeric(j - 1)
      ecsij <- ecsj
      icsij <- icsj
      for(i in 1:(j - 1)){
        if(i < j - 1){
          iota[i] <- sum(ecsij * qe)
          ecsij <- ecsij - ecs[, i]
          zeroinserts <- icsij < 0.00001
          #need a way to avoid counting NAs in notgaps matrix (endchars)
          if(any(zeroinserts)){
            tcsij[1] <- sum(seqweights[!notgaps[, i] & !notgaps[, j] & zeroinserts]) #DD
            tcsij[2] <- sum(seqweights[!notgaps[, i] & notgaps[, j] & zeroinserts]) #DM
            tcsij[4] <- sum(seqweights[notgaps[, i] & !notgaps[, j] & zeroinserts]) #MD
            tcsij[5] <- sum(seqweights[notgaps[, i] & notgaps[, j] & zeroinserts]) #MM
            tcsij[3] <- sum(seqweights[!notgaps[, i] & !zeroinserts]) #DI
            tcsij[6] <- sum(seqweights[notgaps[, i] & !zeroinserts]) #MI
            tcsij[7] <- sum(seqweights[!notgaps[, j] & !zeroinserts]) #ID
            tcsij[8] <- sum(seqweights[notgaps[, j] & !zeroinserts]) #IM
          }else{
            tcsij[c(1, 2, 4, 5)] <- 0
            tcsij[3] <- sum(seqweights[!notgaps[, i]])
            tcsij[6] <- sum(seqweights[notgaps[, i]])
            tcsij[7] <- sum(seqweights[!notgaps[, j]])
            tcsij[8] <- sum(seqweights[notgaps[, j]])
          }
          tcsij[9] <- sum(icsij) - (tcsij[3] + tcsij[6]) #II
          icsij <- icsij - notgaps[, i + 1] * seqweights
          # ends up as vec of zeros at end of each i cycle
        }else{
          tcsij[1] <- sum(seqweights[!notgaps[, i] & !notgaps[, j]]) #DD
          tcsij[2] <- sum(seqweights[!notgaps[, i] & notgaps[, j]]) #DM
          tcsij[3] <- 0 #DI
          tcsij[4] <- sum(seqweights[notgaps[, i] & !notgaps[, j]]) #MD
          tcsij[5] <- sum(seqweights[notgaps[, i] & notgaps[, j]]) #MM
          tcsij[6] <- 0 #MI
          tcsij[7] <- 0 #ID
          tcsij[8] <- 0 #IM
          tcsij[9] <- 0 #II
        }
        cxy <- tcsij + pseudocounts$A
        axy <- cxy/c(rep(sum(cxy[1:3]), 3), rep(sum(cxy[4:6]), 3), rep(sum(cxy[7:9]), 3))
        tau[i] <- sum(cxy * log(axy))
      }
      if(j < L + 2){
        ecsj <- ecsj + ecs[, j - 1]
        icsj <- icsj + seqweights * notgaps[, j]
      }
      tmp <- S[1:(j - 1)] + tau + iota + M[j] + lambda
      sigma[j] <- whichmax(tmp)
      S[j] <- tmp[sigma[j]]
    }
    res <- structure(logical(L + 2), names = 0:(L + 1))
    res[L + 2] <- TRUE
    j <- sigma[L + 2]
    while(j > 0){
      res[j] <- TRUE
      j <- sigma[j]
    }
    res <- res[-(c(1, L + 2))]
  }
  return(res)
}

