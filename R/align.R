#' Multiple sequence alignment in R.
#'
#' \code{align} performs a multiple alignment on a list of
#'   sequences using profile hidden Markov models.
#'
#' @param sequences a list of DNA, amino acid, or other character sequences
#'   consisting of symbols emitted from the chosen residue alphabet.
#'   The vectors can either be of mode "raw" (consistent with the "DNAbin"
#'   or "AAbin" coding scheme set out in the \code{\link[ape]{ape}} package),
#'   or "character", in which case the alphabet should be ideally be specified in
#'   the \code{residues} argument.
#' @param model an optional profile hidden Markov model (a \code{"PHMM"}
#'   object) to align the sequences to. If \code{NULL} a PHMM will
#'   be derived from the list of sequences, and each sequence
#'   will be aligned back to the model to produce the multiple sequence
#'   alignment.
#' @param seqweights either NULL (default; all sequences are given weights
#'   of 1), a numeric vector the same length as \code{x} representing
#'   the sequence weights used to derive the model, or a character string giving
#'   the method to derive the weights from the sequences. Currently only the
#'   \code{"Gerstein"} method is supported (default). For this method, a
#'   tree is first created by k-mer counting (see \code{\link[phylogram]{topdown}}),
#'   and sequence weights are then derived from the tree using the 'bottom up'
#'   algorithm of Gerstein et al. (1994).
#' @param refine the method used to iteratively refine the model parameters
#'   following the initial progressive alignment and model derivation step.
#'   Current supported options are \code{"Viterbi"} (Viterbi training;
#'   the default option), \code{"BaumWelch"} (a modified version of the
#'   Expectation-Maximization algorithm), and "none" (skips the model
#'   refinement step).
#' @param k integer representing the k-mer size to be used for calculating
#'   the distance matrix used in tree-based sequence weighting. Defaults to
#'   5. Note that high values oof k (> 8) may be slow to compute and use a lot
#'   of memory due to the large numbers of calculations required.
#' @param maxiter the maximum number of EM iterations or Viterbi training
#'   iterations to carry out before the cycling process is terminated and
#'   the partially trained model is returned. Defaults to 100.
#' @param maxsize integer giving the upper bound on the number of modules
#'   in the PHMM. If NULL no maximum size is enforced.
#' @param inserts character string giving the model construction method
#'   in which alignment columns
#'   are marked as either match or insert states. Accepted methods include
#'   \code{"threshold"} (only columns with fewer than a specified
#'   proportion of gaps form match states in the model), \code{"map"} (default;
#'   match and insert columns are found using the maximum \emph{a posteriori}
#'   method outlined in Durbin et al. (1998) chapter 5.7), \code{"inherited"}
#'   (match and insert columns are inherited from the "inserts" attribute
#'   of the input alignment), and \code{"none"} (all columns are assigned
#'   match states in the model). Alternatively, insert columns can be
#'   specified manually by providing a logical vector the same length
#'   as the number of columns in the alignment, with \code{TRUE} for insert
#'   columns and \code{FALSE} for match states.
#' @param lambda penalty parameter used to favour models with fewer match
#'   states. Equivalent to the log of the prior probability of marking each
#'   column (Durbin et al. 1998, pg 124). Only applicable when
#'   \code{inserts = "map"}.
#' @param threshold the maximum proportion of gaps for an alignment column
#'   to be considered for a match state in the PHMM (defaults to 0.5).
#'   Only applicable when \code{inserts = "threshold"}.
#'   Note that the maximum \emph{a posteriori}
#'   method works poorly for alignments with few sequences,
#'   so the 'threshold' method is
#'   automatically used when the number of sequences is less than 5.
#' @param deltaLL numeric, the maximum change in log likelihood between EM
#'   iterations before the cycling procedure is terminated (signifying model
#'   convergence). Defaults to 1E-07. Only applicable if
#'   \code{method = "BaumWelch"}.
#' @param DI logical indicating whether delete-insert transitions should be
#'   allowed in the profile hidden Markov model (if applicable). Defaults
#'   to FALSE.
#' @param ID logical indicating whether insert-delete transitions should be
#'   allowed in the profile hidden Markov model (if applicable). Defaults
#'   to FALSE.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Furthermore, the default setting
#'   \code{residues = NULL} will not detect rare residues that are not present
#'   in the sequences, and thus will not assign them emission probabilities
#'   in the model. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix.
#'   Ignored for \code{"DNAbin"} or \code{"AAbin"} objects. Defaults to "-"
#'   otherwise.
#' @param pseudocounts character string, either "background", Laplace"
#'   or "none". Used to account for the possible absence of certain
#'   transition and/or emission types in the input sequences.
#'   If \code{pseudocounts = "background"} (default), pseudocounts
#'   are calculated from the background transition and emission
#'   frequencies in the sequences.
#'   If \code{pseudocounts = "Laplace"} one of each possible transition
#'   and emission type is added to the transition and emission counts.
#'   If \code{pseudocounts = "none"} no pseudocounts are added (not
#'   generally recommended, since low frequency transition/emission types
#'   may be excluded from the model).
#'   Alternatively this argument can be a two-element list containing
#'   a matrix of transition pseudocounts
#'   as its first element and a matrix of emission pseudocounts as its
#'   second.
#' @param qa an optional named 9-element vector of background transition
#'   probabilities with \code{dimnames(qa) = c("DD", "DM", "DI", "MD", "MM",
#'   "MI", "ID", "IM", "II")}, where M, I and D represent match, insert and
#'   delete states, respectively. If \code{NULL}, background transition
#'   probabilities are estimated from the sequences.
#' @param qe an optional named vector of background emission probabilities
#'   the same length as the residue alphabet (i.e. 4 for nucleotides and 20
#'   for amino acids) and with corresponding names (i.e. \code{c("A", "T",
#'   "G", "C")} for DNA). If \code{qe = NULL}, background emission probabilities
#'   are automatically derived from the sequences.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if 'sequences' is not a list.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the number of cores
#'   used is one less than the total number of cores available.
#'   Note that in this case there
#'   may be a tradeoff in terms of speed depending on the number and size
#'   of sequences to be aligned, due to the extra time required to initialize
#'   the cluster.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @param ... aditional arguments to be passed to \code{"Viterbi"} (if
#'   \code{refine = "Viterbi"}) or \code{"forward"} (if
#'   \code{refine = "BaumWelch"}).
#' @return a matrix of aligned sequences, with the same mode and class as the
#'   input sequence list.
#' @details This function uses a progressive alignment procedure similar to
#'   the Clustal Omega algorithm (Sievers et al. 2011). The procedure
#'   involves an initial progressive multiple sequence alignment via a guide tree,
#'   followed by the derivation of a profile hidden Markov model
#'   from the alignment, an
#'   iterative model refinement step, and finally the alignment of the
#'   sequences back to the model. If only two sequences are provided,
#'   a standard pairwise alignment is carried out using the Needleman-Wunch
#'   or Smith-Waterman algorithm.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H,
#'   Remmert M, Soding J, Thompson JD, Higgins DG (2011) Fast, scalable generation
#'   of high-quality protein multiple sequence alignments using Clustal Omega.
#'   \emph{Molecular Systems Biology}, \strong{7}, 539.
#'
#' @seealso \code{\link{unalign}}
#' @examples
#'   ## Protein pairwise alignment example from Durbin et al. (1998) chapter 2.
#'   x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#'   y <- c("P", "A", "W", "H", "E", "A", "E")
#'   sequences <- list(x = x, y = y)
#'   glo <- align(sequences, type = "global")
#'   sem <- align(sequences, type = "semiglobal")
#'   loc <- align(sequences, type = "local")
#'   glo
#'   sem
#'   loc
#' @name align
################################################################################
# align <- function(sequences, model = NULL, seqweights = "Gerstein", k = 5,
#                   refine = "Viterbi", maxiter = 100, maxsize = NULL,
#                   inserts = "map", lambda = 0, threshold = 0.5, deltaLL = 1E-07,
#                   DI = FALSE, ID = FALSE, residues = NULL, gap = "-",
#                   pseudocounts = "background", qa = NULL, qe = NULL,
#                   quiet = FALSE, ...){
#   UseMethod("align")
# }
align <- function(sequences, ...){
  UseMethod("align")
}

#' @rdname align
################################################################################
align.DNAbin <- function(sequences, model = NULL, seqweights = "Gerstein",
                         refine = "Viterbi", k = 5, maxiter = 100,
                         maxsize = NULL, inserts = "map", lambda = 0,
                         threshold = 0.5, deltaLL = 1E-07, DI = FALSE,
                         ID = FALSE, residues = NULL, gap = "-",
                         pseudocounts = "background", qa = NULL, qe = NULL,
                         cores = 1, quiet = FALSE, ...){
  if(is.list(sequences)){
    align.list(sequences, model = model, seqweights = seqweights, refine = refine, k = k,
               maxiter = maxiter, maxsize = maxsize, inserts = inserts, lambda = lambda,
               threshold = threshold, deltaLL = deltaLL, DI = DI, ID = ID,
               residues = residues, gap = gap, pseudocounts = pseudocounts,
               qa = qa, qe = qe, cores = cores, quiet = quiet, ... = ...)
  }else{
    align.default(sequences, model = model, residues = residues, gap = gap,
                  pseudocounts = pseudocounts, maxsize = maxsize, quiet = quiet, ... = ...)
  }
}

#' @rdname align
################################################################################
align.AAbin <- function(sequences, model = NULL, seqweights = "Gerstein",
                        refine = "Viterbi", k = 5, maxiter = 100,
                        maxsize = NULL, inserts = "map", lambda = 0,
                        threshold = 0.5, deltaLL = 1E-07, DI = FALSE,
                        ID = FALSE, residues = NULL, gap = "-",
                        pseudocounts = "background", qa = NULL, qe = NULL,
                        cores = 1, quiet = FALSE, ...){
  if(is.list(sequences)){
    align.list(sequences, model = model, seqweights = seqweights, refine = refine, k = k,
               maxiter = maxiter, maxsize = maxsize, inserts = inserts, lambda = lambda,
               threshold = threshold, deltaLL = deltaLL, DI = DI, ID = ID,
               residues = residues, gap = gap, pseudocounts = pseudocounts,
               qa = qa, qe = qe, cores = cores, quiet = quiet, ... = ...)
  }else{
    align.default(sequences, model = model, residues = residues, gap = gap,
                  pseudocounts = pseudocounts, maxsize = maxsize, quiet = quiet, ... = ...)
  }
}

#' @rdname align
################################################################################
align.list <- function(sequences, model = NULL, seqweights = "Gerstein", k = 5,
                       refine = "Viterbi", maxiter = 100,
                       maxsize = NULL, inserts = "map", lambda = 0,
                       threshold = 0.5, deltaLL = 1E-07,
                       DI = FALSE, ID = FALSE, residues = NULL,
                       gap = "-", pseudocounts = "background",
                       qa = NULL, qe = NULL, cores = 1, quiet = FALSE, ...){
  nseq <- length(sequences)
  DNA <- .isDNA(sequences)
  AA <- .isAA(sequences)
  if(DNA) class(sequences) <- "DNAbin" else if(AA) class(sequences) <- "AAbin"
  residues <- .alphadetect(sequences, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  for(i in 1:length(sequences)) sequences[[i]] <- sequences[[i]][sequences[[i]] != gap]
  if(is.null(model)){
    if(nseq == 2){
      alignment <- align.default(sequences[[1]], sequences[[2]], residues = residues,
                                 gap = gap, pseudocounts = pseudocounts,
                                 quiet = quiet, ... = ...)
      rownames(alignment) <- names(sequences)
      return(alignment)
    }else if(nseq == 1){
      alignment <- matrix(sequences[[1]], nrow = 1)
      rownames(alignment) <- names(sequences)
      return(alignment)
    }
    #if(!quiet) cat("Calculating pairwise distances\n")
    if(nseq > 100){
      nseeds <- ceiling(log(nseq, 2)^2) # LLR algorithm see Blacksheilds et al 2010
      seeds <- sample(1:nseq, size = nseeds)
    }else{
      seeds <- seq_along(sequences)
    }
    if(identical(seqweights, "Gerstein")){
      if(!quiet) cat("Calculating sequence weights\n")
      guidetree <- phylogram::topdown(sequences, k = k, residues = residues, gap = gap)
      myseqweights <- weight.dendrogram(guidetree, method = "Gerstein")[names(sequences)]
    }else if(is.null(seqweights)){
      myseqweights <- rep(1, nseq)
    }else{
      myseqweights <- seqweights
    }
    #if(!quiet) cat("Aligning seed sequences and deriving model\n")
    phmm <- derivePHMM.list(sequences, seeds = seeds, refine = refine, maxiter = maxiter,
                            seqweights = myseqweights, k = k, residues = residues, gap = gap,
                            maxsize = maxsize, inserts = inserts, lambda = lambda,
                            threshold = threshold, deltaLL = deltaLL, DI = DI, ID = ID,
                            pseudocounts = pseudocounts, logspace = TRUE, qa = qa, qe = qe,
                            quiet = quiet, ... = ...)
    if(!quiet) cat("Aligning sequences to model\n")
    res <- align.list(sequences, model = phmm, ... = ...)
    if(!quiet) cat("Produced alignment with", ncol(res), "columns (including inserts)\n")
    return(res)
  }else{
    #note changes here also need apply to 'train'
    stopifnot(inherits(model, "PHMM"))
    l <- model$size
    pathfinder <- function(s, model, ...){
      vit <- Viterbi(model, s, ... = ...)
      res <- c(vit$path, 1) #append the final transition to end state
      # this is just so the c++ function knows where to stop
      attr(res, "score") <- vit$score
      return(res)
    }
    ### insert parLapply code here
    #paths <- lapply(sequences, pathfinder, model, ...)
    if(inherits(cores, "cluster")){
      paths <- parallel::parLapply(cores, sequences, pathfinder, model = model, ...)
    }else if(cores == 1){
      paths <- lapply(sequences, pathfinder, model, ...)
    }else{
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(cores > navailcores) stop("Number of cores is more than number available")
      cl <- parallel::makeCluster(cores)
      paths <- parallel::parLapply(cl, sequences, pathfinder, model = model, ...)
      parallel::stopCluster(cl)
    }
    ###
    fragseqs <- mapply(if(DNA | AA) .fragR else .fragC, sequences, paths, l = l,
                       gap = gap, SIMPLIFY = FALSE)
    odds <- seq(1, 2 * l + 1, by = 2)
    evens <- seq(2, 2 * l, by = 2)
    length(fragseqs)
    inslens <- lapply(fragseqs, function(e) sapply(e[odds], length))
    inslens <- matrix(unlist(inslens), nrow = nseq, byrow = TRUE)
    insmaxs <- apply(inslens, 2, max)
    insappends <- t(insmaxs - t(inslens))
    for(i in 1:nseq){
      needsapp <- insappends[i, ] > 0
      if(any(needsapp)){
        apps <- lapply(insappends[i, needsapp], function(e) rep(gap, e))
        fragseqs[[i]][odds][needsapp] <- mapply(c, fragseqs[[i]][odds][needsapp],
                                                apps, SIMPLIFY = FALSE)
      }
    }
    unfragseqs <- lapply(fragseqs, unlist)
    res <- matrix(unlist(unfragseqs), nrow = nseq, byrow = TRUE)
    score <- sum(sapply(paths, function(p) attr(p, "score")))
    inserts <- vector(length = 2 * l + 1, mode = "list")
    inserts[evens] <- FALSE
    inserts[odds] <- lapply(insmaxs, function(e) rep(TRUE, e))
    inserts <- unlist(inserts)
    resnames <- vector(length = 2 * l + 1, mode = "list")
    resnames[evens] <- paste(1:l)
    resnames[odds] <- lapply(insmaxs, function(e) rep("I", e))
    resnames <- unlist(resnames)
    colnames(res) <- resnames
    rownames(res) <- names(sequences)
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    attr(res, "score") <- score
    attr(res, "inserts") <- inserts
    return(res)
  }
}
################################################################################
#' @rdname align
################################################################################
align.default <- function(sequences, model, pseudocounts = "background",
                          residues = NULL, gap = "-", maxsize = NULL,
                          quiet = FALSE, ...){
  if(is.null(model)) return(sequences)
  DNA <- .isDNA(sequences)
  AA <- .isAA(sequences)
  if(DNA){
    if(!.isDNA(model)) stop("class(sequences) and class(model) must match")
    gap <- as.raw(4)
    # changes here need also apply in Viterbi.default and Viterbi.PHMM
    if(is.list(sequences)){
      if(length(sequences) == 1){
        if(!is.matrix(sequences)) sequences <- matrix(sequences[[1]], nrow = 1,
                                                      dimnames = list(names(sequences), NULL))
        class(sequences) <- "DNAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(sequences)){
        sequences <- matrix(sequences, nrow = 1, dimnames = list(deparse(substitute(sequences)), NULL))
      }
      class(sequences) <- "DNAbin"
    }
    if(is.list(model)){
      if(length(model) == 1){
        if(!is.matrix(model)) model <- matrix(model[[1]], nrow = 1, dimnames = list(names(model), NULL))
        class(model) <- "DNAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(model)) {
        model <- matrix(model, nrow = 1, dimnames = list(deparse(substitute(model)), NULL))
      }
      class(model) <- "DNAbin"
    }
  }else if(AA){
    if(!.isAA(model)) stop("class(sequences) and class(model) must match")
    gap <- as.raw(45)
    # changes here need also apply in Viterbi.default and Viterbi.PHMM
    if(is.list(sequences)){
      if(length(sequences) == 1){
        if(!is.matrix(sequences)) sequences <- matrix(sequences[[1]], nrow = 1, dimnames = list(names(sequences), NULL))
        class(sequences) <- "AAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(sequences)) sequences <- matrix(sequences, nrow = 1, dimnames = list(deparse(substitute(sequences)), NULL))
      class(sequences) <- "AAbin"
    }
    if(is.list(model)){
      if(length(model) == 1){
        if(!is.matrix(model)) model <- matrix(model[[1]], nrow = 1, dimnames = list(names(model), NULL))
        class(model) <- "AAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(model)) model <- matrix(model, nrow = 1, dimnames = list(
        deparse(substitute(model)), NULL))
      class(model) <- "AAbin"
    }
  }else{
    if(!is.matrix(sequences)){
      sequences <- matrix(sequences, nrow = 1, dimnames = list(
        deparse(substitute(sequences)), NULL))
      #rownames(sequences) <-
    }
    if(!is.matrix(model)){
      model <- matrix(model, nrow = 1, dimnames = list(
        deparse(substitute(model)), NULL))
    }
  }
  if(nrow(sequences) == 1 & nrow(model) == 1){
    alig <- Viterbi(sequences, model, ... = ...)
    #, offset = offset) ###not necessary for vec vs vec
    ## TODO local alignment
    if(!any(alig$path == 1)){
      warning("No local alignment found, returning NULL")
      return(NULL)
    }
    xind <- yind <- alig$path
    #xind[alig$path != 2] <- 1:length(sequences)
    xind[alig$path != 2] <- seq(alig$start[1], length.out = sum(alig$path != 2))
    xind[alig$path == 2] <- 0
    newx <- c(gap, as.vector(sequences))[xind + 1]
    # yind[alig$path != 0] <- 1:length(model)
    yind[alig$path != 0] <- seq(alig$start[2], length.out = sum(alig$path != 0))
    yind[alig$path == 0] <- 0
    newy <- c(gap, as.vector(model))[yind + 1]
    res <- rbind(newx, newy)
    rownames(res) <- c(rownames(sequences), rownames(model))
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    return(res)
  }else if(sum(c(nrow(sequences) == 1, nrow(model) == 1)) == 1){
    if(nrow(sequences) == 1){
      tmp <- sequences
      sequences <- model
      model <- tmp
      rm(tmp) # the old switcharoo
    }
    n <- nrow(sequences)
    z <- derivePHMM(sequences, seqweights = "Gerstein", pseudocounts = pseudocounts,
                     residues = residues, logspace = TRUE)
    l <- z$size
    alignment <- Viterbi(z, model, ... = ...) ### changed from alig
    path <- alignment$path
    # lay x out as list with insert elements
    newrow <- vector(length = l * 2 + 1, mode = "list")
    odds <- seq(from = 1, to = length(newrow), by = 2) #insert columns
    evens <- seq(from = 2, to = length(newrow), by = 2) # match columns
    newrow[evens] <- lapply(which(!z$inserts), function(e) e)
    if(any(z$inserts)){
      itp <- apply(rbind(c(FALSE, z$inserts), c(z$inserts, FALSE)), 2, .decimal, from = 2)
      ist <- which(itp == 1)
      ien <- which(itp == 2) - 1
      newrow[odds][which(itp[itp < 2] == 1)] <- mapply(":", ist, ien, SIMPLIFY = FALSE)
    }
    newrow <- lapply(newrow, function(e) if(is.null(e)) 0 else e)
    newx <- lapply(newrow, function(e) sequences[, e, drop = FALSE])
    # analogous list for y but witout insert elements
    newrow <- lapply(1:ncol(model), function(e) e)
    newy <- lapply(newrow, function(e) model[, e, drop = FALSE])
    #
    newxrows <- matrix(gap, nrow = n, ncol = ncol(sequences) + ncol(model))
    rownames(newxrows) <- rownames(sequences)
    newyrow <- matrix(gap, nrow = 1, ncol = ncol(sequences) + ncol(model))
    rownames(newyrow) = rownames(model)
    isinsert <- vector(mode = "logical", length = ncol(sequences) + ncol(model))
    if(DNA){
      class(newxrows) <- "DNAbin"
      class(newyrow) <- "DNAbin"
    }else if(AA){
      class(newxrows) <- "AAbin"
      class(newyrow) <- "AAbin"
    }
    position <- 1
    xcounter <- 2 * alignment$start[1]
    ycounter <- alignment$start[2]
    if(xcounter == 2 & ycounter == 1){ # global and semiglobal alignments
      rightshift <- ncol(newx[[1]])
      if(rightshift > 0){
        newxrows <- .insert(newx[[1]], into = newxrows, at = 1)
        position <- position + rightshift
        isinsert[1:rightshift] <- TRUE
      }
    }
    for(i in seq_along(path)){
      if(path[i] == 1){ #Match state
        rightshift <- 1 + ncol(newx[[xcounter + 1]])
        # match + insert
        newxrows <- .insert(newx[[xcounter]], into = newxrows, at = position)
        newyrow <- .insert(newy[[ycounter]], into = newyrow, at = position)
        position <- position + 1
        if(rightshift > 1){
          newxrows <- .insert(newx[[xcounter + 1]], into = newxrows, at = position)
          isinsert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
        ycounter <- ycounter + 1
      }else if(path[i] == 0){ #Delete state
        rightshift <- 1 + ncol(newx[[xcounter + 1]])
        newxrows <- .insert(newx[[xcounter]], into = newxrows, at = position)
        position <- position + 1
        if(rightshift > 0){
          newxrows <- .insert(newx[[xcounter + 1]], into = newxrows, at = position)
          isinsert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
      }else if(path[i] == 2){
        rightshift <- 1
        newyrow <- .insert(newy[[ycounter]], into = newyrow, at = position)
        position <- position + 1
        ycounter <- ycounter + 1
      }
    }
    position <- position - 1
    newxrows <- newxrows[, 1:position]
    newyrow <- newyrow[, 1:position]
    isinsert <- isinsert[1:position]
    res <- rbind(newxrows, newyrow)
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    # res.list <- unalign(res)
    # res.weights <- weight(res.list, method = "Gerstein", k = 5)
    # res.phmm <- derivePHMM.default(res, seqweights = res.weights, quiet = TRUE)
    # res.phmm <- train(res.phmm, res.list, method = "Viterbi", maxiter = 10,
    #                   logspace = TRUE, quiet = TRUE, ... = ...)
    # res <- align.list(sequences = res.list, model = res.phmm, quiet = TRUE, ... = ...)
    return(res)
  }else if(nrow(sequences) > 1 & nrow(model) > 1){
    #cat(rownames(sequences), "\n")
    #cat(rownames(model), "\n")
    # mat1 <<- sequences
    # mat2 <<- model
    nx <- nrow(sequences)
    ny <- nrow(model)
    zx <- derivePHMM(sequences, seqweights = "Gerstein", pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    zy <- derivePHMM(model,  seqweights = "Gerstein", pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    lx <- zx$size
    ly <- zy$size
    alignment <- Viterbi(zx, zy, ... = ...)
    path <- alignment$path
    # lay x out as list with insert elements
    newrow <- vector(length = lx * 2 + 1, mode = "list")
    odds <- seq(from = 1, to = length(newrow), by = 2) #insert columns
    evens <- seq(from = 2, to = length(newrow), by = 2) # match columns
    newrow[evens] <- lapply(which(!zx$inserts), function(e) e)
    if(any(zx$inserts)){
      itp <- apply(rbind(c(FALSE, zx$inserts), c(zx$inserts, FALSE)), 2, .decimal, from = 2)
      ist <- which(itp == 1)
      ien <- which(itp == 2) - 1
      newrow[odds][which(itp[itp < 2] == 1)] <- mapply(":", ist, ien, SIMPLIFY = FALSE)
    }
    newrow <- lapply(newrow, function(e) if(is.null(e)) 0 else e)
    newx <- lapply(newrow, function(e) sequences[, e, drop = FALSE])
    # lay y out as list with insert elements
    newrow <- vector(length = ly * 2 + 1, mode = "list")
    odds <- seq(from = 1, to = length(newrow), by = 2) #insert columns
    evens <- seq(from = 2, to = length(newrow), by = 2) # match columns
    newrow[evens] <- lapply(which(!zy$inserts), function(e) e)
    if(any(zy$inserts)){
      itp <- apply(rbind(c(FALSE, zy$inserts), c(zy$inserts, FALSE)), 2, .decimal, from = 2)
      ist <- which(itp == 1)
      ien <- which(itp == 2) - 1
      newrow[odds][which(itp[itp < 2] == 1)] <- mapply(":", ist, ien, SIMPLIFY = FALSE)
    }
    newrow <- lapply(newrow, function(e) if(is.null(e)) 0 else e)
    newy <- lapply(newrow, function(e) model[, e, drop = FALSE])
    #create output alignment
    newxrows <- matrix(gap, nrow = nx, ncol = ncol(sequences) + ncol(model))
    rownames(newxrows) <- rownames(sequences)
    newyrows <- matrix(gap, nrow = ny, ncol = ncol(sequences) + ncol(model))
    rownames(newyrows) <- rownames(model)
    isinsert <- vector(mode = "logical", length = ncol(sequences) + ncol(model))
    if(DNA){
      class(newxrows) <- "DNAbin"
      class(newyrows) <- "DNAbin"
    }else if(AA){
      class(newxrows) <- "AAbin"
      class(newyrows) <- "AAbin"
    }
    position <- 1
    xcounter <- 2 * alignment$start[1]
    ycounter <- 2 * alignment$start[2]
    if(xcounter == 2 & ycounter == 2){
      rightshift <- max(c(ncol(newx[[1]]), newy[[1]]))
      if(rightshift > 0){
        newxrows <- .insert(newx[[1]], into = newxrows, at = 1)
        newyrows <- .insert(newy[[1]], into = newyrows, at = 1)
        position <- position + rightshift
        isinsert[1:rightshift] <- TRUE
      }
    }
    for(i in seq_along(path)){
      if(path[i] == 2){ #MM
        rightshift <- 1 + max(c(ncol(newx[[xcounter + 1]]), ncol(newy[[ycounter + 1]])))
        # match + insert
        newxrows <- .insert(newx[[xcounter]], into = newxrows, at = position)
        newyrows <- .insert(newy[[ycounter]], into = newyrows, at = position)
        position <- position + 1
        if(rightshift > 1){
          newxrows <- .insert(newx[[xcounter + 1]], into = newxrows, at = position)
          newyrows <- .insert(newy[[ycounter + 1]], into = newyrows, at = position)
          isinsert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
        ycounter <- ycounter + 2
      }else if(path[i] < 2){
        rightshift <- 1 + ncol(newx[[xcounter + 1]])
        newxrows <- .insert(newx[[xcounter]], into = newxrows, at = position)
        position <- position + 1
        if(rightshift > 0){
          newxrows <- .insert(newx[[xcounter + 1]], into = newxrows, at = position)
          isinsert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
      }else if(path[i] > 2){
        rightshift <- 1 + ncol(newy[[ycounter + 1]])
        newyrows <- .insert(newy[[ycounter]], into = newyrows, at = position)
        position <- position + 1
        if(rightshift > 0){
          newyrows <- .insert(newy[[ycounter + 1]], into = newyrows, at = position)
          isinsert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        ycounter <- ycounter + 2
      }
    }
    position <- position - 1
    newxrows <- newxrows[, 1:position]
    newyrows <- newyrows[, 1:position]
    isinsert <- isinsert[1:position]
    res <- rbind(newxrows, newyrows)
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    res.list <- unalign(res)
    res.weights <- weight(res.list, method = "Gerstein", k = 5)
    # cat("\n", maxsize, "\n")
    # cat(sum(apply(res, 2, function(v) !any(v == gap))), "\n\n")
    newmaxsize <- if(is.null(maxsize)){
      NULL
    }else{
      max(c(sum(apply(res, 2, function(v) !any(v == gap))), maxsize))
    }
    res.phmm <- derivePHMM.default(res, seqweights = res.weights, maxsize = newmaxsize)
    res.phmm <- train(res.phmm, res.list, method = "Viterbi", maxiter = 3,
                      maxsize = newmaxsize, logspace = TRUE, quiet = TRUE, ... = ...)
    res <- align.list(sequences = res.list, model = res.phmm, ... = ...)
    return(res)
  }else{
    stop("invalid arguments provided for sequences and or y")
  }
}
################################################################################
#' Deconstruct an alignment.
#'
#' \code{unalign} deconstructs an alignment to a list of sequences.
#'
#' @param x a matrix of aligned sequences. Accepted modes are "character"
#'   and "raw" (for "DNAbin" and "AAbin" objects).
#' @inheritParams align
#' @return a list of sequences of the same mode and class as the input alignment
#'   (ie "DNAbin", "AAbin", or plain ASCII characters).
#' @details \code{unalign} works in the opposite way to \code{\link{align}},
#'   reducing a matrix of aligned sequences to a list of sequences without gaps.
#'   "DNAbin" and "AAbin" matrix objects are supported (and recommended for
#'   biological sequence data)
#' @author Shaun Wilkinson
#' @seealso \code{\link{align}}.
#' @examples
#' ## Convert the woodmouse alignment in the ape package to a list of
#' ## unaligned sequences
#' library(ape)
#' data(woodmouse)
#' sequences <- unalign(woodmouse)
################################################################################
unalign <- function(x, gap = "-"){
  #x is a matrix representing an alignment
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.list(x)){
    if(length(x) == 1){
      tmpname <- names(x)
      x <- x[[1]]
      if(is.null(dim(x))){
        x <- matrix(x, nrow = 1)
        rownames(x) <- tmpname
      }
    }
  }
  res <- vector(mode = "list", length = nrow(x))
  for(i in 1:nrow(x)) res[[i]] <- x[i, x[i, ] != gap, drop = TRUE]
  if(AA){
    res <- lapply(res, unclass)
    class(res) <- "AAbin"
  }else if(DNA){
    res <- lapply(res, unclass)
    class(res) <- "DNAbin"
  }
  names(res) <- rownames(x)
  return(res)
}
################################################################################
