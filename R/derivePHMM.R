#' Derive a profile hidden Markov model from sequences.
#'
#' \code{derivePHMM} generates a profile HMM from a given multiple sequence alignment
#'   or a list of unaligned sequences.
#'
#' @param x a matrix of aligned sequences or a list of unaligned sequences.
#'   Accepted modes are "character" and "raw" (for "DNAbin" and "AAbin" objects).
#' @param seqweights either NULL (all sequences are given weights
#'   of 1), a numeric vector the same length as \code{x} representing
#'   the sequence weights used to derive the model, or a character string giving
#'   the method to derive the weights from the sequences
#'   (see \code{\link{weight}}).
#' @param wfactor numeric. The factor to multiply the sequence weights by.
#'   Defaults to 1.
#' @param k integer representing the k-mer size to be used in tree-based
#'   sequence weighting (if applicable). Defaults to 5. Note that higher
#'   values of k may be slow to compute and use excessive memory due to
#'   the large numbers of calculations required.
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
#' @param endchar the character used to represent unknown residues in
#'   the alignment matrix (if applicable). Ignored for \code{"DNAbin"} or
#'   \code{"AAbin"} objects. Defaults to "?" otherwise.
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
#' @param logspace logical indicating whether the emission and transition
#'   probabilities in the returned model should be logged. Defaults to TRUE.
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
#' @param maxsize integer giving the upper bound on the number of modules
#'   in the PHMM. If NULL (default) no maximum size is enforced.
#' @param inserts character string giving the model construction method
#'   by which alignment columns
#'   are marked as either match or insert states. Accepted methods include
#'   \code{"threshold"} (only columns with fewer than a specified
#'   proportion of gaps form match states in the model), \code{"map"} (default;
#'   match and insert columns are found using the maximum \emph{a posteriori}
#'   method outlined in Durbin et al (1998) chapter 5.7), \code{"inherited"}
#'   (match and insert columns are inherited from the input alignment),
#'   and \code{"none"} (all columns are assigned
#'   match states in the model). Alternatively, insert columns can be
#'   specified manually by providing a logical vector the same length
#'   as the number of columns in the alignment, with \code{TRUE} for insert
#'   columns and \code{FALSE} for match states.
#' @param threshold the maximum proportion of gaps for an alignment column
#'   to be considered for a match state in the PHMM (defaults to 0.5).
#'   Only applicable when \code{inserts = "threshold"}.
#'   Note that the maximum \emph{a posteriori}
#'   method works poorly for alignments with few sequences,
#'   so the 'threshold' method is
#'   automatically used when the number of sequences is less than 5.
#' @param lambda penalty parameter used to favour models with fewer match
#'   states. Equivalent to the log of the prior probability of marking each
#'   column (Durbin et al 1998, chapter 5.7). Only applicable when
#'   \code{inserts = "map"}.
#' @param DI logical indicating whether delete-insert transitions should be
#'   allowed in the profile hidden Markov model (if applicable). Defaults
#'   to FALSE.
#' @param ID logical indicating whether insert-delete transitions should be
#'   allowed in the profile hidden Markov model (if applicable). Defaults to
#'   FALSE.
#' @param omit.endgaps logical. Should gap characters at each end of the
#'   sequences be ignored when deriving the transition probabilities
#'   of the model? Defaults to FALSE.
#'   Set to TRUE if \code{x} is not a strict global alignment
#'   (i.e. if the alignment contains partial sequences with missing
#'   sections represented with gap characters).
#' @param name an optional character string. The name of the
#'   new profile hidden Markov model.
#' @param description an optional character string. The description of the
#'   new profile hidden Markov model.
#' @param compo logical indicating whether the average emission
#'   probabilities of the model modules should be returned with the
#'   PHMM object.
#' @param consensus placeholder. Consensus sequences will be available in
#'   a future version.
#' @param alignment logical indicating whether the alignment used to
#'   derive the final model (if applicable) should be included as an element of
#'   the returned PHMM object. Defaults to FALSE.
#' @param progressive logical indicating whether the alignment used
#'   to derive the initial model parameters
#'   should be built progressively (assuming input is a list of
#'   unaligned sequences, ignored otherwise).
#'   Defaults to FALSE, in which case the
#'   longest sequence or sequences are used (faster,
#'   but possibly less accurate).
#' @param seeds optional integer vector indicating which sequences should
#'   be used as seeds for building the guide tree for the progressive
#'   alignment (assuming input is a list of unaligned sequences,
#'   and \code{progressive = TRUE}, ignored otherwise).
#'   Defaults to NULL, in which a set of log(n, 2)^2 non-identical
#'   sequences are chosen from the list of sequences by k-means clustering.
#' @param refine the method used to iteratively refine the model parameters
#'   following the initial progressive alignment and model derivation step.
#'   Current supported options are \code{"Viterbi"} (Viterbi training;
#'   the default option), \code{"BaumWelch"} (a modified version of the
#'   Expectation-Maximization algorithm), and "none" (skips the model
#'   refinement step).
#' @param deltaLL numeric, the maximum change in log likelihood between EM
#'   iterations before the cycling procedure is terminated (signifying model
#'   convergence). Defaults to 1E-07. Only applicable if
#'   \code{method = "BaumWelch"}.
#' @param maxiter the maximum number of EM iterations or Viterbi training
#'   iterations to carry out before the cycling process is terminated and
#'   the partially trained model is returned. Defaults to 100.
#' @param cpp logical, indicates whether the dynamic programming matrix
#'   should be filled using compiled C++ functions (default; many times faster).
#'   The FALSE option is primarily retained for bug fixing and experimentation.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if x is not a list.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#'   Note that in this case there
#'   may be a tradeoff in terms of speed depending on the number and size
#'   of sequences to be aligned, due to the extra time required to initialize
#'   the cluster.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @param ... aditional arguments to be passed to \code{"Viterbi"} (if
#'   \code{refine = "Viterbi"}) or \code{"forward"} (if
#'   \code{refine = "BaumWelch"}).
#' @return an object of class \code{"PHMM"}
#' @details
#'   This function performs a similar operation to the  \code{hmmbuild}
#'   function in the \href{http://www.hmmer.org}{HMMER} package, and the
#'   \code{modelfromalign} and \code{buildmodel} functions in the
#'   \href{https://compbio.soe.ucsc.edu/sam.html}{SAM} package.
#'   If the primary input argument is an alignment, the function creates a
#'   profile hidden Markov model (object class:\code{"PHMM"}) using the
#'   method described in Durbin et al (1998) chapter 5.3. Alternatively, if
#'   a list of non-aligned sequences is passed, the sequences are first aligned
#'   using the \code{\link{align}} function before being used to derive the
#'   model.
#'
#'   The function outputs an object of class \code{"PHMM"}, which is a list
#'   consisting of emission and transition probability matrices
#'   (elements named "E" and "A"), vectors of non-position-specific
#'   background emission and transition probabilities
#'   ("qe" and "qa", respectively) and other model metadata including
#'   "name", "description", "size" (the number of modules in the model), and
#'   "alphabet" (the set of symbols/residues emitted by the model).
#'
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#' @seealso \code{\link{deriveHMM}}, \code{\link{map}}
#' @examples
#' ## Small globin alignment data from Durbin et al (1998) Figure 5.3
#' data(globins)
#' ## derive a profile hidden Markov model from the alignment
#' globins.PHMM <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
#' plot(globins.PHMM, main = "Profile HMM for small globin alignment")
#' ##
#' ## derive a profle HMM from the woodmouse dataset in the
#' ## ape package and plot the first 5 modules
#' library(ape)
#' data(woodmouse)
#' woodmouse.PHMM <- derivePHMM(woodmouse)
#' plot(woodmouse.PHMM, from = 0, to = 5, main = "Partial woodmouse profile HMM")
#' @name derivePHMM
################################################################################
derivePHMM <- function(x, ...){
  UseMethod("derivePHMM")
}
################################################################################
#' @rdname derivePHMM
################################################################################
derivePHMM.DNAbin <- function(x, seqweights = "Henikoff", wfactor = 1, k = 5,
                              residues = NULL, gap = "-", endchar = "?",
                              pseudocounts = "background", logspace = TRUE,
                              qa = NULL, qe = NULL, maxsize = NULL,
                              inserts = "map", threshold = 0.5, lambda = 0,
                              DI = FALSE, ID = FALSE,
                              omit.endgaps = FALSE, name = NULL,
                              description = NULL, compo = FALSE,
                              consensus = FALSE, alignment = FALSE,
                              progressive = FALSE, seeds = NULL,
                              refine = "Viterbi",
                              maxiter = 100, deltaLL = 1E-07, cpp = TRUE,
                              cores = 1, quiet = FALSE, ...){
  if(is.list(x)){
    derivePHMM.list(x, progressive = progressive, seeds = seeds,
                    refine = refine, maxiter = maxiter,
                    deltaLL = deltaLL,
                    seqweights = seqweights, wfactor = wfactor,
                    k = k, residues = residues, gap = gap,
                    pseudocounts = pseudocounts,
                    logspace = logspace, qa = qa, qe = qe,
                    maxsize = maxsize,
                    inserts = inserts, lambda = lambda, DI = DI, ID = ID,
                    threshold = threshold,
                    omit.endgaps = omit.endgaps, name = name,
                    description = description, compo = compo,
                    consensus = consensus, alignment = alignment,
                    cpp = cpp, cores = cores, quiet = quiet, ... = ...)
  }else{
    derivePHMM.default(x, seqweights = seqweights, wfactor = wfactor,
                       k = k, residues = residues, gap = gap,
                       endchar = endchar, pseudocounts = pseudocounts,
                       logspace = logspace, qa = qa, qe = qe,
                       maxsize = maxsize, inserts = inserts,
                       threshold = threshold, lambda = lambda,
                       DI = DI, ID = ID, omit.endgaps = omit.endgaps,
                       name = name, description = description,
                       compo = compo, consensus = consensus,
                       alignment = alignment, cpp = cpp, quiet = quiet)
  }
}
################################################################################
#' @rdname derivePHMM
################################################################################
derivePHMM.AAbin <- function(x, seqweights = "Henikoff", wfactor = 1, k = 5,
                             residues = NULL, gap = "-", endchar = "?",
                             pseudocounts = "background", logspace = TRUE,
                             qa = NULL, qe = NULL, maxsize = NULL,
                             inserts = "map", threshold = 0.5, lambda = 0,
                             DI = FALSE, ID = FALSE, omit.endgaps = FALSE,
                             name = NULL, description = NULL, compo = FALSE,
                             consensus = FALSE,  alignment = FALSE,
                             progressive = FALSE, seeds = NULL,
                             refine = "Viterbi", maxiter = 100,
                             deltaLL = 1E-07, cpp = TRUE, cores = 1,
                             quiet = FALSE, ...){
  if(is.list(x)){
    derivePHMM.list(x, progressive = progressive, seeds = seeds,
                    refine = refine, maxiter = maxiter,
                    deltaLL = deltaLL,
                    seqweights = seqweights, wfactor = wfactor,
                    k = k, residues = residues, gap = gap,
                    pseudocounts = pseudocounts,
                    logspace = logspace, qa = qa, qe = qe, maxsize = maxsize,
                    inserts = inserts, lambda = lambda, DI = DI, ID = ID,
                    threshold = threshold, omit.endgaps = omit.endgaps,
                    name = name, description = description, compo = compo,
                    consensus = consensus, alignment = alignment,
                    cpp = cpp, cores = cores, quiet = quiet, ... = ...)
  }else{
    derivePHMM.default(x, seqweights = seqweights, wfactor = wfactor, k = k,
                       residues = residues, gap = gap, endchar = endchar,
                       pseudocounts = pseudocounts, logspace = logspace,
                       qa = qa, qe = qe, maxsize = maxsize, inserts = inserts,
                       threshold = threshold, lambda = lambda, DI = DI, ID = ID,
                       omit.endgaps = omit.endgaps, name = name,
                       description = description, compo = compo,
                       consensus = consensus, alignment = alignment, cpp = cpp,
                       quiet = quiet)
  }
}
################################################################################
#' @rdname derivePHMM
################################################################################
derivePHMM.list <- function(x, progressive = FALSE, seeds = NULL,
                            refine = "Viterbi", maxiter = 100, deltaLL = 1E-07,
                            seqweights = "Henikoff", wfactor = 1, k = 5,
                            residues = NULL, gap = "-",
                            pseudocounts = "background", logspace = TRUE,
                            qa = NULL, qe = NULL, maxsize = NULL,
                            inserts = "map", lambda = 0, DI = FALSE, ID = FALSE,
                            threshold = 0.5, omit.endgaps = FALSE,
                            name = NULL, description = NULL, compo = FALSE,
                            consensus = FALSE, alignment = FALSE, cpp = TRUE,
                            cores = 1, quiet = FALSE, ...){
  nseq <- length(x)
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  for(i in 1:nseq) x[[i]] <- x[[i]][x[[i]] != gap]
  ## set up multithread
  if(inherits(cores, "cluster")){
    para <- TRUE
    stopclustr <- FALSE
  }else if(cores == 1){
    para <- FALSE
    stopclustr <- FALSE
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(cores > 1){
      # if(cores > navailcores) stop("No. cores is more than number available")
      if(!quiet) cat("Multithreading over", cores, "cores\n")
      cores <- parallel::makeCluster(cores)
      para <- TRUE
      stopclustr <- TRUE
    }else{
      para <- FALSE
      stopclustr <- FALSE
    }
  }
  ## calculate sequence weights and initial alignment
  if(nseq > 2){
    if(progressive){
      catchnames <- names(x)
      names(x) <- paste0("S", 1:nseq)
      if(is.null(seeds)) seeds <- seq_along(x)
      stopifnot(
        mode(seeds) %in% c("numeric", "integer"),
        max(seeds) <= nseq,
        min(seeds) > 0
        ## just a check - can remove eventually
      )
      if(is.null(seqweights)){
        seqweights <- rep(1, nseq)
        names(seqweights) <- catchnames
      }else if(identical(seqweights, "Henikoff")){
        if(!quiet) cat("Calculating sequence weights using maximum entropy method\n")
        seqweights <- weight(x, method = "Henikoff", k = k,
                             residues = residues, gap = gap)
        names(seqweights) <- catchnames
      }else if(identical(seqweights, "Gerstein")){
        if(!quiet) cat("Calculating sequence weights using tree-based method\n")
        seqweights <- weight(x, method = "Gerstein", k = k,
                             residues = residues, gap = gap)
        names(seqweights) <- catchnames
      }else{
        stopifnot(mode(seqweights) %in% c("numeric", "integer"),
                  length(seqweights) == nseq)
      }
      guidetree <- kmer::cluster(x[seeds], k = k,
                                      residues = residues, gap = gap)
      attachseqs <- function(tree, sequences){
        if(!is.list(tree)){
          attr(tree, "seqs") <- sequences[[attr(tree, "label")]]
        }
        return(tree)
      }
      guidetree <- dendrapply(guidetree, attachseqs, sequences = x)
      progressive2 <- function(tree, maxsize, ...){
        if(is.list(tree)){
          if(!is.null(attr(tree[[1]], "seqs")) &
             !is.null(attr(tree[[2]], "seqs"))){
            attr(tree, "seqs") <- align.default(attr(tree[[1]], "seqs"),
                                                attr(tree[[2]], "seqs"),
                                                maxsize = maxsize, ... = ...)
            attr(tree[[1]], "seqs") <- attr(tree[[2]], "seqs") <- NULL
          }
        }
        return(tree)
      }
      progressive1 <- function(tree, maxsize, ...){
        tree <- progressive2(tree, maxsize = maxsize, ... = ...)
        if(is.list(tree)) tree[] <- lapply(tree, progressive1,
                                           maxsize = maxsize, ... = ...)
        return(tree)
      }
      if(!quiet) cat("Progressively aligning sequences\n")
      while(is.null(attr(guidetree, "seqs"))){
        guidetree <- progressive1(guidetree, maxsize = maxsize, ... = ...)
      }
      msa1 <- attr(guidetree, "seqs")
      rownames(msa1) <- catchnames[match(rownames(msa1), paste0("S", 1:nseq))]
      names(x) <- catchnames
    }else{
      if(is.null(seqweights)){
        seqweights <- rep(1, nseq)
        names(seqweights) <- names(x)
      }else if(identical(seqweights, "Henikoff") | identical(seqweights, "Gerstein")){
        seqweights <- weight(x, method = seqweights, k = k, residues = residues, gap = gap)
      }else{
        stopifnot(mode(seqweights) %in% c("numeric", "integer"),
                  length(seqweights) == nseq)
      }
      xlengths <- sapply(x, length)
      ## model length is based on max frequency
      lm <- as.numeric(names(sort(table(xlengths), decreasing = TRUE)[1]))
      if(!is.null(maxsize)){
        xlengths2 <- xlengths[xlengths <= maxsize]
        if(length(xlengths2) == 0) stop("maxsize parameter is too low\n")
        if(lm > maxsize) lm <- max(xlengths2)
      }
      longl <- xlengths == lm
      # longl <- xlengths == max(xlengths) #logical
      # changed Aug2017 to reflect frequency rather than length
      seeds <- which.min(seqweights[longl]) #index (or indices)
      if(length(seeds) > 1) seeds <- sample(seeds, size = 1) # length 1 index
      seed <- x[longl][[seeds]]
      msa1 <- matrix(seed, nrow = 1)
      colnames(msa1) <- paste(1:ncol(msa1))
      ## colnames are used when inserts = "inherited"
    }
  }else if(nseq == 2){
    if(!quiet) cat("Aligning seed sequences\n")
    msa1 <- align.default(x[[1]], x[[2]], residues = residues,
                          gap = gap, ... = ...)
    seqweights <- c(1, 1)
    rownames(msa1) <- names(seqweights) <- names(x)
    colnames(msa1) <- seq_len(ncol(msa1))
    seeds <- 1:2
  }else if(nseq == 1){
    msa1 <- matrix(x[[1]], nrow = 1)
    colnames(msa1) <- seq_len(ncol(msa1))
    ## colnames are used when inserts = "inherited"
    seqweights <- 1
    rownames(msa1) <- names(seqweights) <- names(x)
    seeds <- 1
  }else stop("Empty list")
  ## derive PHMM from alignnment
  if(!quiet) cat("Deriving profile HMM\n")
  model <- derivePHMM.default(msa1, seqweights = seqweights[seeds],
                              wfactor = wfactor, k = k, residues = residues,
                              gap = gap, pseudocounts = pseudocounts,
                              logspace = logspace, qa = qa, qe = qe,
                              DI = DI, ID = ID, omit.endgaps = omit.endgaps,
                              maxsize = maxsize, inserts = inserts,
                              lambda = lambda, threshold = threshold,
                              name = name, description = description,
                              compo = compo, consensus = consensus,
                              alignment = alignment, cpp = cpp,
                              quiet = quiet)
  if(nseq < 3){
    if(!quiet) cat("Done\n")
    if(para & stopclustr) parallel::stopCluster(cores)
    return(model)
  }
  ## train model
  if(is.null(refine)) refine <- "none"
  if(refine %in% c("Viterbi", "BaumWelch")){
    if(!quiet) cat("Refining model\n")
    model <- train.PHMM(model, x, seqweights = seqweights, method = refine,
                   maxiter = maxiter, deltaLL = deltaLL,
                   pseudocounts = pseudocounts, maxsize = maxsize,
                   inserts = inserts, lambda = lambda, threshold = threshold,
                   alignment = alignment, cores = cores, quiet = quiet,
                   cpp = cpp, ... = ...)
  }else {
    stopifnot(identical(refine, "none"))
    # Next line communicates with align.list and ensures all seqs get
    # aligned back to the model
    if(length(seeds) < nseq) attr(model, "alignment") <- NULL
  }
  if(para & stopclustr) parallel::stopCluster(cores)
  if(!quiet) cat("Done\n")
  return(model)
}
################################################################################
#'@rdname derivePHMM
################################################################################
derivePHMM.default <- function(x, seqweights = "Henikoff", wfactor = 1, k = 5,
                               residues = NULL, gap = "-", endchar = "?",
                               pseudocounts = "background", logspace = TRUE,
                               qa = NULL, qe = NULL, maxsize = NULL,
                               inserts = "map", lambda = 0, threshold = 0.5,
                               DI = FALSE, ID = FALSE, omit.endgaps = FALSE,
                               name = NULL, description = NULL, compo = FALSE,
                               consensus = FALSE, alignment = FALSE, cpp = TRUE,
                               quiet = FALSE,
                               ...){
  # if(!quiet) cat("Deriving profile HMM from alignment\n")
  #if(!(is.matrix(x))) stop("invalid object type, x must be a matrix")
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  catchnames <- rownames(x)
  rownames(x) <- paste0("S", 1:nrow(x)) ## just ensures names are unique
  DNA <- .isDNA(x) # raw DNA bytes
  AA <- .isAA(x) # raw AA bytes
  gap <- if(DNA) as.raw(4) else if(AA) as.raw(45) else gap
  endchar <- if(DNA) as.raw(2) else if(AA) as.raw(63) else endchar
  if(omit.endgaps) x <- .trim(x, gap = gap, endchar = endchar, DNA = DNA, AA = AA)
  residues <- .alphadetect(x, residues = residues, gap = gap, endchar = endchar)
  if(is.null(name)) name <- unname(sapply(match.call()[2], deparse))
  nres <- length(residues)
  n <- nrow(x)
  m <- ncol(x)
  states <- c("D", "M", "I")
  transitions <- c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")

  if(is.null(seqweights)){
    seqweights <- rep(1, n)
  }else if(identical(seqweights, "Gerstein") | identical(seqweights, "Henikoff")){
    seqlist <- unalign(x, gap = gap)
    seqlengths <- vapply(seqlist, length, 0L)
    if(n > 2 & min(seqlengths) > k + 1L){
      seqweights <- weight(seqlist, method = seqweights, k = k, residues = residues, gap = gap)
    }else{
      if(!quiet) cat("Applying uniform sequence weights\n")
      seqweights <- rep(1, n)
    }
  }else{
    stopifnot(
      length(seqweights) == n,
      !any(is.na(seqweights)),
      mode(seqweights) %in% c("numeric", "integer")
    )
  }
  if(wfactor != 1) seqweights <- seqweights * wfactor
  # background emission probabilities (qe)
  # if(!quiet) cat("Finding background emission probabilities\n")
  if(is.null(qe)){
    allecs <- if(AA){
      apply(x, 2, .tabulateAA, ambiguities = TRUE, seqweights = seqweights)
    }else if(DNA){
      apply(x, 2, .tabulateDNA, ambiguities = TRUE, seqweights = seqweights)
    }else{
      apply(x, 2, .tabulateCH, residues = residues, seqweights = seqweights)
    }
    allecs <- apply(allecs, 1, sum)
    qe <- (allecs + 1)/sum(allecs + 1)
  }else{
    if(!(is.vector(qe) & length(qe) == length(residues))) stop("qe invalid")
    if(is.null(names(qe))) stop("qe argument is missing names attribute")
    if(!(identical(names(qe), residues))) stop("qe names must match residues")
    if(all(qe <= 0)) qe <- exp(qe)
    if(!(round(sum(qe), 2) == 1)) stop("Background emissions (qe) must sum to 1")
    if(any(qe == 0)) warning("At least one background emission probability = 0")
    qe <- qe/sum(qe) #account for rounding errors
  }
  # designate insert-columns
  gaps <- x == gap
  #ends <- x == endchar
  #gapweights <- gaps * seqweights
  sws <- sum(seqweights)
  gpws <- numeric(m)
  for(i in seq_along(gpws)) gpws[i] <- sum(seqweights[gaps[, i]])
  # if(!quiet) cat("Marking insert states\n")
  if(identical(inserts, "none")){
    inserts <- rep(FALSE, m)
  }else if(identical(inserts, "inherited")){
    # inserts <- attr(x, "inserts")
    inserts <- colnames(x) == "I"
    if(length(inserts) == 0L) inserts <- rep(FALSE, ncol(x)) # for 2-row alignments
  }else if(identical(inserts, "threshold")){
    # inserts <- apply(gapweights, 2, sum) > threshold * n
    inserts <- gpws > threshold * sws
  }else if(identical(inserts, "map")){
    if(n < 5 | sum(gaps)/length(gaps) > 0.9){
      # Maximum a posteriori insert assignment unsuitable for
      # fewer than five sequences.
      # Also can use too much memory for very gappy (sparse) alignments
      # inserts <- apply(gapweights, 2, sum) > threshold * n
      inserts <- gpws > threshold * sws
    }else{
      inserts <- !map(x, seqweights = seqweights, residues = residues,
                      gap = gap, endchar = endchar, pseudocounts = pseudocounts,
                      qa = qa, qe = qe, cpp = cpp)
      # if(sum(!inserts) < 3) inserts <- apply(gapweights, 2, sum) > threshold * n
      if(sum(!inserts) < 3) inserts <- gpws > threshold * sws
    }
  }else if(!(mode(inserts) == "logical" & length(inserts) == ncol(x))){
    stop("invalid inserts argument")
  }
  l <- sum(!inserts) # PHMM length (excluding B & E positions)
  if(!is.null(maxsize)){
    maxsize <- as.integer(maxsize)
    if(maxsize < 3) stop("maxsize is too small")
    if(maxsize < l){
      # gapnos <- apply(gapweights[, !inserts, drop = FALSE], 2, sum)
      gapnos <- gpws[!inserts]
      inserts[!inserts][order(gapnos)[seq(maxsize + 1, l)]] <- TRUE
      l <- sum(!inserts)
    }
  }
  # rm(gapweights)
  ### emission counts (redo now that we know insert positions)
  # if(!quiet) cat("Finding emission probabilities\n")
  ecs <- if(AA){
    apply(x[, !inserts, drop = FALSE], 2, .tabulateAA,
          ambiguities = TRUE, seqweights = seqweights)
  }else if(DNA){
    apply(x[, !inserts, drop = FALSE], 2, .tabulateDNA,
          ambiguities = TRUE, seqweights = seqweights)
  }else{
    apply(x[, !inserts, drop = FALSE], 2, .tabulateCH,
          residues = residues, seqweights = seqweights)
  }
  if(length(ecs) > 0){
    dimnames(ecs) <- list(residue = residues, position = 1:l)
  }else ecs <- NULL
  ### transitions
  # if(!quiet) cat("Finding transition probabilities\n")
  xtr <- matrix(NA_integer_, nrow = n, ncol = m)
  insertsn <- matrix(rep(inserts, n), nrow = n, byrow = TRUE)
  xtr[gaps & !insertsn] <- 0L # Delete
  xtr[!gaps & !insertsn] <- 1L # Match
  xtr[!gaps & insertsn] <- 2L # Insert
  rm(gaps)
  rm(insertsn)
  xtr <- cbind(1L, xtr, 1L) # append begin and end match states
  tcs <- .atab(xtr, seqweights = seqweights)
  rm(xtr)
  #gc()
  alltcs <- apply(tcs, 1, sum) + 1 # forced addition of Laplacian pseudos
  ### background transition probs
  if(is.null(qa)){
    if(!DI) alltcs["DI"] <- 0
    if(!ID) alltcs["ID"] <- 0
    qa <- (alltcs)/sum(alltcs)
  }else{
    if(!is.vector(qa) | length(qa) != 9) stop("qa must be a length 9 vector")
    if(all(qa <= 0)) qa <- exp(qa)
    if(round(sum(qa), 2) != 1) stop("qa vector must sum to 1")
    if(!DI & (qa[3] != 0)){
      stop("DI is set to FALSE but delete -> insert transitions are
           assigned non-zero probabilities in qa vector.
           Change DI to TRUE or change qa['DI'] to zero")
    }
    if(!ID & (qa[7] != 0)){
      stop("ID is set to FALSE but insert -> delete transitions are
           assigned non-zero probabilities in qa vector.
           Change ID to TRUE or change qa['ID'] to zero")
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
  }else if(mode(pseudocounts) %in% c("numeric", "integer") & length(pseudocounts) == 1L){
    pseudocounts <- list(A = rep(pseudocounts, 9), E = rep(pseudocounts, nres))
    pseudocounts$A[3] <- pseudocounts$A[3] * DI
    pseudocounts$A[7] <- pseudocounts$A[7] * ID
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
  for(i in c(1, 4, 7)){
    A[, i:(i + 2)] <- A[, i:(i + 2)]/apply(A[, i:(i + 2), drop = FALSE], 1, sum)
  }
  A[1, 1:3] <- 0 # gets rid of NaNs caused by division by zero
  A <- t(A)
  inslens <- .insertlengths(!inserts)
  # which alignment columns correspond to which model positions
  whichcols <- which(!inserts)
  if(length(whichcols) > 0) names(whichcols) <- 1:l
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
  if(!is.null(catchnames)) names(seqweights) <- catchnames
  res <- structure(list(name = name, description = description,
                        size = l, alphabet = alphabet,
                        A = A, E = E, qa = qa, qe = qe, inserts = inserts,
                        insertlengths = inslens,
                        map = whichcols,
                        date = date(), nseq = n, weights = seqweights,
                        reference = reference, mask = mask),
                   class = "PHMM")
  # if(consensus){
  #   res$consensus <- generate(res, size = l * 100, random = FALSE, gap = ".")
  # }
  if(alignment){
    rownames(x) <- catchnames
    res$alignment <- x
  }
  if(compo) res$compo <- log(apply(exp(E), 1, mean))
  #gc()
  # if(!quiet) cat("Done\n")
  return(res)
}
################################################################################
#' Optimized profile HMM construction.
#'
#' Assigns match and insert states to alignment columns using the maximum
#'   \emph{a posteriori} algorithm outlined in Durbin et al (1998) chapter 5.7.
#'
#' @param x a matrix of aligned sequences. Accepted modes are "character"
#'   and "raw" (the latter being used for "DNAbin" and "AAbin" objects).
#' @param seqweights either NULL (default; all sequences are given
#'   weights of 1) or a numeric vector the same length as \code{x} representing
#'   the sequence weights used to derive the model.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Furthermore, the default setting
#'   \code{residues = NULL} will not detect rare residues that are not present
#'   in the sequences, and thus will not assign them emission probabilities
#'   in the model. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param endchar the character used to represent unknown residues in
#'   the alignment matrix (if applicable). Ignored for \code{"DNAbin"} or
#'   \code{"AAbin"} objects. Defaults to "?" otherwise.
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
#' @param lambda penalty parameter used to favour models with fewer match
#'   states. Equivalent to the log of the prior probability of marking each
#'   column (Durbin et al 1998, chapter 5.7).
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
#' @param cpp logical, indicates whether the dynamic programming matrix
#'   should be filled using compiled C++ functions (default; many times faster).
#'   The FALSE option is primarily retained for bug fixing and experimentation.
#' @return a logical vector with length = ncol(x) indicating the columns to be
#'   assigned as match states (\code{TRUE}) and those assigned as inserts
#'   (\code{FALSE}).
#' @details see Durbin et al (1998) chapter 5.7 for details of
#'   the maximum \emph{a posteriori} algorithm for optial match and insert
#'   state assignment.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @seealso \code{\link{derivePHMM}}
#' @examples
#' ## Maximum a posteriori assignment of match states to the small
#' ## alignment example in Figure 5.3, Durbin et al (1998)
#' data(globins)
#' map(globins)
################################################################################
map <- function(x, seqweights = NULL, residues = NULL,
                gap = "-", endchar = "?", pseudocounts = "background",
                lambda = 0, qa = NULL, qe = NULL, cpp = TRUE){
  if(!is.matrix(x)) stop("x must be a matrix")
  L <- ncol(x)
  n <- nrow(x)
  if(n < 4) return(structure(rep(TRUE, L), names = paste(1:L)))
  AA <- .isAA(x)
  DNA <- .isDNA(x)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  endchar <- if(DNA) as.raw(2) else if(AA) as.raw(63) else endchar
  residues <- .alphadetect(x, residues = residues, gap = gap, endchar = endchar)
  nres <- length(residues)
  transitions = c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II")
  S <- sigma <- c(0, rep(NA, L + 1))
  if(is.null(seqweights)) seqweights <- rep(1, n)
  if(length(seqweights) != n) stop("invalid seqweights argument")
  ecs <- if(AA){
    apply(x, 2, .tabulateAA, ambiguities = TRUE, seqweights = seqweights)
  }else if(DNA){
    apply(x, 2, .tabulateDNA, ambiguities = TRUE, seqweights = seqweights)
  }else{
    apply(x, 2, .tabulateCH, residues = residues, seqweights = seqweights)
  }
  # ecs <- t(t(ecs) * seqweights)
  allecs <- apply(ecs, 1, sum)
  if(is.null(qe)) {
    qe <- log((allecs + 1)/sum(allecs + 1))
  }else if(all(qe >= 0 & qe <= 1) & round(sum(qe), 2) == 1){
    qe <- log(qe)
  }else if(any(qe > 0) | round(sum(exp(qe)), 2) != 1) stop("invalid qe")
  gaps <- x == gap
  #ends <- x == endchar
  notgaps <- cbind(TRUE, !gaps, TRUE)
  #if(!DI) alltcs[c(3, 7)] <- 0 # TODO
  if(is.null(qa)){
    gapweights <- gaps * seqweights
    inserts <- apply(gapweights, 2, sum) > 0.5 * nrow(x)
    xtr <- matrix(nrow = nrow(x), ncol = ncol(x))
    insertsn <- matrix(rep(inserts, n), nrow = n, byrow = TRUE)
    xtr[gaps & !insertsn] <- 0L # Delete
    xtr[!gaps & !insertsn] <- 1L # Match
    xtr[!gaps & insertsn] <- 2L # Insert
    xtr <- cbind(1L, xtr, 1L) # append begin and end match states
    tcs <- .atab(xtr, seqweights) # modules = sum(!inserts) + 2)
    alltcs <- apply(tcs, 1, sum)
    qa <- log((alltcs + 1)/sum(alltcs + 1)) # force addition of Lapl pseudos
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
  }else if(mode(pseudocounts) %in% c("numeric", "integer") & length(pseudocounts) == 1L){
    pseudocounts <- list(A = rep(pseudocounts, 9), E = rep(pseudocounts, nres))
  }else if(identical(pseudocounts, "background")){
    pseudocounts <- list(A = exp(qa) * 9, E = exp(qe) * nres)
  }else if(identical(pseudocounts, "Laplace")){
    pseudocounts <- list(A = rep(1, 9), E = rep(1, nres))
  }else if(identical(pseudocounts, "none")){
    pseudocounts <- list(A = rep(0, 9), E = rep(0, nres))
  }else stop("invalid pseudocounts argument")
  if(cpp){
    res <- .map(ecs, notgaps, pseudocounts, seqweights, qe, lambda)
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
          # avoid counting NAs in notgaps matrix (endchars)
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
      sigma[j] <- .whichmax(tmp)
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
  #gc()
  return(res)
}
################################################################################
