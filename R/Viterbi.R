#' The Viterbi algorithm.
#'
#' The \code{Viterbi} function finds the optimal path of a sequence through a HMM
#'   or PHMM and returns its full (log) probability or log-odds score.
#'
#' @param x an object of class \code{HMM} or \code{PHMM}. Optionally, both x and
#'   y can be sequences (character vectors or DNAbin/AAbin objects), in which
#'   case the operation becomes either the Needleman-Wunch (global algnment) or
#'   Smith-Waterman (local alignment) algorithm.
#' @param y a vector of mode "character" or "raw" (a "DNAbin" or "AAbin"
#'   object) representing a single sequence hypothetically emitted by
#'   the model in \code{x}. Optionally, both x and
#'   y can be profile hidden Markov models (object class "PHMM"), in which
#'   case the sum of log-odds algorithm of Soding (2005) is used.
#' @param qe an optional named vector of background residue frequencies (only
#'   applicable if x is a PHMM). If \code{qe = NULL} the function looks for
#'   a qe vector as an attribute of the PHMM. If these are not available
#'   equal background residue frequencies are assumed.
#' @param logspace logical indicating whether the emission and transition
#'   probabilities of x are logged. If \code{logspace = "autodetect"}
#'   (default setting), the function will automatically detect
#'   if the probabilities are logged, returning an error if
#'   inconsistencies are found. Note that choosing the latter option
#'   increases the computational overhead; therefore specifying
#'   \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param type character string indicating whether insert and delete states
#'   at the beginning and end of the path should count towards the final score
#'   ('global'; default), or not ('semiglobal'), or whether the highest scoring
#'   sub-path should be returned ('local').
#' @param odds logical, indicates whether the returned scores
#'   should be odds ratios (TRUE) or full logged probabilities (FALSE).
#' @param offset column score offset to specify level of greediness. Defaults to
#'   -0.1 bits for PHMM x PHMM alignments (as recommended by Soding (2005)), and 0
#'   otherwise.
#' @param d gap opening penalty (in bits) for sequence vs. sequence alignment.
#'   Defaults to 8.
#' @param e gap extension penalty (in bits) for sequence vs. sequence alignment.
#'   Defaults to 2.
#' @param S an optional scoring matrix with rows and columns named according
#'   to the residue alphabet. Only applicable when both x and y are sequences
#'   (Needleman-Wunch or Smith-Waterman alignments).
#'   Note that for Smith-Waterman local alignments, scores for
#'   mismatches should generally take negative values to avoid spurious
#'   alignments. If NULL default settings are used. Default scoring matrices are
#'   'NUC.4.4' for For DNAbin objects, and 'MATCH' (matches are scored 1 and
#'   mismatches are scored -1) for AAbin objects and character sequences.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors.
#' @param DI logical indicating whether delete-insert transitions should be
#'   allowed in the profile hidden Markov model (if applicable). Defaults
#'   to FALSE.
#' @param ID logical indicating whether insert-delete transitions should be
#'   allowed in the profile hidden Markov model (if applicable). Defaults to
#'   FALSE.
#' @param windowspace a two-element integer vector providing the search space for
#'   dynamic programming (see Wilbur & Lipman 1983 for details). The first element
#'   should be negative, and represent the lowermost diagonal of the
#'   dynammic programming array, and the second element should be positive,
#'   representing the leftmost diagonal. Alternatively, if the the character
#'   string "all" is passed (the default setting) the entire dynamic programming
#'   array will be computed.
#' @param cpp logical, indicates whether the dynamic programming matrix
#'   should be filled using compiled C++ functions (default; many times faster).
#'   The FALSE option is primarily retained for bug fixing and experimentation.
#' @param ... additional arguments to be passed between methods.
#' @return an object of class \code{"DPA"}, which is a list including the
#'   score, the dynammic programming array, and the optimal path (an integer
#'   vector, see details section).
#' @details
#'   This function is a wrapper for a compiled C++ function that recursively
#'   fills a dynamic programming matrix with probabilities, and
#'   calculates the (logged) probability and optimal path of a sequence
#'   through a HMM or PHMM.
#'
#'   If x is a PHMM and y is a sequence, the path is represented as
#'   an integer vector containing zeros, ones and twos, where a zero represents
#'   a downward transition, a one represents a diagonal transition downwards and
#'   left, and a two represents a left transition in the dynamic programming
#'   matrix (see Durbin et al (1998) chapter 2.3). This translates to
#'   0 = delete state, 1 = match state and 2 = insert state.
#'
#'   If x and y are both sequences, the function implements the
#'   Needleman-Wunch or Smith Waterman algorithm depending on
#'   the type of alignment specified. In this case, a zero
#'   in the path refers to x aligning to a gap in y, a one refers
#'   to a match, and a two refers to y aligning to a gap in x.
#'
#'   If x is a standard hidden Markov model (HMM) and y is a sequence,
#'   each integer in the path represents a state in the model.
#'   Note that the path elements can take values between 0 and
#'   one less than number of states, as in the C/C++ indexing
#'   style rather than R's.
#'
#'   For a thorough explanation of the backward, forward and Viterbi
#'   algorithms, see Durbin et al (1998) chapters 3.2 (HMMs) and 5.4 (PHMMs).
#'
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Soding J (2005) Protein homology detection by HMM-HMM comparison.
#'   \emph{Bioinformatics}, \strong{21}, 951-960.
#'
#'   Wilbur WJ, Lipman DJ (1983) Rapid similarity searches of nucleic acid and
#'   protein data banks. \emph{Proc Natl Acad Sci USA}, \strong{80}, 726-730.
#'
#' @seealso \code{\link{backward}}, \code{\link{forward}}, \code{\link{align}}
#' @examples
#'   ## Viterbi algorithm for standard HMMs:
#'   ## The dishonest casino example from Durbin et al (1998) chapter 3.2
#'   states <- c("Begin", "Fair", "Loaded")
#'   residues <- paste(1:6)
#'   ### Define the transition probability matrix
#'   A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9), nrow = 3)
#'   dimnames(A) <- list(from = states, to = states)
#'   ### Define the emission probability matrix
#'   E <- matrix(c(rep(1/6, 6), rep(1/10, 5), 1/2), nrow = 2, byrow = TRUE)
#'   dimnames(E) <- list(states = states[-1], residues = residues)
#'   ### Build and plot the HMM object
#'   x <- structure(list(A = A, E = E), class = "HMM")
#'   plot(x, main = "Dishonest casino HMM")
#'   ### Find optimal path of sequence
#'   data(casino)
#'   casino.DPA <- Viterbi(x, casino)
#'   casino.DPA$score # full (log) prob of sequence given model = -538.8109
#'   ### Show optinal path path as indices
#'   casino.DPA$path
#'   ### Show optimal path as character strings
#'   rownames(x$E)[casino.DPA$path + 1]
#'   ##
#'   ## Needleman-Wunch pairwise sequence alignment:
#'   ## Pairwise protein alignment example from Durbin et al (1998) chapter 2.3
#'   x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#'   y <- c("P", "A", "W", "H", "E", "A", "E")
#'   Viterbi(x, y,  d = 8, e = 2, type = "global")
#'   ###
#'   ## Viterbi algorithm for profile HMMs:
#'   ## Small globin alignment data from Durbin et al (1998) Figure 5.3
#'   data(globins)
#'   ### Derive a profile hidden Markov model from the alignment
#'   globins.PHMM <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
#'   plot(globins.PHMM, main = "Profile hidden Markov model for globins")
#'   ### Simulate a random sequence from the model
#'   suppressWarnings(RNGversion("3.5.0"))
#'   set.seed(999)
#'   simulation <- generate(globins.PHMM, size = 20)
#'   simulation ## "F" "S" "A" "N" "N" "D" "W" "E"
#'   ### Calculate the odds of the optimal path of the sequence given the model
#'   x <- Viterbi(globins.PHMM, simulation, odds = FALSE)
#'   x # -23.07173
#'   ### Show dynammic programming array
#'   x$array
#'   ### Show the optimal path as an integer vector
#'   x$path
#'   ### Show the optimal path as either delete states, matches or insert states
#'   c("D", "M", "I")[x$path + 1]
#'   ### Correctly predicted the actual path:
#'   names(simulation)
#' @name Viterbi
################################################################################
Viterbi <- function(x, y, ...){
  UseMethod("Viterbi")
}
################################################################################
#' @rdname Viterbi
################################################################################
Viterbi.PHMM <- function(x, y, qe = NULL, logspace = "autodetect",
                         type = "global", odds = TRUE, offset = 0,
                         windowspace = "all", DI = FALSE, ID = FALSE,
                         cpp = TRUE, ...){
  if(type != "global" & !odds) {
    stop("Non-odds option only available for global alignment")
  }
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  if(!(type %in% c("global", "semiglobal", "local"))) stop("Invalid type")
  pp <- inherits(y, "PHMM")
  pd <- .isDNA(y)
  pa <- .isAA(y)
  pc <- !pp & !pd & !pa
  if(!pp & is.list(y)) if(length(y) == 1) y <- y[[1]] else stop("y is invalid")
  if(pd){
    rownames(x$E) <- toupper(rownames(x$E))
    if("U" %in% rownames(x$E)) rownames(x$E)[rownames(x$E) == "U"] <- "T"
    NUCorder <- sapply(rownames(x$E), match, c("A", "T", "G", "C"))
    x$E <- x$E[NUCorder, ]
    if(!(identical(rownames(x$E), c("A", "T", "G", "C")))){
      stop("invalid model for DNA, residue alphabet does not correspond to
            nucleotide alphabet")
    }
    class(y) <- "DNAbin"
    y.DNAbin <- y
    y <- .encodeDNA(y, arity = 15, na.rm = TRUE)
  }else if(pa){
    rownames(x$E) <- toupper(rownames(x$E))
    PFAMorder <- sapply(rownames(x$E), match, LETTERS[-c(2, 10, 15, 21, 24, 26)])
    x$E <- x$E[PFAMorder, ]
    if(!(identical(rownames(x$E), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
      stop("Invalid model for AA, residue alphabet does not correspond to
            20-letter amino acid alphabet")
    }
    class(y) <- "AAbin"
    y.AAbin <- y
    y <- .encodeAA(y, arity = 27, na.rm = TRUE)
  }else if(pc){
    if(mode(y) == "character"){
      y <- match(y, rownames(x$E)) - 1
      if(any(is.na(y))){
        warning("Residues in sequence(s) are missing from the model")
        y <- y[!is.na(y)]
      }
    }
  }
  n <- ncol(x$E) + 1
  m <- if(pp) ncol(y$E) + 1 else length(y) + 1
  states <- if(pp) c("MI", "DG", "MM", "GD", "IM") else c("D", "M", "I")
  # background emission probabilities
  if(!(is.null(qe))){
    if(all(qe >= 0) & all(qe <= 1)) qe <- log(qe)
  }else if(pp){
    if(!is.null(x$qe) | !is.null(y$qe)){
      if(!is.null(x$qe) & !is.null(y$qe)){
        qe <- if(logspace) (exp(x$qe) + exp(y$qe))/2 else (x$qe + y$qe)/2
        qe <- log(qe)
      }else if(!is.null(x$qe)){
        qe <- if(logspace) x$qe else log(x$qe)
      }else{
        qe <- if(logspace) y$qe else log(y$qe)
      }
    }else{
      qe <- log(rep(1/nrow(x$E), nrow(x$E)))
      names(qe) <- rownames(x$E)
    }
  }else{
    if(!is.null(x$qe)){
      qe <- if(logspace) x$qe else log(x$qe)
    }else{
      qe <- log(rep(1/nrow(x$E), nrow(x$E)))
      names(qe) <- rownames(x$E)
    }
  }
  type = switch(type, "global" = 0L, "semiglobal" = 1L, "local" = 2L,
                stop("invalid type"))
  V <- array(-Inf, dim = c(n, m, length(states)))
  dimnames(V) <- list(x = 0:(n - 1), y = 0:(m - 1), state = states)
  # pointer array
  P <- V + NA
  # fill scoring and pointer arrays. Key: D = 1, M = 2, I = 3
  if(pp){
    if(!odds) stop("Full probabilities for PHMM-PHMM alignment are not supported")
    if(.logdetect(y) != logspace) {
      stop("Both models must be in the same logspace format")
    }
    if(identical(windowspace, "WilburLipman")){
      warning("No Wilbur Lipman method available for PHMM-PHMM comparison yet")
      windowspace <- c(-x$size, y$size)
      # TODO
      # xseq <- generate.PHMM(x, size = 10 * ncol(x$A), random = FALSE)
      # xseq <- setNames(seq_along(rownames(x$E)) - 1,  rownames(x$E))[xseq]
      # yseq <- generate.PHMM(y, size = 10 * ncol(y$A), random = FALSE)
      # yseq <- setNames(seq_along(rownames(x$E)) - 1,  rownames(x$E))[yseq]
      # #x$A to ensure order is same
      # windowspace <- .streak(xseq, yseq, arity = nrow(x$E), k = 2)
    }else if(identical(windowspace, "all")){
      windowspace <- c(-x$size, y$size)
    }else if(length(windowspace) != 2) stop("invalid windowspace argument")
    Ax <- if(logspace) x$A else log(x$A)
    Ay <- if(logspace) y$A else log(y$A)
    Ex <- if(logspace) x$E else log(x$E)
    Ey <- if(logspace) y$E else log(y$E)
    if(cpp){
      res <- .ViterbiPP(Ax, Ay, Ex, Ey, qe, type, windowspace, offset)
      V[, , 1] <- res$MImatrix
      V[, , 2] <- res$DGmatrix
      V[, , 3] <- res$MMmatrix
      V[, , 4] <- res$GDmatrix
      V[, , 5] <- res$IMmatrix
      P[, , 1] <- res$MIpointer
      P[, , 2] <- res$DGpointer
      P[, , 3] <- res$MMpointer
      P[, , 4] <- res$GDpointer
      P[, , 5] <- res$IMpointer
      res$array <- V
      res$pointer <- P
      res$MImatrix <- NULL
      res$DGmatrix <- NULL
      res$MMmatrix <- NULL
      res$GDmatrix <- NULL
      res$IMmatrix <- NULL
      res$MIpointer <- NULL
      res$DGpointer <- NULL
      res$MMpointer <- NULL
      res$GDpointer <- NULL
      res$IMpointer <- NULL
    }else{
      fun <- Vectorize(function(a, b) logsum((Ex[, a] + Ey[, b]) - qe))
      Saa <- outer(1:(n - 1), 1:(m - 1), fun)
      if(type == 0){
        V[-1, 1, "MI"] <- cumsum(c(Ax["MM", 1] + Ay["MI", 1], Ax["MM", 2:(n - 1)] + Ay["II", 1]))
        P[-1, 1, "MI"] <- c(3, rep(1, n - 2))
        V[-1, 1, "DG"] <- cumsum(c(Ax["MD", 1], Ax["DD", 2:(n - 1)]))
        P[-1, 1, "DG"] <- c(3, rep(2, n - 2))
        V[1, 1, "MM"] <- 0
        V[1, -1, "GD"] <- cumsum(c(Ay["MD", 1], Ay["DD", 2:(m - 1)]))
        P[1, -1, "GD"] <- c(3, rep(4, m - 2))
        V[1, -1, "IM"] <- cumsum(c(Ax["MI", 1] + Ay["MM", 1], Ax["II", 1] + Ay["MM", 2:(m - 1)]))
        P[1, -1, "IM"] <- c(3, rep(5, m - 2))
      }else{
        V[1, , "MM"] <- V[, 1, "MM"] <- 0
      }
      for(i in 2:n){
        for(j in 2:m){
          if(j - i >= windowspace[1] & j - i <= windowspace[2]){
            sij <- Saa[i - 1, j - 1] + offset
            MIcdt <- c(V[i - 1, j, "MM"] + Ax["MM", i - 1] + Ay["MI", j],
                       V[i - 1, j, "MI"] + Ax["MM", i - 1] + Ay["II", j])
            DGcdt <- c(V[i - 1, j, "MM"] + Ax["MD", i - 1],
                       V[i - 1, j, "DG"] + Ax["DD", i - 1])
            MMcdt <- c(V[i - 1, j - 1, "MI"] + Ax["MM", i - 1] + Ay["IM", j - 1] + sij,
                       V[i - 1, j - 1, "DG"] + Ax["DM", i - 1] + Ay["MM", j - 1] + sij,
                       V[i - 1, j - 1, "MM"] + Ax["MM", i - 1] + Ay["MM", j - 1] + sij,
                       V[i - 1, j - 1, "GD"] + Ax["MM", i - 1] + Ay["DM", j - 1] + sij,
                       V[i - 1, j - 1, "IM"] + Ax["IM", i - 1] + Ay["MM", j - 1] + sij,
                       if(type == 2) 0 else NULL)
            GDcdt <- c(V[i, j - 1, "MM"] + Ay["MD", j - 1],
                       V[i, j - 1, "GD"] + Ay["DD", j - 1])
            IMcdt <- c(V[i, j - 1, "MM"] + Ax["MI", i] + Ay["MM", j - 1],
                       V[i, j - 1, "IM"] + Ax["II", i] + Ay["MM", j - 1])
            MImax <- .whichmax(MIcdt)
            DGmax <- .whichmax(DGcdt)
            MMmax <- .whichmax(MMcdt)
            GDmax <- .whichmax(GDcdt)
            IMmax <- .whichmax(IMcdt)
            V[i, j, "MI"] <- MIcdt[MImax]
            V[i, j, "DG"] <- DGcdt[DGmax]
            V[i, j, "MM"] <- MMcdt[MMmax]
            V[i, j, "GD"] <- GDcdt[GDmax]
            V[i, j, "IM"] <- IMcdt[IMmax]
            P[i, j, "MI"] <- c(3, 1)[MImax]
            P[i, j, "DG"] <- c(3, 2)[DGmax]
            P[i, j, "MM"] <- c(1:6)[MMmax]
            P[i, j, "GD"] <- c(3, 4)[GDmax]
            P[i, j, "IM"] <- c(3, 5)[IMmax]
          }
        }
      }
      path <- c()
      progression <- matrix(nrow = 2, ncol = 0)
      if(type == 0){
        LLcdt <- c(V[n, m, "MI"] + Ax["MM", n] + Ay["IM", m],
                   V[n, m, "DG"] + Ax["DM", n] + Ay["MM", m],
                   V[n, m, "MM"] + Ax["MM", n] + Ay["MM", m],
                   V[n, m, "GD"] + Ax["MM", n] + Ay["DM", m],
                   V[n, m, "IM"] + Ax["IM", n] + Ay["MM", m])
        LLptr <- .whichmax(LLcdt)
        score <- LLcdt[LLptr]
        z <- c(n, m, LLptr)
        while(z[1] > 1 | z[2] > 1){
          path <- c(z[3], path)
          progression <- cbind(z[1:2], progression)
          z[3] <- P[z[1], z[2], z[3]]
          z <- z - switch(path[1], c(1,0,0), c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,0))
        }
        startposition = c(1, 1)
      }else if (type == 1){
        tmp <- V[, , "MM"]
        tmp[1:(n - 1), 1:(m - 1)] <- -Inf
        ind <- which(tmp == max(tmp), arr.ind = TRUE)
        if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
        z <- c(ind, 3)
        score <- V[z[1], z[2], z[3]]
        if(z[1] < n){
          path <- rep(2, n - z[1])
          progression <- rbind((z[1] + 1):n, rep(z[2], n - z[1]))
        }else if(z[2] < m){
          path <- rep(4, m - z[2])
          progression <- rbind((z[2] + 1):m, rep(z[1], m - z[2]))
        }
        while(z[1] > 1 & z[2] > 1){
          path <- c(z[3], path)
          progression <- cbind(z[1:2], progression)
          z[3] <- P[z[1], z[2], z[3]]
          z <- z - switch(path[1], c(1,0,0), c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,0))
        }
        if(z[1] > 1){
          path <- c(rep(2, z[1] - 1), path)
          progression <- cbind(rbind(1:(z[1] - 1), rep(1, z[1] - 1)), progression)
        }else if(z[2] > 1){
          path <- c(rep(4, z[2] - 1), path)
          progression <- cbind(rbind(rep(1, z[2] - 1), 1:(z[2] - 1)), progression)
        }
        startposition = c(1, 1)
      }else if(type == 2){
        ind <- which(V[, , "MM"] == max(V[, , "MM"]), arr.ind = TRUE)
        if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
        z <- c(ind, 3)
        score <- V[z[1], z[2], z[3]]
        P[1, 1, "MM"] <- 6
        while(P[z[1], z[2], z[3]] != 6){
          path <- c(z[3], path)
          progression <- cbind(z[1:2], progression)
          z[3] <- P[z[1], z[2], z[3]]
          z <- z - switch(path[1], c(1,0,0), c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,0))
        }
        startposition = z[1:2]
      }else(stop("Invalid alignment type"))
      #key <- "MI = 1, DG = 2, MM = 3, GD = 4, IM = 5"
      rownames(progression) <- c(deparse(substitute(x)), deparse(substitute(y)))
      #progression <- progression - 1 #to account for 0 row
      path <- path - 1
      res <- structure(list(score = score,
                            path = path,
                            #progression = progression,
                            start = startposition,
                            array = V,
                            pointer = P,
                            Saa = Saa),
                       class = "DPA")
    }
  }else{
    if(identical(windowspace, "WilburLipman")){
      xseq <- generate.PHMM(x, size = 10 * ncol(x$A), random = FALSE, AA = pa, DNA = pd)
      if(pd){
        xqt  <- match(xseq, as.raw(c(136, 24, 72, 40))) - 1
        yqt <- .encodeDNA(y.DNAbin, arity = 4, na.rm = TRUE)
        windowspace <- .streak(xqt, yqt, arity = 4, k = 5)
      }else if(pa){
        y.comp <- .encodeAA(y.AAbin, arity = 6, na.rm = TRUE)
        xseq.comp <- .encodeAA(xseq, arity = 6, na.rm = TRUE)
        windowspace <- .streak(xseq.comp, y.comp, arity = 6, k = 5)
      }else{
        xseq <- match(xseq, rownames(x$E)) - 1
        windowspace <- .streak(xseq, y, arity = nrow(x$E), k = 3)
      }
    }else if(identical(windowspace, "all")){
      windowspace <- c(-x$size, length(y))
    }else if(length(windowspace) != 2) stop("invalid windowspace argument")
    qey <- if(odds){
      rep(0, m - 1)
    }else if(pd){
      sapply(y, .probDNA, qe)
    }else if(pa){
      sapply(y, .probAA, qe)
    }else qe[y + 1]
    A <- if(logspace) x$A else log(x$A)
    E <- if(logspace) x$E else log(x$E)
    if(n == 1){
      V[1, 1, "M"] <- 0
      if(m > 1) V[1, 2, "I"] <- A["MI", 1]
      if(m > 2) for(j in 3:m) V[1, j, "I"] <- V[1, j - 1, "I"] + A["II", 1] + qey[j - 1]
      P[, -1, "I"] <- c(1, rep(2, m - 2))
      score <- if(m > 1) V[1, m, "I"] + A["IM", 1] else 0
      res <- structure(list(score = score,
                            path = rep(2, m - 1),
                            start = c(1, 1),
                            array = V,
                            pointer = P),
                       class = "DPA")
      return(res)
    }
    if(m == 1){
      V[1, 1, "M"] <- 0
      if(n > 1) V[2, 1, "D"] <- A["MD", 1]
      if(n > 2) for(i in 3:n) V[i, 1, "D"] <- V[i - 1, 1, "D"] + A["DD", i - 1]
      P[-1, , "D"] <- c(1, rep(0, n - 2))
      score <- if(n > 1) V[n, 1, "D"] + A["DM", n] else 0
      res <- structure(list(score = score,
                            path = rep(0, n - 1),
                            start = c(1, 1),
                            array = V,
                            pointer = P),
                       class = "DPA")
      return(res)
    }
    if(odds) E <- E - qe
    if(cpp){
      res <- .ViterbiP(y, A, E, qe, qey, type, windowspace, offset, DI, ID, DNA = pd, AA = pa)
      V[, , 1] <- res$Dmatrix
      V[, , 2] <- res$Mmatrix
      V[, , 3] <- res$Imatrix
      P[, , 1] <- res$Dpointer
      P[, , 2] <- res$Mpointer
      P[, , 3] <- res$Ipointer
      res$array <- V
      res$pointer <- P
      res$Dmatrix <- NULL
      res$Mmatrix <- NULL
      res$Imatrix <- NULL
      res$Dpointer <- NULL
      res$Mpointer <- NULL
      res$Ipointer <- NULL
    }else{
      P[, 1, 1] <- c(NA, 2, rep(1, n - 2))
      P[1, , 3] <- c(NA, 2, rep(3, m - 2))
      V[1, 1, 2] <- 0
      if(type == 0){
        V[-1, 1, 1] <- cumsum(c(0, A["DD", 2:(n - 1)])) + A["MD", 1]
        V[1, 2, "I"] <- A["MI", 1] + qey[1]
        for(j in 3:m) V[1, j, "I"] <- V[1, j - 1, "I"] + A["II", 1] + qey[j - 1]
      }else{
        V[-1, 1, 1] <- 0
        V[1, -1, 3] <- cumsum(qey)
      }
      for(i in 2:n){
        for(j in 2:m){
          if(j - i >= windowspace[1] & j - i <= windowspace[2]){
            if(pd){
              sij <- .probDNA(y[j - 1], E[, i - 1]) + offset
            }else if(pa){
              sij <- .probAA(y[j - 1], E[, i - 1]) + offset
            }else{
              sij <- E[y[j - 1] + 1, i - 1] + offset
            }
            Dcdt <- c(V[i - 1, j, "D"] + A["DD", i - 1],
                      V[i - 1, j, "M"] + A["MD", i - 1],
                      if(ID) V[i - 1, j, "I"] + A["ID", i - 1] else -Inf)
            Mcdt <- c(V[i - 1, j - 1, "D"] + A["DM", i - 1] + sij,
                      V[i - 1, j - 1, "M"] + A["MM", i - 1] + sij,
                      V[i - 1, j - 1, "I"] + A["IM", i - 1] + sij,
                      if(type == 2) 0 else NULL)
            Icdt <- c(if(DI) V[i, j - 1, "D"] + A["DI", i] else -Inf,
                      V[i, j - 1, "M"] + A["MI", i],
                      V[i, j - 1, "I"] + A["II", i])
            Dmax <- .whichmax(Dcdt)
            Mmax <- .whichmax(Mcdt)
            Imax <- .whichmax(Icdt)
            V[i, j, "D"] <- Dcdt[Dmax]
            V[i, j, "M"] <- Mcdt[Mmax]
            V[i, j, "I"] <- Icdt[Imax] + qey[j - 1]
            P[i, j, "D"] <- Dmax
            P[i, j, "M"] <- Mmax
            P[i, j, "I"] <- Imax
          }
        }
      }
      path <- c()
      progression <- c()
      if(type == 0){
        LLcdt <- c(V[n, m, "D"] + A["DM", n],
                   V[n, m, "M"] + A["MM", n],
                   V[n, m, "I"] + A["IM", n])
        LLptr <- .whichmax(LLcdt)
        score <- LLcdt[LLptr]
        z <- c(n, m, LLptr)
        while(z[1] > 1 | z[2] > 1){
          path <- c(z[3], path)
          progression <- c(z[1], progression)
          z[3] <- P[z[1], z[2], z[3]]
          z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
        }
        startposition = c(1, 1)
      }else if (type == 1){
        tmp <- V[, , "M"]
        tmp[1:(n - 1), 1:(m - 1)] <- -Inf
        ind <- which(tmp == max(tmp), arr.ind = TRUE)
        if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1), ]
        z <- c(ind, 2)
        score <- V[z[1], z[2], z[3]]
        if(z[1] < n){
          path <- rep(1, n - z[1])
          progression <- (z[1] + 1):n
        }else if(z[2] < m){
          path <- rep(3, m - z[2])
          progression <- rep(n, m - z[2])
        }
        while(z[1] > 1 & z[2] > 1){
          path <- c(z[3], path)
          progression <- c(z[1], progression)
          z[3] <- P[z[1], z[2], z[3]]
          z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
        }
        if(z[1] > 1){
          path <- c(rep(1, z[1] - 1), path)
          progression <- c(1:(z[1] - 1), progression)
        }else if(z[2] > 1){
          path <- c(rep(3, z[2] - 1), path)
          progression <- c(rep(1, z[2] - 1), progression)
        }
        startposition = c(1, 1)
      }else if(type == 2){
        ind <- which(V[, , "M"] == max(V[, , "M"]), arr.ind = TRUE)
        if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
        if(V[ind[1], ind[2], "M"] == 0) stop("Alignment failed, try increasing
                                           offset or changing type")
        z <- c(ind, 2)
        score <- V[z[1], z[2], z[3]]
        P[1, 1, "M"] <- 4
        while(P[z[1], z[2], z[3]] != 4){
          path <- c(z[3], path)
          progression <- c(z[1], progression)
          z[3] <- P[z[1], z[2], z[3]]
          z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
        }
        startposition = z[1:2]
      }
      #key: "1 = delete, 2 = match, 3 = insert"
      #progression <- progression - 1 #to account for 0 row
      path <- path - 1
      res <- structure(list(score = score,
                            path = path,
                            #progression = progression,
                            start = startposition,
                            array = V,
                            pointer = P),
                       class = "DPA")
    }
  }
  return(res)
}
################################################################################
#' @rdname Viterbi
################################################################################
Viterbi.HMM <- function (x, y, logspace = "autodetect", cpp = TRUE, ...){
  if(identical(logspace, 'autodetect')) logspace <- .logdetect(x)
  DNA <- .isDNA(y)
  AA <- .isAA(y)
  if(DNA){
    colnames(x$E) <- toupper(colnames(x$E))
    NUCorder <- sapply(colnames(x$E), match, c("A", "T", "G", "C"))
    x$E <- x$E[, NUCorder]
    if(!(identical(colnames(x$E), c("A", "T", "G", "C")))){
      stop("invalid model for DNA, residue alphabet does not correspond to
           nucleotide alphabet")
    }
    if(is.list(y)){
      if(length(y) == 1){
        y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "DNAbin"
      }else stop("Invalid input object y: multi-sequence list")
    }
    y <- .encodeDNA(y, arity = 15, na.rm = TRUE)
  }else if(AA){
    colnames(x$E) <- toupper(colnames(x$E))
    PFAMorder <- sapply(colnames(x$E), match, LETTERS[-c(2, 10, 15, 21, 24, 26)])
    x$E <- x$E[, PFAMorder]
    if(!(identical(colnames(x$E), LETTERS[-c(2, 10, 15, 21, 24, 26)]))){
      stop("invalid model for AA, residue alphabet does not correspond to
           20-letter amino acid alphabet")
    }
    if(is.list(y)){
      if(length(y) == 1){
        y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "AAbin"
      }else stop("Invalid input object y: multi-sequence list")
    }
    y <- .encodeAA(y, arity = 27, na.rm = TRUE)
  }else{
    if(is.list(y)){
      if(length(y) == 1){
        y <- y[[1]]
      }else stop("Invalid input object y: multi-sequence list")
    }
    if(mode(y) == "character") y <- match(y, colnames(x$E)) - 1
  }
  E <- if(logspace) x$E else log(x$E)
  A <- if(logspace) x$A else log(x$A)
  states <- rownames(E)
  n <- length(y)
  if(cpp){
    res <- .ViterbiH(y, A, E, DNA, AA)
  }else{
    H <- length(states) # not including Begin state
    path <- integer(n)
    V <- array(-Inf, dim = c(H, n), dimnames = list(state = states, position = 1:n))
    P <- V + NA # pointer array
    tmp <- matrix(nrow = H, ncol = H)
    if(DNA){
      for(k in 1:H) V[k, 1] <- .probDNA(y[1], E[k, ]) + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- V[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H){
          P[k, i] <- .whichmax(tmp[, k])
          V[k, i] <- tmp[P[k, i], k] + .probDNA(y[i], E[k, ])
        }
      }
    }else if(AA){
      for(k in 1:H) V[k, 1] <- .probAA(y[1], E[k, ]) + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- V[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H) {
          P[k, i] <- .whichmax(tmp[, k])
          V[k, i] <- tmp[P[k, i], k] + .probAA(y[i], E[k, ])
        }
      }
    }else{
      y <- y + 1 # convert to R's indexing style
      for(k in 1:H) V[k, 1] <- E[k, y[1]] + A[1, k + 1]
      for (i in 2:n){
        for(k in 1:H){
          for(l in 1:H){
            tmp[k, l] <- V[k, i - 1] + A[k + 1, l + 1]
          }
        }
        for(k in 1:H) {
          P[k, i] <- .whichmax(tmp[, k])
          V[k, i] <- tmp[P[k, i], k] + E[k, y[i]]
        }
      }
    }
    ak0 <- if(any(is.finite(A[-1, 1]))) A[-1, 1] else rep(0, H)
    endstate <- .whichmax(V[, n] + ak0)
    maxLL <- V[endstate, n] + ak0[endstate]
    path[n] <- endstate
    tmp <- path[n]
    for(i in n:2){
      path[i - 1] <- P[tmp, i]
      tmp <- path[i - 1]
    }
    res <- structure(list(score = maxLL,
                          path = path - 1,
                          array = V,
                          pointer = P - 1),
                     class = "DPA")
  }
  return(res)
}
################################################################################
#' @rdname Viterbi
################################################################################
Viterbi.default <- function(x, y, type = "global", d = 8, e = 2,
                            residues = NULL, S = NULL, windowspace = "all",
                            offset = 0, cpp = TRUE, ...){
  if(!(type %in% c("global", "semiglobal", "local"))) stop("invalid type")
  if(is.list(x)) if(length(x) == 1) x <- x[[1]] else stop("x is invalid")
  if(is.list(y)) if(length(y) == 1) y <- y[[1]] else stop("y is invalid")
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA){
    if(!.isDNA(y)) stop("x is a DNAbin object but y is not")
    class(x) <- class(y) <- "DNAbin"
    if(is.null(S)){
      S <- aphid::substitution$NUC.4.4
    }else{
      IUPAC <- c("A", "T", "G", "C", "S", "W", "R", "Y", "K", "M", "B", "V", "H", "D", "N")
      if(!(identical(rownames(S), IUPAC) & identical(colnames(S), IUPAC))) {
        stop("For DNAbin objects, the scoring matrix (S) should have
             15 rows and 15 columns, ordered and named according to the IUPAC
             ambiguity codes: A, T, G, C, S, W, R, Y, K, M, B, V, H, D, N")
      }
    }
    if(identical(windowspace, "WilburLipman")){
      windowspace <- .streak(.encodeDNA(x, arity = 4, na.rm = TRUE),
                            .encodeDNA(y, arity = 4, na.rm = TRUE), arity = 4)
    }else if(identical(windowspace, "all")){
      windowspace <- c(-length(x), length (y))
    }else if(length(windowspace) != 2) stop("invalid windowspace argument")
    x <- .encodeDNA(x, arity = 15, na.rm = TRUE)
    y <- .encodeDNA(y, arity = 15, na.rm = TRUE)
  }else if(AA){
    if(!.isAA(y)) stop("x is an AAbin object but y is not")
    class(x) <- class(y) <- "AAbin"
    if(is.null(S)){
      S <- aphid::substitution$MATCH
    }else{
      IUPAC <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                 "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*")
      if(!(identical(rownames(S), IUPAC) & identical(colnames(S), IUPAC))) {
        stop("For AAbin objects, the scoring matrix (S) should have
              24 rows and 24 columns, ordered and named according to the IUPAC codes:
              A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, Z, X, *")
      }
    }
    if(identical(windowspace, "WilburLipman")){
      windowspace <- .streak(.encodeAA(x, arity = 6, na.rm = TRUE),
                            .encodeAA(y, arity = 6, na.rm = TRUE), arity = 6)
    }else if(identical(windowspace, "all")){
      windowspace <- c(-length(x), length (y))
    }else if(length(windowspace) != 2) stop("invalid windowspace argument")
    if(identical(S, aphid::substitution$GONNET)){
      x <- .encodeAA(x, arity = 22, na.rm = TRUE)
      y <- .encodeAA(y, arity = 22, na.rm = TRUE)
    }else{
      x <- .encodeAA(x, arity = 24, na.rm = TRUE)
      y <- .encodeAA(y, arity = 24, na.rm = TRUE)
    }
  }else{
    if(mode(x) != mode(y)) stop("x and y have different modes")
    if(is.null(S)){
      residues <- .alphadetect(list(x, y), residues = residues)
      S <- diag(2, nrow = length(residues)) - 1
      dimnames(S) <- list(residues, residues)
    }else{
      residues <- rownames(S)
      if(is.null(residues)) stop("Scoring matrix S must have 'dimnames' attributes")
    }
    # code x and y vectors as integers in with arity = length(residues) - 1
    x  <- .encodeCH(x, residues = residues, na.rm = TRUE)
    y  <- .encodeCH(y, residues = residues, na.rm = TRUE)
    if(identical(windowspace, "WilburLipman")) {
      windowspace <- .streak(x, y, arity = length(residues))
    }else if(identical(windowspace, "all")){
      windowspace <- c(-length(x), length (y))
    }else if(length(windowspace) != 2) stop("invalid windowspace argument")
  }
  type = switch(type, "global" = 0L, "semiglobal" = 1L, "local" = 2L, stop("invalid type"))
  n <- length(x) + 1
  m <- length(y) + 1
  # initialize scoring and pointer arrays (M and P)
  M <- array(-Inf, dim = c(n, m, 3))
  P <- M + NA
  if(cpp){
    res <- .ViterbiD(x, y, type, d, e, S, windowspace, offset)
    M[, , 1] <- res$IXmatrix
    M[, , 2] <- res$MMmatrix
    M[, , 3] <- res$IYmatrix
    P[, , 1] <- res$IXpointer
    P[, , 2] <- res$MMpointer
    P[, , 3] <- res$IYpointer
    res$array <- M
    res$pointer <- P
    res$IXmatrix <- NULL
    res$MMmatrix <- NULL
    res$IYmatrix <- NULL
    res$IXpointer <- NULL
    res$MMpointer <- NULL
    res$IYpointer <- NULL
  }else{
    x <- x + 1
    y <- y + 1 # convert to R's indexing style
    M[1, 1, 2] <- 0
    if(type == 0){
      M[, 1, 1] <- c(0, seq(from = -d, to = (- d + (n - 2) * -e), by = -e))
      M[1, , 3] <- c(0, seq(from = -d, to = (- d + (m - 2) * -e), by = -e))
      P[2, 1, 1] <- P[1, 2, 1] <- 2
      P[3:n, 1, 1] <- 1
      P[1, 3:m, 3] <- 3
    }else if (type == 1){
      M[2:n, 1, 1] <- M[1, 2:m, 3] <- 0
      P[2, 1, 1] <- P[1, 2, 1] <- 2
      P[3:n, 1, 1] <- 1
      P[1, 3:m, 3] <- 3
    }else{
      M[2:n, 1, 1] <- M[1, 2:m, 3] <- 0
      P[1, , ] <- P[, 1, ]  <- 4
    }
    for(i in 2:n){
      for(j in 2:m){
        if(j - i >= windowspace[1] & j - i <= windowspace[2]){
          sij <- S[x[i - 1], y[j - 1]] + offset
          Ixcdt <- c(M[i - 1, j, 1] - e, M[i - 1, j, 2] - (d + e))#x alig to gap in y
          Mcdt <- c(M[i - 1, j - 1, 1] + sij,
                    M[i - 1, j - 1, 2] + sij,
                    M[i - 1, j - 1, 3] + sij,
                    if(type == 2) 0 else NULL)
          Iycdt <- c(-Inf, M[i, j - 1, 2] - (d + e), M[i, j - 1, 3] - e)
          Ixmax <- .whichmax(Ixcdt)
          Mmax <- .whichmax(Mcdt)
          Iymax <- .whichmax(Iycdt)
          M[i, j, 1] <- Ixcdt[Ixmax]
          M[i, j, 2] <- Mcdt[Mmax]
          M[i, j, 3] <- Iycdt[Iymax]
          P[i, j, 1] <- Ixmax
          P[i, j, 2] <- Mmax
          P[i, j, 3] <- Iymax
        }
      }
    }
    path <- c()
    progression <- matrix(nrow = 2, ncol = 0)
    if(type == 0){
      # find highest score in bottom right corner of scoring array M
      z <- c(n, m, .whichmax(M[n, m, ]))
      score <- M[z[1], z[2], z[3]]
      while(z[1] > 1 | z[2] > 1){
        path <- c(z[3], path)
        progression <- cbind(z[1:2], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
      }
      startposition <- c(1, 1)
    }else if(type == 1){
      # find highest score on bottom row or right column of scoring array M
      border <- rbind(M[n, , ], M[-n, m, ]) ### consider M[,,2] only?
      ind <- which(border == max(border), arr.ind = TRUE)
      if(nrow(ind) > 1) ind <- ind[sample(1:nrow(ind), 1),]
      z <- if(ind[1] <= m) c(n, ind) else c(ind[1] - m, m, ind[2])
      if(is.na(P[z[1], z[2], z[3]])) stop("Alignment failed, try increasing offset")
      score <- M[z[1], z[2], z[3]]
      if(z[1] < n){
        path <- rep(1, n - z[1])
        progression <- rbind((z[1] + 1):n, rep(z[2], n - z[1]))
      }else if(z[2] < m){
        path <- rep(3, m - z[2])
        progression <- rbind((z[2] + 1):m, rep(z[1], m - z[2]))
      }
      while(z[1] > 1 & z[2] > 1){
        path <- c(z[3], path)
        progression <- cbind(z[1:2], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
      }
      if(z[1] > 1){
        path <- c(rep(1, z[1] - 1), path)
        progression <- cbind(rbind(1:(z[1] - 1), rep(1, z[1] - 1)), progression)
      }else if(z[2] > 1){
        path <- c(rep(3, z[2] - 1), path)
        progression <- cbind(rbind(rep(1, z[2] - 1), 1:(z[2] - 1)), progression)
      }
      startposition <- c(1, 1)
    }else if(type == 2){
      # find highest score in scoring array M
      ind <- which(M[, , 2] == max(M[, , 2]), arr.ind = TRUE)
      z <- c(ind, 2)
      score <- M[z[1], z[2], z[3]]
      P[1, 1, 2] <- 4
      while(P[z[1], z[2], z[3]] != 4){
        path <- c(z[3], path)
        progression <- cbind(z[1:2], progression)
        z[3] <- P[z[1], z[2], z[3]]
        z <- z - switch(path[1], c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
      }
      startposition <- z[1:2]
    }
    #key <- "1: x aligns to gap in y, 2: match, 3: y aligns to gap in x"
    #progression <- progression - 1
    path <- unname(path) - 1
    #rownames(progression) <- c(deparse(substitute(x)), deparse(substitute(y)))
    res <- structure(list(score = score,
                          path = path,
                          start = startposition,
                          #progression = progression,
                          array = M,
                          pointer = P),
                     class = "DPA")
  }
  return(res)
}
################################################################################
