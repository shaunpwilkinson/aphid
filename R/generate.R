#' Generate random sequences from a model.
#'
#' The \code{generate} function outputs a random sequence from a HMM or PHMM.
#'
#' @param x an object of class \code{'HMM'} or \code{'PHMM'}.
#' @param size a non-negative integer representing the length of the output
#'   sequence if x is a \code{"HMM"} object with zero probability of
#'   transitioning to the begin/end state, or the maximum length of the
#'   output sequence otherwise (this acts as a safeguard against overflow).
#' @param logspace logical indicating whether the emission and transition
#'   probabilities of x are logged. If \code{logspace = "autodetect"}
#'   (the default setting), the function will automatically detect
#'   if the probabilities are logged, returning an error if
#'   inconsistencies are found. Note that choosing the latter option
#'   increases the computational overhead; therefore specifying
#'   \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param gap the character used to represent gaps (delete states)
#'   in the output sequence (only applicable for \code{PHMM} objects).
#' @param random logical indicating whether residues should be emitted randomly
#'   with probabilities defined by the emission probabilities in the model
#'   (TRUE; default), or deterministically, whereby each residue is emitted
#'   and each transition taken based on the maximum emission/transition
#'   probability in the current state.
#' @param DNA logical indicating whether the returned sequence should be a
#'   \code{"DNAbin"} object. Only applicable if the matrix of emission
#'   probabilities in the model has four residues corresponding to the nucleotide
#'   alphabet (A, T, G, and C).
#' @param AA logical indicating whether the returned sequence should be a
#'   \code{"AAbin"} object. Only applicable if the matrix of emission
#'   probabilities in the model has 20 residues corresponding to the amino acid
#'   alphabet.
#' @param ... additional arguments to be passed between methods.
#' @return a named vector giving the sequence of residues emitted by the model,
#'  with the "names" attribute representing the hidden states.
#' @details
#'   This simple function generates a single sequence
#'   from a HMM or profile HMM by recursively simulating a path through
#'   the model. The function is fairly slow in its current state, but a
#'   faster C++ function may be made available in a future version depending
#'   on demand.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @examples
#'   ## Generate a random sequence from a standard HMM
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
#'   ### Generate a random sequence from the model
#'   generate(x, size = 300)
#'   ##
#'   ## Generate a random sequence from a profile HMM:
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
#'   ### Names attribute indicates that all residues came from "match" states
#' @name generate
################################################################################
generate <- function(x, size, ...){
  UseMethod("generate")
}
################################################################################
#' @rdname generate
################################################################################
generate.HMM <- function (x, size, logspace = "autodetect", random = TRUE, ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  A <- if(logspace) exp(x$A) else x$A
  E <- if(logspace) exp(x$E) else x$E
  states <- rownames(A)
  residues <- colnames(E)
  hidden <- character(size)
  emitted <- character(size)
  if(random){
    state <- sample(states, size = 1, prob = A["Begin", ])
  }else{
    state = states[.whichmax(A["Begin", ])]
  }
  if(state == "Begin") return(character(0))
  counter <- 1L
  while(state != "Begin" & counter <= size){
    hidden[counter] <- state
    if(random){
      emitted[counter] <- sample(residues, size = 1, prob = E[state, ])
      state <- sample(states, size = 1, prob = A[state,])
    }else{
      emitted[counter] <- residues[.whichmax(E[state, ])]
      state <- states[.whichmax(A[state,])]
    }
    counter <- counter + 1L
  }
  hidden <- hidden[1:(counter - 1)]
  emitted <- emitted[1:(counter - 1)]
  names(emitted) <- hidden
  return(emitted)
}
################################################################################
#' @rdname generate
################################################################################
generate.PHMM <- function (x, size, logspace = "autodetect", gap = "-",
                           random = TRUE, DNA = FALSE, AA = FALSE, ...){
  if(identical(logspace, "autodetect")) logspace <- .logdetect(x)
  A <- if(logspace) exp(x$A) else x$A
  E <- if(logspace) exp(x$E) else x$E
  qe <- if(is.null(x$qe)) rep(1/nrow(E), nrow(E)) else if(logspace) exp(x$qe) else x$qe
  #### condition if names qe and rownames E mismatch
  stopifnot(!(DNA & AA))
  if(DNA) gap <- as.raw(4) else if(AA) gap <- as.raw(45)
  emitted <- if(DNA | AA) raw(size) else character(size)
  hidden <- integer(size)
  states <- c("D", "M", "I")
  if(DNA) rownames(E)[toupper(rownames(E)) == "U"] <- "T"
  residues <- if(DNA){
    as.raw(c(136, 40, 72, 24))[sapply(toupper(rownames(E)), match, c("A", "C", "G", "T"))]
  }else if(AA){
    as.raw(65:89)[-c(2, 10, 15, 21, 24)]
  }else rownames(E)
  position <- 0
  state <- "M"
  counter <- 1
  while(counter <= size & position <= x$size){
    if(random){
      state <- sample(states, size = 1, prob = A[paste0(state, states), paste(position)])
    }else{
      state <- states[.whichmax(A[paste0(state, states), paste(position)])]
    }
    if(state == states[2]){
      position <- position + 1
      if(position > x$size) break
      if(random){
        emitted[counter] <- sample(residues, size = 1, prob = E[, position])
      }else{
        emitted[counter] <- residues[.whichmax(E[, position])]
      }
    }else if(state == states[3]){
      if(random){
        emitted[counter] <- sample(residues, size = 1, prob = qe)
      }else{
        emitted[counter] <- residues[.whichmax(qe)]
      }
    }else{
      position <- position + 1
      emitted[counter] <- gap
    }
    hidden[counter] <- state
    counter <- counter + 1
  }
  emitted <- emitted[1:(counter - 1)]
  hidden <- hidden[1:(counter - 1)]
  names(emitted) <- hidden
  class(emitted) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
  return(emitted)
}
################################################################################
