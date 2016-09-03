#' Generate random sequences from a model.
#'
#' \code{generate} outputs an instance of a sequence of residues randomly emitted by a
#' given hidden Markov model.
#'
#' @param x an object of class \code{'HMM'} or \code{'PHMM'}.
#' @param size a non-negative integer representing the length of the output sequence
#' (if x is a HMM with zero probability of transitioning to the BeginEnd state),
#' or the maximum length of the output sequence (acts as a safeguard against overflow).
#' @param logspace logical argument indicating whether the emission and transition
#' probabilities of x are logged (base e; TRUE) or raw (FALSE). Alternatively, if
#' \code{logspace = "autodetect"} (default), the function will automatically detect
#' if the probabilities are in log space, returning an error if inconsistencies are found.
#' Note that choosing the latter option increases the computational
#' overhead; therefore specifying \code{TRUE} or \code{FALSE} can reduce the running time.
#' @param gapchar the character used to represent gaps in the output sequence
#' (\code{generate.PHMM} only).
#' @return a named character vector giving the sequence of residues emitted by the model,
#' with "names" attribute representing the hidden states.
#' @name generate
#'
generate <- function(x, size, logspace = "autodetect", gapchar = "-"){
  UseMethod("generate")
}

#' @rdname generate
generate.HMM <- function (x, size, logspace = "autodetect"){
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  A <- if(logspace) exp(x$A) else x$A
  E <- if(logspace) exp(x$E) else x$E
  states <- rownames(A)
  residues <- colnames(E)
  hidden <- character(size)
  emitted <- character(size)
  state <- sample(states, size = 1, prob = A["BeginEnd", ])
  if(state == "BeginEnd") return(character(0))
  counter <- 1L
  while(state != "BeginEnd" & counter <= size){
    hidden[counter] <- state
    emitted[counter] <- sample(residues, size = 1, prob = E[state, ])
    state <- sample(states, size = 1, prob = A[state,])
    counter <- counter + 1L
  }
  hidden <- hidden[1:(counter - 1)]
  emitted <- emitted[1:(counter - 1)]
  names(emitted) <- hidden
  return(emitted)
}

#' @rdname generate
generate.PHMM <- function (x, size, logspace = "autodetect", gapchar = "-"){
  if(identical(logspace, "autodetect")) logspace <- logdetect(x)
  A <- if(logspace) exp(x$A) else x$A
  E <- if(logspace) exp(x$E) else x$E
  qe <- if(logspace) exp(x$qe) else x$qe
  emitted <- character(size)
  hidden <- integer(size)
  states <- c("D", "M", "I")
  position <- 0
  state <- "M"
  counter <- 1
  while(counter <= size & position <= x$size){
    state <- sample(states, size = 1, prob = A[paste0(state, states), paste(position)])
    if(state == states[2]){
      position <- position + 1
      if(position > x$size) break
      emitted[counter] <- sample(rownames(E), size = 1, prob = E[, position])
    } else if(state == states[3]){
      emitted[counter] <- sample(rownames(E), size = 1, prob = qe)
    } else{
      position <- position + 1
      emitted[counter] <- gapchar
    }
    hidden[counter] <- state
    counter <- counter + 1
  }
  emitted <- emitted[1:(counter - 1)]
  hidden <- hidden[1:(counter - 1)]
  names(emitted) <- hidden
  return(emitted)
}

