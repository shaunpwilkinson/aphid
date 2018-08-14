#' Export profile hidden Markov models as text.
#'
#' \code{writePHMM} takes an object of class \code{"PHMM"} and writes it to a
#'   text file in HMMER3 format.
#'
#' @param x an object of class \code{"PHMM"}.
#' @param file the name of the file to write the model to.
#' @param append logical indicating whether the model text should be appended
#'   below any existing text in the output file, or whether any existing text
#'   should be overwritten. Defaults to FALSE.
#' @param form character string indicating the format in which to write the model.
#'   Currently only HMMER3f is supported.
#' @param vers character string indicating the version of version of the format
#'   in which to write the model. Currently only "f" is supported.
#' @return NULL (invisibly)
#' @details
#'   This function writes an object of class \code{"PHMM"} to a
#'   HMMER3/f text file. Note that unlike HMMER, the \pkg{aphid}
#'   package does not currently support position-specific background
#'   emission probabilities.
#' @author Shaun Wilkinson
#' @references
#'   Finn, RD, Clements J & Eddy SR (2011) HMMER web server: interactive sequence
#'   similarity searching. \emph{Nucleic Acids Research}. \strong{39} W29-W37.
#'   \url{http://hmmer.org/}.
#'
#'   HMMER: biosequence analysis using profile hidden Markov models.
#'   \url{http://www.hmmer.org}.
#'
#' @seealso \code{\link{readPHMM}} to parse a PHMM object from a HMMER3 text file.
#' @examples
#'   ## Derive a profile hidden Markov model from the small globin alignment
#'   data(globins)
#'   x <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
#'   x
#'   fl <- tempfile()
#'   writePHMM(x, file = fl)
#'   readPHMM(fl)
#'   ##
#'   ## Derive a PHMM for the woodmouse data and write to file
#'   \donttest{
#'     library(ape)
#'     data(woodmouse)
#'     woodmouse.PHMM <- derivePHMM(woodmouse)
#'     tmpf <- tempfile(fileext = ".hmm")
#'     writePHMM(woodmouse.PHMM, file = tmpf)
#'   }
################################################################################
writePHMM <- function(x, file = "", append = FALSE, form = "HMMER3", vers = "f"){
  options(digits = 7, scipen = 10)
  stopifnot(form == "HMMER3" & vers == "f")
  if(!(inherits(x, "PHMM"))) stop("Input object must be of class 'PHMM'")
  cat(paste0(form, "/", vers), file = file, sep = "\n", append = append)
  ##[produced by R::profile?]
  cat(paste0("NAME  ", x$name), file = file, sep = "\n", append = TRUE)
  if(!is.null(x$accession)) cat(paste0("ACC   ", x$accession), file = file,
                                sep = "\n", append = TRUE)
  if(!is.null(x$description)) cat(paste0("DESC  ", x$description), file = file,
                                  sep = "\n", append = TRUE)
  cat(paste0("LENG  ", x$size), file = file, sep = "\n", append = TRUE)
  if(!is.null(x$maxlength)) cat(paste0("MAXL  ", x$maxlength), file = file,
                                sep = "\n", append = TRUE)
  cat(paste0("ALPH  ", x$alphabet), file = file, sep = "\n", append = TRUE)

  rfl <- if(is.null(x$reference)) "RF    no" else "RF    yes"
  mml <- if(is.null(x$mask))      "MM    no" else "MM    yes"
  cl <- if(is.null(x$consensus))  "CONS  no" else "CONS  yes"
  csl <- if(is.null(x$construct)) "CS    no" else "CS    yes"
  mpl <- if(is.null(x$map)) "MAP   no" else "MAP   yes"
  cat(rfl, mml, cl, csl, mpl, file = file, sep = "\n", append = TRUE)
  if(!is.null(x$date)) cat(paste0("DATE  ", x$date), file = file,
                           append = TRUE, sep = "\n")
  if(!is.null(x$nseq)) cat(paste0("NSEQ  ", x$nseq), file = file,
                           append = TRUE, sep = "\n")
  if(!is.null(x$effn)) cat(paste0("EFFN  ", round(x$effn, 6)), file = file,
                           append = TRUE, sep = "\n")
  if(!is.null(x$checksum)) cat(paste0("CKSUM ", x$checksum), file = file,
                               append = TRUE, sep = "\n")
  ### placeholder for stats, cutoffs, etc
  residues <- rownames(x$E)
  cat(c("HMM          ", paste0(residues, "        "), "\n"), file = file,
      append = TRUE, sep = "")
  cat("            m->m     m->i     m->d     i->m     i->i     d->m     d->d",
      file = file, append = TRUE, sep = "\n")
  if(!is.null(x$compo)){
    compoline <- formatC(-x$compo, digits = 5, format = "f")
    compoline[trimws(compoline) == "Inf"] <- "      *"
    compoline <- gsub("(.......).+", "\\1", compoline)
    # in case any logprobs are > 10
    compoline = paste0("  COMPO   ", paste(compoline, collapse = "  "))
    cat(compoline, file = file, sep = "\n", append = TRUE)
  }
  qe <- formatC(-x$qe, digits = 5, format = "f")
  qe[trimws(qe) == "Inf"] <- "      *"
  qe <- gsub("(.......).+", "\\1", qe)
  qeline = paste0("          ", paste(qe, collapse = "  "))
  cat(qeline, file = file, sep = "\n", append = TRUE)
  A <- formatC(-x$A[c(5, 6, 4, 8, 9, 2, 1), ], digits = 6, format = "f")
  A[trimws(A) == "Inf"] <- "      *"
  A <- gsub("(.......).+", "\\1", A)
  dim(A) <- c(7, x$size + 1)
  A[6, 1] <- A[6, x$size + 1] <- "0.00000"
  E <- formatC(-x$E, digits = 6, format = "f")
  E[trimws(E) == "Inf"] <- "      *"
  E <- gsub("(.......).+", "\\1", E)
  dim(E) <- c(length(residues), x$size)
  trans0line = paste("         ", paste(A[, 1], collapse = "  "))
  cat(trans0line, file = file, sep = "\n", append = TRUE)
  for(i in 1:x$size){
    emline <- paste0(if(i < 10) "      " else if(i < 100) "     " else "    ",
                     i, "   ")
    emline <- paste0(emline, paste(E[, i], collapse = "  "))
    if(!is.null(x$map)){
      tmp <- x$map[i]
      emline <- paste0(emline, if(tmp < 10) "      " else if(tmp < 100) "     " else  "    ")
      emline <- paste0(emline, tmp)
    }else{
      emline <- paste0(emline, "      -")
    }
    emline <- paste0(emline, " ", if(!is.null(x$consensus)) x$consensus[i] else "-")
    emline <- paste0(emline, " ", if(!is.null(x$reference)) x$reference[i] else "-")
    emline <- paste0(emline, " ", if(!is.null(x$mask)) x$mask[i] else "-")
    emline <- paste0(emline, " ", if(!is.null(x$construct)) x$construct[i] else "-")
    cat(emline, file = file, sep = "\n", append = TRUE)
    cat(qeline, file = file, sep = "\n", append = TRUE)
    transline <- paste0("          ", paste(A[, i + 1], collapse = "  "))
    cat(transline, file = file, sep = "\n", append = TRUE)
  }
  cat("//", file = file, sep = "", append = TRUE)
}
################################################################################
