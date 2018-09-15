#' The \pkg{aphid} package for analysis with profile hidden Markov models.
#'
#' \pkg{aphid} is an R package for the development and application of
#'   hidden Markov models and profile HMMs for biological sequence analysis.
#'   Functions are included for multiple and pairwise sequence alignment,
#'   model construction and parameter optimization, calculation of conditional
#'   probabilities (using the forward, backward and Viterbi algorithms),
#'   tree-based sequence weighting, sequence simulation, and file import/export
#'   compatible with the \href{http://www.hmmer.org}{HMMER} software package.
#'   The package has a wide variety of uses including database searching,
#'   gene-finding and annotation, phylogenetic analysis and sequence classification.
#'
#' @details
#' The \pkg{aphid} package is based on the algorithms outlined in the book 'Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids' by Richard Durbin,
#'   Sean Eddy, Anders Krogh and Graeme Mitchison. This book is highly recommended
#'   for those wishing to develop a better understanding of HMMs and PHMMs, regardless of
#'   prior experience. Many of the examples in the function help pages are taken directly
#'   from the book, so that readers can learn to use the package as they work through the
#'   chapters.
#'
#'   There are also excellent rescources available for those wishing to use profile hidden
#'   Markov models outside of the R environment. The \pkg{aphid} package maintains
#'   compatibility with the \href{http://www.hmmer.org}{HMMER} software suite
#'   through the file input and output functions \code{\link{readPHMM}} and
#'   \code{\link{writePHMM}}. Those interested are further encouraged to check out the
#'   \href{https://compbio.soe.ucsc.edu/sam.html}{SAM} software package,
#'   which also features a comprehensive suite of functions and tutorials.
#'
#'   The \pkg{aphid} package is designed to work in conjunction with the "DNAbin"
#'   and "AAbin" object types produced by the \code{\link[ape]{ape}} package
#'   (Paradis et al 2004, 2012). This is an essential piece of software for those
#'   using R for biological sequence analysis, and provides a binary coding format
#'   for nucleotides and amino acids that maximizes memory and speed efficiency.
#'   While \pkg{aphid} also works with standard character vectors and matrices,
#'   it may not recognize the DNA and amino acid amibguity codes and therefore is not
#'   guaranteed to treat them appropriately.
#'
#'   To maximize speed, the low-level dynamic programming functions such
#'   as \code{\link{Viterbi}}, \code{\link{forward}} and \code{\link{backward}}
#'   are written in C++ with the help of the \code{\link[Rcpp]{Rcpp}}
#'   package (Eddelbuettel & Francois 2011).
#'   Note that R versions of these functions are also maintained
#'   for the purposes of debugging, experimentation and code interpretation.
#'
#' @section Classes:
#'   The \pkg{aphid} package creates two primary object classes, \code{"HMM"}
#'   (hidden Markov models) and \code{"PHMM"} (profile hidden Markov models)
#'   with the functions \code{\link{deriveHMM}} and \code{\link{derivePHMM}}, respectively.
#'   These objects are lists consisting of emission and transition probability matrices
#'   (denoted E and A), vectors of non-position-specific background emission and transition
#'   probabilies (denoted qe and qa) and other model metadata.
#'   Objects of class \code{"DPA"} (dynammic programming array) are also generated
#'   by the Viterbi and forward/backward functions.
#'   These are primarily created for succinct console printing.
#'
#' @section Functions:
#' A breif description of the primary \pkg{aphid} functions are provided with links
#' to their help pages below.
#'
#' @section File import and export:
#' \itemize{
#' \item \code{\link{readPHMM}} parses a \href{http://www.hmmer.org}{HMMER} text file
#'   into R and creates an object of class \code{"PHMM"}
#' \item \code{\link{writePHMM}} writes a \code{"PHMM"} object to a text file in
#'   \href{http://www.hmmer.org}{HMMER} v3 format
#' }
#'
#' @section Visualization:
#' \itemize{
#' \item \code{\link{plot.HMM}} plots a \code{"PHMM"} object as a cyclic directed graph
#' \item \code{\link{plot.PHMM}} plots a \code{"PHMM"} object as a directed graph with
#'  sequential modules consisting of match, insert and delete states
#' }
#'
#' @section Model building and training:
#' \itemize{
#' \item \code{\link{deriveHMM}} builds a \code{"HMM"} object from a list of training
#' sequences
#' \item \code{\link{derivePHMM}} builds a \code{"PHMM"} object from a multiple sequence
#' alignment or a list of non-aligned sequences
#' \item \code{\link{map}} optimizes profile hidden Markov model construction
#' using the maximum \emph{a posteriori} algorithm
#' \item \code{\link{train}} optimizes the parameters of a \code{"HMM"} or
#' \code{"PHMM"} object using a list of training sequences
#' }
#'
#' @section Sequence alignment and weighting:
#' \itemize{
#' \item \code{\link{align}} performs a multiple sequence alignment
#' \item \code{\link{weight}} assigns weights to sequences
#' }
#'
#' @section Conditional probabilities:
#' \itemize{
#' \item \code{\link{Viterbi}} finds the optimal path of a sequence through a HMM
#' or PHMM, and returns its log odds or probability given the model
#' \item \code{\link{forward}} finds the full probability of a sequence
#' given a HMM or PHMM using the forward algorithm
#' \item \code{\link{backward}} finds the full probability of a sequence
#' given a HMM or PHMM using the backward algorithm
#' \item \code{\link{posterior}} finds the position-specific posterior probability
#' of a sequence given a HMM or PHMM
#' }
#'
#' @section Sequence simulation:
#' \itemize{
#' \item \code{\link{generate.HMM}} simulates a random sequence from an HMM
#' \item \code{\link{generate.PHMM}} simulates a random sequence from a PHMM
#' }
#'
#' @section Datasets:
#' \itemize{
#' \item \code{\link{substitution}} a collection of DNA and amino acid
#' substitution matrices from \href{ftp://ftp.ncbi.nih.gov/blast/matrices/}{NCBI}
#' including the PAM, BLOSUM, GONNET, DAYHOFF and NUC matrices
#' \item \code{\link{casino}} data from the dishonest casino example of
#' Durbin et al (1998) chapter 3.2
#' \item \code{\link{globins}} Small globin alignment data from
#' Durbin et al (1998) Figure 5.3
#' }
#'
#' @author Shaun Wilkinson
#'
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Eddelbuettel D, Francois R (2011) Rcpp: seamless R and C++ integration.
#'   \emph{Journal of Statistical Software} \strong{40}, 1-18.
#'
#'   Finn RD, Clements J & Eddy SR (2011) HMMER web server: interactive sequence
#'   similarity searching.
#'   \emph{Nucleic Acids Research}. \strong{39}, W29-W37. \url{http://hmmer.org/}.
#'
#'   HMMER: biosequence analysis using profile hidden Markov models.
#'   \url{http://www.hmmer.org}.
#'
#'   NCBI index of substitution matrices.
#'   \url{ftp://ftp.ncbi.nih.gov/blast/matrices/}.
#'
#'   Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
#'   and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.
#'
#'   Paradis E (2012) Analysis of Phylogenetics and Evolution with R
#'   (Second Edition). Springer, New York.
#'
#' @docType package
#' @name aphid
################################################################################
NULL

