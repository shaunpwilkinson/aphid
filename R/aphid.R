#' The \pkg{aphid} package for analysis with profile hidden Markov models.
#'
#' The \pkg{aphid} package contains functions for building and using profile hidden
#'   Markov models for biological sequence analysis. Functions are included for
#'   multiple and pairwise sequence alignment, model construction and parameter optimization,
#'   calculation of conditional probabilities (with the forward and Viterbi algorithms),
#'   tree-based sequence weighting, k-mer distance calculation, sequence simulation,
#'   and file import/export compatible with the \href{http://www.hmmer.org}{HMMER}
#'   software package.
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
#'   and "AAbin" object types produced by the \code{\link[ape]{ape}} package.
#'   This is an essential piece of software for those
#'   using R for biological sequence analysis, since it provides a binary coding format
#'   for nucleotides and amino acids that maximizes memory and speed efficiency. While
#'   \pkg{aphid} also works with standard character vectors and matrices, it will
#'   not recognize the DNA and amino acid amibguity codes and therefore is not
#'   guaranteed to treat them appropriately.
#'
#'   To maximize speed, the low-level dynamic programming functions such
#'   as \code{\link{Viterbi}}, \code{\link{forward}} and \code{\link{backward}}
#'   are written in C++ with the help of the \code{\link[Rcpp]{Rcpp}}
#'   package. Thus those wishing to build the package from source will need a C/C++
#'   compiler such as \href{https://clang.llvm.org/}{clang}.
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
#'   Objects of class \code{"Viterbi"} and \code{"fullprob"} are also generated
#'   by the Viterbi and forward/backward functions, respectively. Functions in the
#'   \pkg{aphid} package are compatible with "DNAbin" and "AAbin" objects generated
#'   by the \code{\link[ape]{ape}} package.
#'   These object types, in which sequences are represented in a bit-level
#'   coding scheme, are recommended for maximizing memory and speed effficiency.
#'
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
#' alignment or a list of unaligned sequences
#' \item \code{\link{map}} optimizes profile hidden Markov model construction
#' using the maximum \emph{a posteriori} algorithm
#' \item \code{\link{train.HMM}} optimizes the parameters of a \code{"HMM"} object
#' using a list of training sequences
#' \item \code{\link{train.PHMM}} optimizes the parameters of a \code{"PHMM"} object
#' using a list of training sequences
#' }
#'
#' @section Sequence alignment and weighting:
#' \itemize{
#' \item \code{\link{align}} performs a PHMM-based progressive multiple sequence
#' alignment
#' \item \code{\link{weight}} assigns weights to sequences based on a tree
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
#' \item \code{\link{generate}} simulates random sequences from a HMM
#' or PHMM
#' }
#'
#' @section Datasets:
#' \itemize{
#' \item \code{\link{substitution}} DNA and protein substitution
#' matrices including PAM, BLOSUM, GONNET, DAYHOFF, NUC.4.2 and NUC.4.4
#' \item \code{\link{casino}} the dishonest casino example from Durbin et al.
#' (1998) chapter 3.2
#' \item \code{\link{globins}} Small globin alignment data from
#' Durbin et al. (1998) Figure 5.3
#' }
#'
#' @author Shaun Wilkinson
#'
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#'
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Edgar RC (2004) Local homology recognition and distance measures in
#'   linear time using compressed amino acid alphabets.
#'   \emph{Nucleic Acids Research}, \strong{32}, 380-385.
#'
#'   Finn, RD, Clements J & Eddy SR (2011) HMMER web server: interactive sequence
#'   similarity searching.
#'   \emph{Nucleic Acids Research}. \strong{39}:W29-W37. \url{http://hmmer.org/}.
#'
#'   Gerstein M, Sonnhammer ELL, Chothia C (1994) Volume changes in protein evolution.
#'   \emph{Journal of Molecular Biology}, \strong{236}, 1067-1078.
#'
#'   HMMER: biosequence analysis using profile hidden Markov models.
#'   \url{http://www.hmmer.org}.
#'
#'   NCBI index of substitution matrices.
#'   \url{ftp://ftp.ncbi.nih.gov/blast/matrices/}.
#'
#'   Juang B-H, Rabiner LR (1990) The segmental K-means
#'   algorithm for estimating parameters of hidden Markov models.
#'   \emph{IEEE Transactions on Acoustics, Speech, and Signal Processing},
#'   \strong{38}, 1639-1641.
#'
#'   Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H,
#'   Remmert M, Soding J, Thompson JD, Higgins DG (2011) Fast, scalable generation
#'   of high-quality protein multiple sequence alignments using Clustal Omega.
#'   \emph{Molecular Systems Biology}, \strong{7}, 539.
#'
#'   Soding J (2005) Protein homology detection by HMM-HMM comparison.
#'   \emph{Bioinformatics}, \strong{21}, 951-960.
#'
#'   Wilbur WJ, Lipman DJ (1983) Rapid similarity searches of nucleic acid and
#'   protein data banks. \emph{Proc Natl Acad Sci USA}, \strong{10}, 197-206.
#'
#'   Yang K, Zhang L (2008) Performance comparison between k-tuple distance
#'   and four model-based distances in phylogenetic tree reconstruction.
#'   \emph{Nucleic Acids Research}, \strong{36}, e33.
#'
#' @docType package
#' @name aphid
################################################################################
NULL

