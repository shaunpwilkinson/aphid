#' Substutution matrices.
#'
#' A dataset containing several of the more popular substitution scoring matrices
#' for DNA and amino acids.
#'
#' @format A list of 71 matrices, most of which have 24 rows and 24 columns
#' corresponding to the amino acid alphabet plus the ambiguity codes B, Z,
#' X and *:
#' \describe{
#'   \item{PAM}{the PAM matrices from PAM10 to PAM500.}
#'   \item{BLOSUM}{the BLOSUM matrices from BLOSUM30 to BLOSUM100.}
#'   \item{others}{also included are the DAYHOFF, GONNET, IDENTITY
#'   and MATCH substitution matrices for proteins
#'   and the NUC.4.2 and NUC.4.4 substitution matrices for DNA.}
#' }
#' @source \url{ftp://ftp.ncbi.nih.gov/blast/matrices/}
"substitution_matrices"



#' Dishonest casino example.
#'
#' The 'dishonest casino' example from Durbin et al. (1998) chapter 3.2.
#'
#' @format A named character vector showing the result of 300 rolls of a die
#'   that switches from "Fair" to "Loaded" with a probability
#'   of 0.05 and back to "Fair" with a probability of 0.1. In the Fair
#'   state each outcome from 1 to 6 has an equal probability of occurring,
#'   while in the Loaded state the probability of rolling a "6" increases
#'   to 0.5 (with the remaining five probabilities reduced to 0.1).
#'   The elements of the vector are the outcomes of the 300 rolls
#'   ("1", "2", "3", "4", "5", or "6") and the "names" attribute
#'   represents the underlying Markov states ("Fair" or "Loaded").
#' @source
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
"casino"


#' Globin protein alignment.
#'
#' The small globin protein alignment from figure 5.3 of
#' Durbin et al. (1998).
#'
#' @format a 7 x 10 character matrix showing ten columns from a multiple
#'   alignment of globin proteins shown in Durbin et al. (1998) chapter
#'   5.3.
#' @source
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
"globins"
