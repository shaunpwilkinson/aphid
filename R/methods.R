# ViterbiC <- function(x, y, qe = NULL, logspace = FALSE,
#                     type = "semiglobal", offset = -0.1, d = 8,
#                     e = 2, S = NULL, itertab = NULL){
#
#   UseMethod("ViterbiC")
# }

kdist <- function(x, k = 5, asmatrix = FALSE){
  UseMethod("kdist")
}

