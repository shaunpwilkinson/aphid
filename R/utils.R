#' Utilities.
#'
# x is an integer vector in 'from' numbering system (eg for binary, from = 2)
# decimal(x, from) is the same as convert(x, from, to = 10, collapse = FALSE)
decimal <- function(x, from) sum(x * from^rev(seq_along(x) - 1))

convert <- function(x, from = 10, to = 2, collapse = FALSE){
  if(to %% 1 > 0) stop("Non-integers are not supported yet")
  if(length(x) == 1) x <- as.integer(strsplit(paste(x), split = "")[[1]])
  x <- sum(x * from^rev(seq_along(x) - 1))
  if(to == 10){
    if(collapse) x else as.integer(strsplit(paste(x), split = "")[[1]])
  }
  dividend <- as.integer(to)
  quotient <- floor(x/dividend)
  result <- x %% dividend
  while(quotient > 0){
    remainder <- quotient %% dividend
    result = c(remainder, result)
    quotient <- floor(quotient/dividend)
  }
  if(collapse) result <- as.integer(paste(result, collapse = ""))
  return(result)
}

tab <- function(v, residues){
  res <- structure(integer(length(residues)), names = residues)
  for(i in residues) res[i] <- sum(v == i, na.rm = TRUE)
  return(res)
}


whichismax <- function(v){
  ind <- which(v == max(v, na.rm = TRUE))
  if(length(ind) > 1) ind <- sample(ind, 1)
  ind
}



