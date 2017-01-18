#' Tree-based sequence weighting.
#'
#' Takes a dendogram object and returns a vector of sequence weights following
#' the method of
#' Gerstein et al. (1994).
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param method a character string indicating the weighting method to be used.
#'     Currently only that of Gerstein et al. (1994) is supported
#'     (\code{method = "Gerstein"}).
#' @return a named vector of weights, the sum of which is equal to
#'    the total number of sequences.
#' @examples
#' data(woodmouse)
#' wood_dist <- kdistance(unalign(woodmouse))
#' wood_den <- as.dendrogram(hclust(wood_dist, method = "average"))
#' wood_weights <- weight(wood_den, method = "Gerstein")
#' @name weight
#' @export
#'
weight <- function(x, method = "Gerstein", k = 5, residues = NULL, gapchar = "-"){
  UseMethod("weight")
}


#' @rdname weight
#' @export
#'
weight.DNAbin <- function(x, method = "Gerstein", k = 5){
  if(is.list(x)){
    weight.list(x, method = method, k = k)
  }else{
    x <- unalign(x, gapchar = as.raw(4))
    weight.list(x, method = method, k = k)
  }
}


#' @rdname weight
#' @export
#'
weight.AAbin <- function(x, method = "Gerstein", k = 5){
  if(is.list(x)){
    weight.list(x, method = method, k = k)
  }else{
    x <- unalign(x, gapchar = as.raw(45))
    weight.list(x, method = method, k = k)
  }
}


#' @rdname weight
#' @export
#'
weight.list <- function(x, method = "Gerstein", k = 5, residues = NULL, gapchar = "-"){
  nsq <- length(x)
  DNA <- is.DNA(x)
  AA <- is.AA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- alphadetect(x, residues = residues, gapchar = gapchar)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  for(i in 1:nsq) x[[i]] <- x[[i]][x[[i]] != gapchar]
  # cache names for later
  tmpnames <- names(x)
  names(x) <- paste0("S", 1:nsq)
  if(nsq > 2){
    qds <- kdistance(x, k = k, alpha = if(AA) "Dayhoff6" else if(DNA) NULL else residues)
    guidetree <- as.dendrogram(hclust(qds, method = "average"))
    res <- weight.dendrogram(guidetree, method = "Gerstein")[names(x)]
  }else if(nsq == 2){
    res <- c(0.5, 0.5)
  }else if(nsq == 1){
    res <- 1
  }else{
    res <- numeric(0)
  }
  names(res) <- tmpnames
  return(res)
}


#' @rdname weight
#' @export
#'
weight.dendrogram <- function(x, method = "Gerstein"){
  if(!identical(method, "Gerstein")) stop("Only Gerstein et al. 1994 method supported")
  acal <- function(d) !any(sapply(d, is.list)) # all children are leaves?
  md <- function(d) all(sapply(d, acal)) & !acal(d) # mergable dendro?
  ch <- function(d) sapply(d, attr, "height") # child heights
  Gerstein <- function(x){ # x is a dendrogram
    ngrandchildren <- sapply(x, length)
    childisdendro <- ngrandchildren > 1
    if(md(x)){
      childheights <- ch(x[childisdendro])
      childedges <- attr(x, "height") - childheights # ch only works on dendro lists
      grandchildheights <- lapply(x[childisdendro], ch) #  list same length as childedges
      grandchildedges <- mapply("-", childheights, grandchildheights, SIMPLIFY = FALSE)
      grandchildedges <- lapply(grandchildedges, function(e) e + 0.0000001) # this just
      # safeguards against 0 denominators (but is a bit of a hack)
      ratios <- lapply(grandchildedges, function(v) v/sum(v))
      inheritances <- mapply("*", childedges, ratios, SIMPLIFY = FALSE)
      newgrandchildedges <- mapply("+", grandchildedges, inheritances, SIMPLIFY = FALSE)
      lcounter <- 1 #leaf counter
      dcounter <- 1 #dendro counter
      tmp <- x
      for(i in seq_along(ngrandchildren)){
        if(ngrandchildren[i] > 1){
          for(j in 1:ngrandchildren[i]){
            leafj <- tmp[[i]][[j]]
            attr(leafj, "height") <- attr(tmp, "height") - newgrandchildedges[[dcounter]][j]
            x[[lcounter]] <- leafj
            lcounter <- lcounter + 1
          }
          dcounter <- dcounter + 1
        }else{
          x[[lcounter]] <- tmp[[i]]
          lcounter <- lcounter + 1
        }
      }
    }else{
      x[childisdendro] <- lapply(x[childisdendro], Gerstein)
    }
    x
  }
  while(any(sapply(x, is.list))) x <- Gerstein(x)
  res <- sapply(x, function(d) attr(x, "height") - attr(d, "height"))
  res <- res * length(res)/sum(res)
  names(res) <- sapply(x, function(d) attr(d, "label"))
  res
}


#' @rdname weight
#' @export
#'
weight.default <- function(x, method = "Gerstein", k = 5, residues = NULL, gapchar = "-"){
  x <- unalign(x, gapchar = gapchar)
  weight.list(x, method = method, k = k, residues = residues, gapchar = gapchar)
}

