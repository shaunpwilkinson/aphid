#' Tree-based sequence weighting.
#'
#' Takes a dendogram object and returns a vector of sequence weights following
#' the method of
#' Gerstein et al. (1994).
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param method a character string indicating the weighting method to be used.
#' Currently only that of Gerstein et al. (1994) is supported
#' (\code{method = "Gerstein"}).
#' @return a named vector of weights, the sum of which is equal to
#' the total number of sequences.
#' @examples
#' data(woodmouse)
#' wood_dist <- kdistance(unalign(woodmouse))
#' wood_den <- as.dendrogram(hclust(wood_dist, method = "average"))
#' wood_weights <- weight(wood_den, method = "Gerstein")
#'
weight <- function(x, method = "Gerstein"){
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
      # safeguards against 0 denominators (but its a bit of a hack)
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



