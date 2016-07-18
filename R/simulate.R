#' Simulate the evolution of DNA
#' 
#' \code{basesim}, \code{seqsim}, and \code{treesim} simulate the evolution of DNA for a single base, 
#' a single sequence, and a set of sequences given a phylogenetic tree, respectively.
#' 
#' @param x a vector of mode "character" containing the values "a", "c", "g", and "t". 
#' For \code{basesim}, the length of x must be one.
#' 
#' @param K the mutation rate in number of mutations per base pair per generation. 
#' Must be between 0 and 1.
#' 
#' @param n the number of generations to run the simulation for.
#' 
#' @param d an object of class \code{"dendrogram"}.
#' 
#' @return an object of class \code{"dendrogram"} with a 'sequence' attribute at each node.
#' 
#' @examples 
#' set.seed(999)
#' myseq <- randomSequence(20)
#' mynewick <- "(A:1,(B:0.5,C:0.2):0.6);"
#' mydendrogram <- read.dendrogram(text = mynewick)
#' mysimulation <- treesim(myseq, K = 0.01, n = 100, d = mydendrogram)
#' mysimmatrix <- treesim.as.matrix(mysimulation)
#' myDNAbin <- ape::as.DNAbin(mysimmatrix)
#' mydistmat <- ape::dist.dna(myDNAbin, model = "JC69")
#' myhc <- hclust(mydistmat, "ave")
#' plot(as.dendrogram(myhc))
#' 
basesim <- function(x, K, alphabet = c("a", "c", "g", "t")){
  if(length(x) != 1) stop("the input sequence x must be of length one")
  if((K < 0)|(K > 1)) stop("the mutation rate K must be between zero and one")
  if(!(x %in% alphabet)) stop("the initial state x must be in the alphabet")
  sample(c(x, alphabet[alphabet != x]), 1, 
         prob = c(1 - K, rep(K/(length(alphabet) - 1), length(alphabet) - 1)))
}

seqsim <- function (x, K, n, alphabet = c("a", "c", "g", "t")){ 
  # x is a sequence vector of mode character 
  # K is the mutation rate (mutations per basepair per generation)
  # n is the number of generations to run the simulation for
  if(!(all(x %in% alphabet))) 
    stop("all characters in the initial sequence x must be in the alphabet")
  n <- as.integer(n)
  stopifnot(n > 0)
  for (i in 1:n) {
    x <- sapply(x, basesim, K = K, alphabet = alphabet)
  }
  unname(x)
}

treesim <- function(x, K, n, d, alphabet = c("a", "c", "g", "t")){
  if (!(inherits(d, "dendrogram"))) 
    stop("input object d must be of class 'dendrogram'")
  max.height <- attr(d, 'height')
  min.height <- min(unlist(dendrapply(d, attr, 'height')))
  scale.factor <- max.height - min.height
  d <- dendrapply(d, function(y){
    attr(y, 'height') <- as.integer(attr(y, 'height') * n/scale.factor)
    y
  })
  attr(d, 'sequence') <- as.vector(x, mode = mode(alphabet))
  evolver <- function(y){
    if(is.list(y)){
      y[] <- lapply(y, function(z){
        edge.length <- as.integer(attr(y, 'height') - attr(z, 'height'))
        attr(z, 'sequence') <- seqsim(attr(y, 'sequence'), 
                                      K = K, n = edge.length, alphabet = alphabet)
        z
      })
      y[] <- lapply(y, evolver)
    }
    y
  }
  d <- evolver(d)
  return(d)
}

treesim.as.matrix <- function(x){ # an object of class "dendrogram" with 'sequence' attributes
  taxa <- unlist(dendrapply(x, attr, 'label'))
  res <- matrix(unlist(dendrapply(x, function(y) if(is.leaf(y)) attr(y, 'sequence'))),
    nrow = length(taxa), byrow = TRUE, dimnames = list(taxa))
  res
}

randomSequence <- function(m, alphabet = c("a", "c", "g", "t")){
  res <- rep(alphabet[1], m)
  res <- sapply(res, basesim, K = (length(alphabet) - 1)/length(alphabet), 
                alphabet = alphabet)
  unname(res)
} 


seqsim2 <- function (x, K, n, alphabet = c("a", "c", "g", "t")){ 
  # x is a sequence vector of mode character 
  # K is the mutation rate (mutations per basepair per generation)
  # n is the number of generations to run the simulation for
  n <- as.integer(n)
  stopifnot(n > 0)
  number.x1.mutations <- 0
  number.x2.mutations <- 0
  x1 <- x
  x2 <- x
  for (i in 1:n) {
    x1old <- x1
    x2old <- x2
    x1 <- sapply(x1, basesim, K = K, alphabet = alphabet)
    x2 <- sapply(x2, basesim, K = K, alphabet = alphabet)
    number.x1.mutations <- number.x1.mutations + sum(x1 != x1old)
    number.x2.mutations <- number.x2.mutations + sum(x2 != x2old)
  }
  list(unname(x1), unname(x2), number.x1.mutations, number.x2.mutations)
}




