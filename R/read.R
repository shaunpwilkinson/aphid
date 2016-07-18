#' Read dendrogram from parenthetic text
#' 
#' Reads a file or character string in Newick (New Hampshire) format into 
#' an object of class \code{dendrogram}.
#' 
#' @param file the name of the file to read the data from.
#' @param text character string: if a text argument is provided instead of a 
#' file path then data are read from the value of text via a text connection.
#' @param strip.edges a logical value indicating whether edge weights 
#' provided in the Newick string should be ignored.
#' @param ... further arguments to be passed to \code{scan}.
#'  
#' @return an object of class \code{"dendrogram"}.
#' 
#' @seealso \code{\link{write.dendrogram}} to write an object of 
#' class \code{"dendrogram"} to a text string.
#' 
#' @examples mynewick <- "(A:1,(B:0.5,C:0.2):0.6);"
#' mydendrogram <- read.dendrogram(text = mynewick)
#' plot(mydendrogram) 
#' 
read.dendrogram <- function(file = "", text = NULL, strip.edges = FALSE, ...){
  if(!is.null(text)){
    if(!is.character(text)) 
      stop("argument 'text' must be of mode character.")
    x <- text
  }else{
    x <- scan(file = file, what = "", sep = "\n", quiet = TRUE, ...)
  }
  if (identical(x, character(0))) {
    warning("empty character string.")
    return(NULL)
  }
  x <- paste0(x, collapse = "")
  if(strip.edges) x <- gsub(":([0-9.]+)", "", x)
  newick.has.edgeweights <- grepl(":", x)
  if(newick.has.edgeweights){ 
    tmp <- gsub("([[:alnum:]]+):", "\\(\"\\1\"):", x)    
    tmp <- gsub(";", ":1", tmp)
  }else{
    tmp <- gsub("([[:alnum:]]+)", "\\(\"\\1\")", x)
    tmp <- gsub("(\\))","\\1:1", tmp)
    tmp <- gsub(";", "", tmp)
  }
  tmp <- gsub("\\(", "structure\\(\\(", tmp)
  tmp <- gsub("ture\\(\\(struc", "ture\\(list\\(struc", tmp)
  tmp <- gsub(":([0-9.]+)", ",edge = \\1)", tmp)
  tmp <- eval(parse(text = tmp))
  attr(tmp, 'edge') <- 0
  setAttributes <- function(x){ # x is a list with 'edge' attributes
    if(is.list(x)){
      clade_sizes <- sapply(x, function(y) length(unlist(y)))
      attr(x, 'members') <- sum(clade_sizes)
      attr(x, 'midpoint') <- ((clade_sizes[1] - 1)/2 + 
                                (clade_sizes[1] + (clade_sizes[2] - 1)/2))/2
      if(is.null(attr(x, 'height'))) attr(x, 'height') <- 0
      if(!(exists("leaf.values"))) leaf.values <- cbind(unlist(x), 1:length(unlist(x)))
      x[] <- lapply(x, function(y){
        attr(y, 'height') <- attr(x, 'height') - attr(y, 'edge')
        attr(y, 'edge') <- NULL
        if(!(is.list(y))){
          attr(y, "label") <- y[1]
          attr(y, "leaf") <- TRUE
          attr(y, 'members') <- 1
          y[1] <- leaf.values[which(leaf.values[,1] == y), 2]
          mode(y) <- 'integer'
        }
        y
      })
      x[] <- lapply(x, setAttributes)
    }
    x
  }
  res <- setAttributes(tmp)
  attr(res, 'edge') <- NULL
  attr(res, 'class') <- 'dendrogram'
  min.height <- min(unlist(dendrapply(res, attr, 'height')))
  res <- dendrapply(res, function(y){
    attr(y, 'height') <- attr(y, 'height') - min.height
    y
  })
  if(!(newick.has.edgeweights)){
    res <- dendrapply(res, function(y){
      if(is.leaf(y)){
        attr(y, 'height') <- 0
      }
      y
    })
  }  
  res
}
