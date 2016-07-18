#' Write dendrogram to Newick string
#' 
#' Writes a dendrogram object to a file or connection in Newick 
#' (New Hampshire) format.
#' 
#' @param x an object of class \code{"dendrogram"}.
#' 
#' @param file a character string naming a file or connection to write the output to.
#' If no file path is specified or \code{file = ""} the result will be printed to the console.
#' 
#' @param append a logical value indicating whether the output should be appended to the file.
#' If \code{append = FALSE} the contents of the file will be overwritten (the default setting).
#'  
#' @param strip.edges a logical value indicating whether edge weights should be 
#' removed from the output.
#' 
#' @param dec.places an integer incidating how many decimal places the edge weights should 
#' be rounded to.
#' 
#' @seealso \code{\link{read.dendrogram}} to create a \code{"dendrogram"} object from a 
#' text file.
#' 
#' @examples arrests.hc <- hclust(dist(USArrests[1:6,]), "ave")
#' arrests.den <- as.dendrogram(arrests.hc)
#' write.dendrogram(arrests.den)
#' 
write.dendrogram <- function(x, file = "", append = FALSE, 
                             strip.edges = FALSE, dec.places = 2){
  if(!(inherits(x, 'dendrogram'))) stop("input object must be of class 'dendrogram'")
  renameLeaves <- function(y){
    if(is.leaf(y)){
      y[1] <- attr(y, 'label')
    }
    y
  }
  x <- dendrapply(x, renameLeaves)
  x <- dendrapply(x, unclass)
  addEdges <- function(y){
    if(is.list(y)){
      y[] <- lapply(y, function(z){
        attr(z, 'edge') <- attr(y, 'height') - attr(z, 'height')
        z
      })
      attributes(y)[names(attributes(y)) != 'edge'] <- NULL
      y[] <- lapply(y, addEdges)
    }else{
      attributes(y)[names(attributes(y)) != 'edge'] <- NULL
    }
    y
  }
  x <- addEdges(x)
  attr(x, 'edge') <- 0
  tmp <- deparse(x)
  tmp <- paste0(tmp, collapse = "")
  tmp <- gsub(" ", "", tmp)
  fun <- function(taxon, edge=1) paste0("(", taxon, ":", edge, ")")
  tmp <- gsub("structure", "fun", tmp)
  tmp <- gsub("list", "paste", tmp)
  res <- eval(parse(text = tmp))
  res <- gsub(") \\(", ",", res)
  res <- substr(res, start = 2, stop = nchar(res)-3)
  res <- paste0(res, ";")
  if(strip.edges){
    res <- gsub(":[0-9.]+", "", res)
  }else{
    target <- paste0(c("(\\.", rep("[0-9]", dec.places), ")[0-9]+"), collapse = "")
    res <- gsub(target, "\\1", res)
    if(dec.places == 0) res <- gsub("\\.", "", res)
  }
  if(file == ""){
    return(res)
  }else{
    cat(res, file = file, append = append, sep = "\n")
  }
}

