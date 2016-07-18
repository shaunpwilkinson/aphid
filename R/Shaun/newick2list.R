list2newick<-function(x, strip.edges = FALSE){
  fun <- function(taxon, edge=1) paste0("(", taxon, ":", edge, ")")
  tmp <- deparse(x)
  tmp <- paste0(tmp, collapse = "")
  tmp <- gsub(" ", "", tmp)
  tmp <- gsub("structure", "fun", tmp)
  tmp <- gsub("list", "paste", tmp)
  res <- eval(parse(text = tmp))
  res <- gsub(") \\(", ",", res)
  res <- substr(res, start = 2, stop = nchar(res)-3)
  res <- paste0(res, ";")
  if(strip.edges) res <- gsub(":[0-9.]+", "", res)
  res
}

newick2list <- function(x){
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
  tmp
}

list2dendrogram <- function(x) {
  if (is.list(x)) {
    clade_sizes <- sapply(x, function(y) length(unlist(y)))
    attr(x, "members") <- sum(clade_sizes)
    attr(x, "midpoint") <- ((clade_sizes[1] - 1)/2 + 
                              (clade_sizes[1] + (clade_sizes[2] - 1)/2))/2
    if (is.null(attr(x, "height"))) 
      attr(x, "height") <- 0
    if (!(exists("leaf.values"))) 
      leaf.values <- cbind(unlist(x), 1:length(unlist(x)))
    x[] <- lapply(x, function(y) {
      attr(y, "height") <- attr(x, "height") - attr(y,"edge")
      attr(y, "edge") <- NULL
      if (!(is.list(y))) {
        attr(y, "label") <- y[1]
        attr(y, "leaf") <- TRUE
        attr(y, "members") <- 1
        y[1] <- leaf.values[which(leaf.values[, 1] == 
                                    y), 2]
        mode(y) <- "integer"
      }
      y
    })
    x[] <- lapply(x, list2dendrogram)
  }
  x
}# need to manually set class attribute to 'dendrogram' after running function


