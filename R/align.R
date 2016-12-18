#' Multiple sequence alignment.
#'
#' \code{align} finds the optimal alignment for a list of sequences using a hybrid
#' algorithm involving a progressive multiple sequence alignment, the generation of a profile
#' HMM, an iterative model refinement step, and finally the alignment of the sequences
#' back to the model.
#'
#' @param sequences a list of character vectors consisting of symbols from
#' the residue alphabet.
#' @param model an optional profile hidden Markov model (object of class \code{"PHMM"})
#' to align the sequences to.
#' @param refine the method used to iteratively train the model.
#' @param gapchar the character used to represent gaps in the alignment matrix.
#' @param residues either NULL (default; emitted residues are automatically
#' detected from the list of sequences), or a case sensitive character vector specifying the
#' residue alphabet (e.g. c(A, C, G, T) for DNA). The character strings "RNA", "DNA", "AA",
#' and "AMINO"are also accepted.
#' Note that the default option can be slow for large lists of character vectors;
#' therefore specifying the residue alphabet can increase speed in these cases.
#' Also note that the default setting \code{residues = NULL} will not
#' detect rare residues that are not present in the sequence list, and thus will
#' not assign them emission probabilities.
#' @param refine the method used to refine the model following progressive alignment.
#' Current supported methods are Viterbi training (\code{refine = "Viterbi"}) and
#' Baum Welch parameter optimization (\code{refine = "BaumWelch"}).
#' @param DI logical indicating whether delete -> insert transitions should be allowed
#' in the profile HMM. Defaults to FALSE, and not recommended for small training sets when
#' \code{refine = "BaumWelch"} due
#' to the tendency to converge to suboptimal local optima.
#' @param ID logical indicating whether insert -> delete transitions should be allowed
#' in the profile HMM. Defaults to FALSE, and similar to DI, not recommended for small
#' training sets when \code{refine = "BaumWelch"} due
#' to the tendency to converge to suboptimal local optima.
#' @param quiet logical argument indicating whether feedback should be printed
#' to the console.
#' @param ... aditional arguments to pass to \code{"Viterbi"} (for Viterbi training refinement) or
#' \code{"forward"} (for Baum Welch model training).
#' @return a character matrix of aligned sequences.
#' @references Soding...
#' @examples
#' x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#' y <- c("P", "A", "W", "H", "E", "A", "E")
#' z <- align(x, y)
#'
#'
#' @name Viterbi
align <- function(sequences, model = NULL, refine = "Viterbi", inserts = "map", lambda = 0,
                  threshold = 0.5, quiet = FALSE, residues = NULL, gapchar = "-", k = 5,
                  pseudocounts = "background", ...){
  UseMethod("align")
}

#' @rdname align
align.DNAbin <- function(sequences, model = NULL, refine = "Viterbi", inserts = "map",
                         lambda = 0, threshold = 0.5, quiet = FALSE, residues = NULL,
                         gapchar = "-", k = 5, pseudocounts = "background", ...){
  if(is.list(sequences)){
    align.list(sequences, model = model, refine = refine, inserts = inserts,
               lambda = lambda, threshold = threshold, quiet = quiet, residues = residues,
               gapchar = gapchar, k = k, pseudocounts = pseudocounts, ... = ...)
  }else{
    align.default(sequences, model = model, residues = residues, quiet = quiet,
                  gapchar = gapchar, pseudocounts = pseudocounts, ... = ...)
  }
}

#' @rdname align
align.AAbin <- function(sequences, model = NULL, refine = "Viterbi", inserts = "map",
                         lambda = 0, threshold = 0.5, quiet = FALSE, residues = NULL,
                         gapchar = "-", k = 5, pseudocounts = "background", ...){
  if(is.list(sequences)){
    align.list(sequences, model = model, refine = refine, inserts = inserts,
               lambda = lambda, threshold = threshold, quiet = quiet, residues = residues,
               gapchar = gapchar, k = k, pseudocounts = pseudocounts, ... = ...)
  }else{
    align.default(sequences, model = model, residues = residues, quiet = quiet,
                  gapchar = gapchar, pseudocounts = pseudocounts, ... = ...)
  }
}

#' @rdname align
align.list <- function(sequences, model = NULL, refine = "Viterbi", inserts = "map",
                       lambda = 0, threshold = 0.5, quiet = FALSE, residues = NULL,
                       gapchar = "-", k = 5, pseudocounts = "background", ...){
  nsq <- length(sequences)
  if(nsq == 2 & is.null(model)){
    align.default(sequences, model = NULL, residues = residues, quiet = quiet,
                  gapchar = gapchar, pseudocounts = pseudocounts, ... = ...)
  }
  DNA <- is.DNA(sequences)
  AA <- is.AA(sequences)
  residues <- alphadetect(sequences, residues = residues, gapchar = gapchar)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  for(i in 1:length(sequences)) sequences[[i]] <- sequences[[i]][sequences[[i]] != gapchar]
  if(is.null(model)){
    if(!quiet) cat("Calculating pairwise distances\n")
    tmp <- names(sequences)
    names(sequences) <- paste0("S", 1:nsq)
    if(DNA) class(sequences) <- "DNAbin" else if(AA) class(sequences) <- "AAbin"
    qds <- kdistance(sequences, k = k, alpha = if(AA) "Dayhoff6" else if(DNA) NULL else residues)
    if(!quiet) cat("Building guide tree\n")
    guidetree <- as.dendrogram(hclust(qds, method = "average"))
    if(!quiet) cat("Calculating sequence weights\n")
    seqweights <- weight(guidetree, method = "Gerstein")
    ## add embed step here
    newick <- write.dendrogram(guidetree, strip.edges = TRUE)
    newick <- gsub(";", "", newick)
    newick <- gsub("\\(", "align\\(", newick)
    newick <- gsub("\\)", ", ... = ...\\)", newick)
    #if(type == "global") newick <- gsub("\\)", ", type='global'\\)", newick)
    if(!quiet) cat("Building initial alignment\n")

    ### gets messy here
    msa1 <- with(sequences, eval(parse(text = newick)))
    if(!quiet) cat("Deriving profile hidden Markov model\n")
    omniphmm <- derive.PHMM(msa1, seqweights = seqweights, pseudocounts = pseudocounts,
                            inserts = inserts, lambda = lambda, threshold = threshold)
    if(refine %in% c("Viterbi", "BaumWelch")){
      if(!quiet) cat("Refining model\n")
      finalphmm <- train(omniphmm, sequences, method = refine, quiet = quiet, ... = ...)
    }else if (refine == "none"){
      finalphmm <- omniphmm
    }else stop("Argument 'refine' must be set to either 'Viterbi', 'BaumWelch' or 'none' (case sensitive).")
    if(!quiet) cat("Aligning sequences to model\n")
    if(!is.null(tmp)) names(sequences) <- tmp
    res <- align(sequences, model = finalphmm)
    if(!quiet) cat("Done\n")
    return(res)
  }else{
    #note changes here also need apply to 'train'
    stopifnot(class(model) == "PHMM")
    DNA <- is.DNA(sequences)
    AA <- is.AA(sequences)
    gapchar <- if(DNA) as.raw(4) else if(AA) as.raw(45) else gapchar
    if(!is.list(sequences)){
      if(is.null(dim(sequences))){
        seqname <- deparse(substitute(sequences))
        sequences <- list(sequences)
        names(sequences) <- seqname
      }else{
        sequences <- unalign(sequences)
      }
    }
    l <- model$size
    nseq <- length(sequences)
    out <- matrix(nrow = nseq, ncol = 2 * l + 1)
    rownames(out) <- attr(sequences, "names")
    colnames(out) <- rep("I", 2 * l + 1)
    if(l > 0) colnames(out)[seq(2, 2 * l, by = 2)] <- 1:l
    score <- 0
    for(i in 1:nseq){
      alignment <- Viterbi(model, sequences[i], ... = ...)
      score <- score + alignment$score
      newrow <- rep(gapchar, length(alignment$path))
      if(DNA | AA){
        newrow <- rep("-", length(alignment$path))
        newrow[alignment$path > 0] <- if(DNA){
          ape::as.character.DNAbin(sequences[[i]])
        }else{
          ape::as.character.AAbin(sequences[[i]])
        }
      }else{
        newrow <- rep(gapchar, length(alignment$path))
        newrow[alignment$path > 0] <- sequences[[i]]
      }
      inserts <- alignment$path == 2
      if(any(inserts)){
        itp <- apply(rbind(c(F, inserts), c(inserts, F)), 2, decimal, 2)
        ist <- which(itp == 1)
        ien <- which(itp == 2) - 1
        runs <- mapply(":", ist, ien, SIMPLIFY = FALSE)
        newrow[ist] <- sapply(runs, function(e) paste0(newrow[e], collapse = ""))
        itp <- itp[seq_along(newrow)] # removes final transition which can't be a 1(MI) or a 3(II)
        insflag <- itp == 1
        newrow <- newrow[itp < 3] # now newrow is x$size + number of MI transitions
        insflag <- insflag[itp < 3] #now same length as newrow
        insflag2 <- c(FALSE, insflag, FALSE, FALSE)
        tuples <- rbind(insflag2[-(length(insflag2))], insflag2[-1])
        decs <- apply(tuples, 2, decimal, 2)
        tuples <- tuples[, decs != 2, drop = FALSE] #remove true->falses
        tuples <- tuples[,-ncol(tuples), drop = FALSE]
        tuples[1, -1] <- TRUE
      }else{
        tuples <- rbind(c(FALSE, rep(TRUE, l)), rep(FALSE, l + 1))
      }
      indices <- as.vector(tuples)[-1]
      newrow2 <- rep(if(DNA | AA) "-" else gapchar, 2 * l + 1)
      names(newrow2) <- c(rep(c("I", "M"), l), "I") #seq(0.5, x$size + 0.5, by = 0.5)
      newrow2[indices] <- newrow
      out[i,] <- newrow2
    }
    discardcols <- apply(out, 2, function(v) all(v == if(DNA | AA) "-" else gapchar))
    matchcols <- c(FALSE, rep(c(TRUE, FALSE), l))
    out <- out[, matchcols | !discardcols, drop = FALSE]
    out <- as.list(as.data.frame(out, stringsAsFactors = F))
    fun <- function(e){
      ee <- strsplit(e, split = "")
      elengths <- sapply(ee, length)
      if(all(elengths == 1)) return(unlist(ee))
      maxlen <- max(elengths)
      no.gapstoappend <- maxlen - elengths
      appges <- lapply(no.gapstoappend, function(v) rep(if(DNA | AA) "-" else gapchar, v))
      return(t(mapply(c, ee, appges)))
    }
    out[names(out) == "I"] <- lapply(out[names(out) == "I"], fun)
    ncols <- sapply(out, function(e) if(is.null(dim(e))) 1 else ncol(e))
    totcols <- sum(ncols)
    res <- matrix(nrow = nseq, ncol = totcols)
    rownames(res) <- names(sequences)
    rescolnames <- vector(mode = "character", length = totcols)
    counter <- 1
    for(i in seq_along(ncols)){
      if(ncols[i] == 1){
        res[, counter] <- out[[i]]
        rescolnames[counter] <- names(out)[i]
      }else{
        endrange <- counter + ncols[i] - 1
        res[, counter:endrange] <- out[[i]]
        rescolnames[counter:endrange] <- rep("I", ncols[i])
      }
      counter <- counter + ncols[i]
    }
    colnames(res) <- rescolnames
    inserts <- sapply(colnames(res), function(v) grepl("I", v))
    #colnames(res)[inserts] <- "I"
    #rownames(res) <- names(sequences)
    if(DNA) res <- ape::as.DNAbin(res) else if(AA) res <- ape::as.AAbin(res)
    attr(res, "score") <- score
    attr(res, "inserts") <- unname(inserts)
    return(res)
  }
}


#' @rdname align
align.default <- function(sequences, model, pseudocounts = "background",
                          residues = NULL, gapchar = "-", quiet = FALSE, ...){
  if(is.null(model)) return(sequences)
  DNA <- is.DNA(sequences)
  AA <- is.AA(sequences)
  if(DNA){
    if(!is.DNA(model)) stop("class(sequences) and class(model) must match")
    gapchar <- as.raw(4)
    # changes here need also apply in Viterbi.default and Viterbi.PHMM
    if(is.list(sequences)){
      if(length(sequences) == 1){
        if(!is.matrix(sequences)) sequences <- matrix(sequences[[1]], nrow = 1, dimnames = list(names(sequences), NULL))
        class(sequences) <- "DNAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(sequences)) {
        sequences <- matrix(sequences, nrow = 1, dimnames = list(deparse(substitute(sequences)), NULL))
      }
      class(sequences) <- "DNAbin"
    }
    if(is.list(model)){
      if(length(model) == 1){
        if(!is.matrix(model)) model <- matrix(model[[1]], nrow = 1, dimnames = list(names(model), NULL))
        class(model) <- "DNAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(model)) {
        model <- matrix(model, nrow = 1, dimnames = list(deparse(substitute(model)), NULL))
      }
      class(model) <- "DNAbin"
    }
  }else if(AA){
    if(!is.AA(model)) stop("class(sequences) and class(model) must match")
    gapchar <- as.raw(45)
    # changes here need also apply in Viterbi.default and Viterbi.PHMM
    if(is.list(sequences)){
      if(length(sequences) == 1){
        if(!is.matrix(sequences)) sequences <- matrix(sequences[[1]], nrow = 1, dimnames = list(names(sequences), NULL))
        class(sequences) <- "AAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(sequences)) sequences <- matrix(sequences, nrow = 1, dimnames = list(deparse(substitute(sequences)), NULL))
      class(sequences) <- "AAbin"
    }
    if(is.list(model)){
      if(length(model) == 1){
        if(!is.matrix(model)) model <- matrix(model[[1]], nrow = 1, dimnames = list(names(model), NULL))
        class(model) <- "AAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(model)) model <- matrix(model, nrow = 1, dimnames = list(deparse(substitute(model)), NULL))
      class(model) <- "AAbin"
    }
  }else{
    if(!is.matrix(sequences)){
      sequences <- matrix(sequences, nrow = 1, dimnames = list(deparse(substitute(sequences)), NULL))
      #rownames(sequences) <-
    }
    if(!is.matrix(model)){
      model <- matrix(model, nrow = 1, dimnames = list(deparse(substitute(model)), NULL))
    }
  }
  if(nrow(sequences) == 1 & nrow(model) == 1){
    alig <- Viterbi(sequences, model, ... = ...)
    #, offset = offset) ###not necessary for vec vs vec
    xind <- yind <- alig$path
    xind[alig$path != 2] <- 1:length(sequences)
    xind[alig$path == 2] <- 0
    newx <- c(gapchar, as.vector(sequences))[xind + 1]
    yind[alig$path != 0] <- 1:length(model)
    yind[alig$path == 0] <- 0
    newy <- c(gapchar, as.vector(model))[yind + 1]
    res <- rbind(newx, newy)
    rownames(res) <- c(rownames(sequences), rownames(model))
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    return(res)
  }else if(sum(c(nrow(sequences) == 1, nrow(model) == 1)) == 1){
    if(nrow(sequences) == 1){
      tmp <- sequences
      sequences <- model
      model <- tmp
      rm(tmp) # the old switcharoo
    }
    n <- nrow(sequences)
    z <- derive.PHMM(sequences, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    l <- z$size
    alignment <- Viterbi(z, model, ... = ...) ### changed from alig
    path <- alignment$path
    # lay x out as list with insert elements
    newrow <- vector(length = l * 2 + 1, mode = "list")
    odds <- seq(from = 1, to = length(newrow), by = 2) #insert columns
    evens <- seq(from = 2, to = length(newrow), by = 2) # match columns
    newrow[evens] <- lapply(which(!z$inserts), function(e) e)
    if(any(z$inserts)){
      itp <- apply(rbind(c(FALSE, z$inserts), c(z$inserts, FALSE)), 2, decimal, from = 2)
      ist <- which(itp == 1)
      ien <- which(itp == 2) - 1
      newrow[odds][which(itp[itp < 2] == 1)] <- mapply(":", ist, ien, SIMPLIFY = FALSE)
    }
    newrow <- lapply(newrow, function(e) if(is.null(e)) 0 else e)
    newx <- lapply(newrow, function(e) sequences[, e, drop = FALSE])
    # analogous list for y but witout insert elements
    newrow <- lapply(1:ncol(model), function(e) e)
    newy <- lapply(newrow, function(e) model[, e, drop = FALSE])
    #
    newxrows <- matrix(gapchar, nrow = n, ncol = ncol(sequences) + ncol(model))
    rownames(newxrows) <- rownames(sequences)
    newyrow <- matrix(gapchar, nrow = 1, ncol = ncol(sequences) + ncol(model))
    rownames(newyrow) = rownames(model)
    is.insert <- vector(mode = "logical", length = ncol(sequences) + ncol(model))
    if(DNA){
      class(newxrows) <- "DNAbin"
      class(newyrow) <- "DNAbin"
    }else if(AA){
      class(newxrows) <- "AAbin"
      class(newyrow) <- "AAbin"
    }
    position <- 1
    xcounter <- 2 * alignment$start[1]
    ycounter <- alignment$start[2]
    if(xcounter == 2 & ycounter == 1){ # global and semiglobal alignments
      rightshift <- ncol(newx[[1]])
      if(rightshift > 0){
        newxrows <- insert(newx[[1]], into = newxrows, at = 1)
        position <- position + rightshift
        is.insert[1:rightshift] <- TRUE
      }
    }
    for(i in seq_along(path)){
      if(path[i] == 1){ #Match state
        rightshift <- 1 + ncol(newx[[xcounter + 1]])
        # match + insert
        newxrows <- insert(newx[[xcounter]], into = newxrows, at = position)
        newyrow <- insert(newy[[ycounter]], into = newyrow, at = position)
        position <- position + 1
        if(rightshift > 1){
          newxrows <- insert(newx[[xcounter + 1]], into = newxrows, at = position)
          is.insert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
        ycounter <- ycounter + 1
      }else if(path[i] == 0){ #Delete state
        rightshift <- 1 + ncol(newx[[xcounter + 1]])
        newxrows <- insert(newx[[xcounter]], into = newxrows, at = position)
        position <- position + 1
        if(rightshift > 0){
          newxrows <- insert(newx[[xcounter + 1]], into = newxrows, at = position)
          is.insert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
      }else if(path[i] == 2){
        rightshift <- 1
        newyrow <- insert(newy[[ycounter]], into = newyrow, at = position)
        position <- position + 1
        ycounter <- ycounter + 1
      }
    }
    position <- position - 1
    newxrows <- newxrows[, 1:position]
    newyrow <- newyrow[, 1:position]
    is.insert <- is.insert[1:position]
    res <- rbind(newxrows, newyrow)
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    return(res)
  }else if(nrow(sequences) > 1 & nrow(model) > 1){
    nx <- nrow(sequences)
    ny <- nrow(model)
    zx <- derive.PHMM(sequences, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    zy <- derive.PHMM(model, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    lx <- zx$size
    ly <- zy$size
    alignment <- Viterbi(zx, zy, ... = ...)
    path <- alignment$path
    # lay x out as list with insert elements
    newrow <- vector(length = lx * 2 + 1, mode = "list")
    odds <- seq(from = 1, to = length(newrow), by = 2) #insert columns
    evens <- seq(from = 2, to = length(newrow), by = 2) # match columns
    newrow[evens] <- lapply(which(!zx$inserts), function(e) e)
    if(any(zx$inserts)){
      itp <- apply(rbind(c(FALSE, zx$inserts), c(zx$inserts, FALSE)), 2, decimal, from = 2)
      ist <- which(itp == 1)
      ien <- which(itp == 2) - 1
      newrow[odds][which(itp[itp < 2] == 1)] <- mapply(":", ist, ien, SIMPLIFY = FALSE)
    }
    newrow <- lapply(newrow, function(e) if(is.null(e)) 0 else e)
    newx <- lapply(newrow, function(e) sequences[, e, drop = FALSE])
    # lay y out as list with insert elements
    newrow <- vector(length = ly * 2 + 1, mode = "list")
    odds <- seq(from = 1, to = length(newrow), by = 2) #insert columns
    evens <- seq(from = 2, to = length(newrow), by = 2) # match columns
    newrow[evens] <- lapply(which(!zy$inserts), function(e) e)
    if(any(zy$inserts)){
      itp <- apply(rbind(c(FALSE, zy$inserts), c(zy$inserts, FALSE)), 2, decimal, from = 2)
      ist <- which(itp == 1)
      ien <- which(itp == 2) - 1
      newrow[odds][which(itp[itp < 2] == 1)] <- mapply(":", ist, ien, SIMPLIFY = FALSE)
    }
    newrow <- lapply(newrow, function(e) if(is.null(e)) 0 else e)
    newy <- lapply(newrow, function(e) model[, e, drop = FALSE])
    #create output alignment
    newxrows <- matrix(gapchar, nrow = nx, ncol = ncol(sequences) + ncol(model))
    rownames(newxrows) <- rownames(sequences)
    newyrows <- matrix(gapchar, nrow = ny, ncol = ncol(sequences) + ncol(model))
    rownames(newyrows) <- rownames(model)
    is.insert <- vector(mode = "logical", length = ncol(sequences) + ncol(model))
    if(DNA){
      class(newxrows) <- "DNAbin"
      class(newyrows) <- "DNAbin"
    }else if(AA){
      class(newxrows) <- "AAbin"
      class(newyrows) <- "AAbin"
    }
    position <- 1
    xcounter <- 2 * alignment$start[1]
    ycounter <- 2 * alignment$start[2]
    if(xcounter == 2 & ycounter == 2){
      rightshift <- max(c(ncol(newx[[1]]), newy[[1]]))
      if(rightshift > 0){
        newxrows <- insert(newx[[1]], into = newxrows, at = 1)
        newyrows <- insert(newy[[1]], into = newyrows, at = 1)
        position <- position + rightshift
        is.insert[1:rightshift] <- TRUE
      }
    }
    for(i in seq_along(path)){
      if(path[i] == 2){ #MM
        rightshift <- 1 + max(c(ncol(newx[[xcounter + 1]]), ncol(newy[[ycounter + 1]])))
        # match + insert
        newxrows <- insert(newx[[xcounter]], into = newxrows, at = position)
        newyrows <- insert(newy[[ycounter]], into = newyrows, at = position)
        position <- position + 1
        if(rightshift > 1){
          newxrows <- insert(newx[[xcounter + 1]], into = newxrows, at = position)
          newyrows <- insert(newy[[ycounter + 1]], into = newyrows, at = position)
          is.insert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
        ycounter <- ycounter + 2
      }else if(path[i] < 2){
        rightshift <- 1 + ncol(newx[[xcounter + 1]])
        newxrows <- insert(newx[[xcounter]], into = newxrows, at = position)
        position <- position + 1
        if(rightshift > 0){
          newxrows <- insert(newx[[xcounter + 1]], into = newxrows, at = position)
          is.insert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        xcounter <- xcounter + 2
      }else if(path[i] > 2){
        rightshift <- 1 + ncol(newy[[ycounter + 1]])
        newyrows <- insert(newy[[ycounter]], into = newyrows, at = position)
        position <- position + 1
        if(rightshift > 0){
          newyrows <- insert(newy[[ycounter + 1]], into = newyrows, at = position)
          is.insert[position:(position + rightshift - 1)] <- TRUE
          position <- position + rightshift - 1
        }
        ycounter <- ycounter + 2
      }
    }
    position <- position - 1
    newxrows <- newxrows[, 1:position]
    newyrows <- newyrows[, 1:position]
    is.insert <- is.insert[1:position]
    res <- rbind(newxrows, newyrows)
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    return(res)
  }else{
    stop("invalid arguments provided for sequences and or y")
  }
}


#' Deconstruct an alignment to its component sequences.
#'
#'
unalign <- function(x, gapchar = "-"){
  #x is a matrix representing an alignment
  DNA <- is.DNA(x)
  AA <- is.AA(x)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  res <- vector(mode = "list", length = nrow(x))
  for(i in 1:nrow(x)) res[[i]] <- x[i, x[i, ] != gapchar, drop = TRUE]
  if(AA){
    res <- lapply(res, unclass)
    class(res) <- "AAbin"
  }else if(DNA){
    res <- lapply(res, unclass)
    class(res) <- "DNAbin"
  }
  names(res) <- rownames(x)
  return(res)
}

