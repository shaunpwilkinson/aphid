#' Multiple sequence alignment.
#'
#' \code{align} finds the optimal alignment for a list of sequences involving a
#' progressive multiple sequence alignment, followed by the generation of a profile
#' HMM, an iterative model refinement step, and finally the alignment of the sequences
#' back to the model.
#'
#' @param sequences a list of vectors consisting of symbols emitted from
#' the residue alphabet. The vectors can be of mode "character", or raw bytes in the "DNAbin" or
#' "AAbin" coding scheme (Paradis 2007).
#' @param model an optional profile hidden Markov model (object of class \code{"PHMM"})
#' to align the sequences to.
#' @param refine the method used to iteratively train the model following the
#' initial progressive alignment step.
#' Current supported methods are \code{refine = "Viterbi"} (Viterbi training)
#' and \code{refine = "BaumWelch"}
#' (parameter optimization with the Expectation-Maximization algorithm).
#' @param gapchar the character used to represent gaps in the alignment matrix. Ignored for
#' "DNAbin" and "AAbin" objects.
#' @param residues either NULL (default; emitted residues are automatically
#' detected from the list of sequences), or a case sensitive character vector specifying the
#' residue alphabet (e.g. c(A, C, G, T) for DNA). The character strings "RNA", "DNA", "AA",
#' and "AMINO"are also accepted.
#' Note that the default setting \code{residues = NULL} will not
#' detect rare residues that are not present in the sequence list, and thus will
#' not assign them emission probabilities.
#' Also note that the default option can be slow for large lists of character vectors;
#' therefore specifying the residue alphabet can increase speed in these cases.
#' @param quiet logical argument indicating whether feedback should be printed
#' to the console.
#' @param ... aditional arguments to pass to \code{"Viterbi"} for Viterbi training refinement or
#' \code{"forward"} if using Baum Welch model training.
#' @return a character matrix of aligned sequences.
#' @references Durbin et al, Soding, Paradis 2007.
#' @examples
#' x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#' y <- c("P", "A", "W", "H", "E", "A", "E")
#' z <- align(x, y)
#'
#' @name align
#' @export
#'
align <- function(sequences, model = NULL, seqweights = "Gerstein", refine = "Viterbi", k = 5,
                  maxiter = if(refine == "Viterbi") 10 else 100, maxsize = NULL,
                  inserts = "map", lambda = 0, threshold = 0.5, deltaLL = 1E-07,
                  DI = FALSE, ID = FALSE, residues = NULL, gapchar = "-",
                  pseudocounts = "background", qa = NULL, qe = NULL, quiet = FALSE, ...){
  UseMethod("align")
}

#' @rdname align
#' @export
#'
align.DNAbin <- function(sequences, model = NULL, seqweights = "Gerstein", refine = "Viterbi",
                         k = 5, maxiter = if(refine == "Viterbi") 10 else 100, maxsize = NULL,
                        inserts = "map", lambda = 0, threshold = 0.5, deltaLL = 1E-07,
                        DI = FALSE, ID = FALSE, residues = NULL, gapchar = "-",
                        pseudocounts = "background", qa = NULL, qe = NULL, quiet = FALSE, ...){
  if(is.list(sequences)){
    align.list(sequences, model = model, seqweights = seqweights, refine = refine, k = k,
               maxiter = maxiter, maxsize = maxsize, inserts = inserts, lambda = lambda,
               threshold = threshold, deltaLL = deltaLL, DI = DI, ID = ID,
               residues = residues, gapchar = gapchar, pseudocounts = pseudocounts,
               qa = qa, qe = qe, quiet = quiet, ... = ...)
  }else{
    align.default(sequences, model = model, residues = residues, gapchar = gapchar,
                  pseudocounts = pseudocounts, maxsize = maxsize, quiet = quiet, ... = ...)
  }
}

#' @rdname align
#' @export
#'
align.AAbin <- function(sequences, model = NULL, seqweights = "Gerstein", refine = "Viterbi", k = 5,
                        maxiter = if(refine == "Viterbi") 10 else 100, maxsize = NULL,
                        inserts = "map", lambda = 0, threshold = 0.5, deltaLL = 1E-07,
                        DI = FALSE, ID = FALSE,
                        residues = NULL, gapchar = "-", pseudocounts = "background",
                        qa = NULL, qe = NULL, quiet = FALSE, ...){
  if(is.list(sequences)){
    align.list(sequences, model = model, seqweights = seqweights, refine = refine, k = k,
               maxiter = maxiter, maxsize = maxsize, inserts = inserts, lambda = lambda,
               threshold = threshold, deltaLL = deltaLL, DI = DI, ID = ID,
               residues = residues, gapchar = gapchar, pseudocounts = pseudocounts,
               qa = qa, qe = qe, quiet = quiet, ... = ...)
  }else{
    align.default(sequences, model = model, residues = residues, gapchar = gapchar,
                  pseudocounts = pseudocounts, maxsize = maxsize, quiet = quiet, ... = ...)
  }
}

#' @rdname align
#' @export
#'
align.list <- function(sequences, model = NULL, seqweights = "Gerstein", k = 5,
                       refine = "Viterbi", maxiter = if(refine == "Viterbi") 10 else 100,
                       maxsize = NULL, inserts = "map", lambda = 0, threshold = 0.5, deltaLL = 1E-07,
                       DI = FALSE, ID = FALSE, residues = NULL, gapchar = "-",
                       pseudocounts = "background", qa = NULL, qe = NULL, quiet = FALSE, ...){
  nseq <- length(sequences)
  DNA <- is.DNA(sequences)
  AA <- is.AA(sequences)
  if(DNA) class(sequences) <- "DNAbin" else if(AA) class(sequences) <- "AAbin"
  residues <- alphadetect(sequences, residues = residues, gapchar = gapchar)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  for(i in 1:length(sequences)) sequences[[i]] <- sequences[[i]][sequences[[i]] != gapchar]
  if(is.null(model)){
    if(nseq == 2){
      return(align.default(sequences, model = NULL, residues = residues,  gapchar = gapchar,
                    pseudocounts = pseudocounts, quiet = quiet,... = ...))
    }else if(nseq == 1){
      return(matrix(sequences[[1]], nrow = 1))
    }
    if(!quiet) cat("Calculating pairwise distances\n")
    if(nseq > 100){
      nseeds <- min(nseq, 100 + 2 * ceiling(log(nseq, 2)))
      seeds <- sample(1:length(sequences), size = nseeds)
    }else seeds <- seq_along(sequences)
    if(identical(seqweights, "Gerstein")){
      if(!quiet) cat("Calculating sequence weights\n")
      qds <- kdistance(sequences, k = k, alpha = if(AA) "Dayhoff6" else if(DNA) NULL else residues)
      guidetree <- as.dendrogram(hclust(qds, method = "average"))
      myseqweights <- weight(guidetree, method = "Gerstein")[names(sequences)]
    }else if(is.null(seqweights)){
      myseqweights <- rep(1, nseq)
    }
    phmm <- derive.PHMM.list(sequences, seeds = seeds, refine = refine, maxiter = maxiter,
                        seqweights = myseqweights, k = k, residues = residues, gapchar = gapchar,
                        maxsize = maxsize, inserts = inserts, lambda = lambda,
                        threshold = threshold, deltaLL = deltaLL, DI = DI, ID = ID,
                        pseudocounts = pseudocounts, logspace = TRUE, qa = qa, qe = qe,
                        quiet = quiet, ... = ...)
    if(!quiet) cat("Aligning sequences to model\n")
    # if(!is.null(tmp)) names(sequences) <- tmp
    res <- align.list(sequences, model = phmm, ... = ...)
    if(!quiet) cat("Produced alignment with", ncol(res), "columns (including inserts)\n")
    return(res)
  }else{
    #note changes here also need apply to 'train'
    stopifnot(class(model) == "PHMM")
    # DNA <- is.DNA(sequences)
    # AA <- is.AA(sequences)
    # gapchar <- if(DNA) as.raw(4) else if(AA) as.raw(45) else gapchar
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
    # if(!quiet) cat("Model size:", l, "internal modules\n")
    #nseq <- length(sequences)
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
#' @export
#'
align.default <- function(sequences, model, pseudocounts = "background",
                          residues = NULL, gapchar = "-", maxsize = NULL, quiet = FALSE, ...){
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
    z <- derive.PHMM(sequences, seqweights = "Gerstein", pseudocounts = pseudocounts,
                     residues = residues, logspace = TRUE)
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
    # res.list <- unalign(res)
    # res.weights <- weight(res.list, method = "Gerstein", k = 5)
    # res.phmm <- derive.PHMM.default(res, seqweights = res.weights, quiet = TRUE)
    # res.phmm <- train(res.phmm, res.list, method = "Viterbi", maxiter = 10,
    #                   logspace = TRUE, quiet = TRUE, ... = ...)
    # res <- align.list(sequences = res.list, model = res.phmm, quiet = TRUE, ... = ...)
    return(res)
  }else if(nrow(sequences) > 1 & nrow(model) > 1){
    nx <- nrow(sequences)
    ny <- nrow(model)
    zx <- derive.PHMM(sequences, seqweights = "Gerstein", pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    zy <- derive.PHMM(model,  seqweights = "Gerstein", pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
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
    res.list <- unalign(res)
    res.weights <- weight(res.list, method = "Gerstein", k = 5)
    # cat("\n", maxsize, "\n")
    # cat(sum(apply(res, 2, function(v) !any(v == gapchar))), "\n\n")
    newmaxsize <- if(is.null(maxsize)){
      NULL
    }else{
      max(c(sum(apply(res, 2, function(v) !any(v == gapchar))), maxsize))
    }
    res.phmm <- derive.PHMM.default(res, seqweights = res.weights, maxsize = newmaxsize)
    res.phmm <- train(res.phmm, res.list, method = "Viterbi", maxiter = 3,
                      maxsize = newmaxsize, logspace = TRUE, quiet = TRUE, ... = ...)
    res <- align.list(sequences = res.list, model = res.phmm, ... = ...)
    return(res)
  }else{
    stop("invalid arguments provided for sequences and or y")
  }
}


#' Deconstruct an alignment.
#'
#' \code{unalign} deconstructs an alignment to its component sequences.
#'
#' @param x a matrix consisting of aligned sequences
#' @inheritParams align
#' @return a list of sequences in the same coding scheme as the input alignment (ie DNAbin,
#' AAbin, or plain ASCII characters).
#' @export
#'
unalign <- function(x, gapchar = "-"){
  #x is a matrix representing an alignment
  DNA <- is.DNA(x)
  AA <- is.AA(x)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  if(is.list(x)){
    if(length(x) == 1){
      tmpname <- names(x)
      x <- x[[1]]
      if(is.null(dim(x))){
        x <- matrix(x, nrow = 1)
        rownames(x) <- tmpname
      }
    }
  }
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

