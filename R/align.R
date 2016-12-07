#' Multiple sequence alignment.
#'
#' \code{align} finds the optimal alignment for a list of sequences using a hybrid
#' algorithm involving a progressive multiple sequence alignment, the generation of a profile
#' HMM, an iterative model refinement step, and finally the alignment of the sequences
#' back to the model.
#'
#' @param sequences a list of character vectors consisting of symbols from
#' the residue alphabet.
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
#' @param ... aditional arguments to pass to \code{"train"}.
#'
#'
align <- function(sequences, type = "semiglobal", residues = NULL,
                  gapchar = "-", k = 5, refine = "Viterbi", DI = FALSE, ID = FALSE,
                  quiet = FALSE, cpp = TRUE, ...){
  if(!(is.list(sequences))) stop("invalid 'sequences' argument")
  nsq <- length(sequences)
  #if(is.null(attr(sequences, "names")))
  DNA <- is.DNA(sequences)
  AA <- is.AA(sequences)
  residues <- alphadetect(sequences, residues = residues, gapchar = gapchar)
  gapchar <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gapchar
  for(i in 1:length(sequences)) sequences[[i]] <- sequences[[i]][sequences[[i]] != gapchar]
  if(!quiet) cat("calculating pairwise distances\n")
  tmp <- names(sequences)
  names(sequences) <- paste0("S", 1:nsq)
  qds <- kdistance(sequences, k = k, alpha = if(AA) "Dayhoff6" else if(DNA) NULL else residues)
  if(!quiet) cat("building guide tree\n")
  guidetree <- as.dendrogram(hclust(qds, method = "average"))
  if(!quiet) cat("calculating sequence weights\n")
  seqweights <- weight(guidetree, method = "Gerstein")
  newick <- write.dendrogram(guidetree, strip.edges = TRUE)
  newick <- gsub(";", "", newick)
  newick <- gsub("\\(", "alignpair\\(", newick)
  if(type == "global") newick <- gsub("\\)", ", type='global'\\)", newick)
  if(!quiet) cat("building initial alignment\n")
  msa1 <- with(sequences, eval(parse(text = newick)))
  if(!quiet) cat("deriving profile hidden Markov model\n")
  omniphmm <- derive.PHMM(msa1, seqweights = seqweights, DI = DI, ID = ID,
                          pseudocounts = "background")
  if(refine == "Viterbi"){
    if(!quiet) cat("refining model\n")
    finalphmm <- train(omniphmm, sequences, method = refine,
                       DI = DI, ID = ID,
                       quiet = quiet, cpp = cpp,
                       ... = ...)
  }else if(refine == "BaumWelch"){
    if(!quiet) cat("refining model\n")
    finalphmm <- train(omniphmm, sequences, method = refine,
                       DI = DI, ID = ID,
                       quiet = quiet, cpp = cpp, ... = ...) ### condense
  }else if (refine == "none"){
    finalphmm <- omniphmm
  }else stop("argument 'refine' must be set to either 'Viterbi', 'BaumWelch' or 'none'")
  if(!quiet) cat("aligning sequences to model\n")
  if(!is.null(tmp)) names(sequences) <- tmp
  res <- align2phmm(sequences, model = finalphmm)
  if(!quiet) cat("done\n")
  return(res)
}

#' Pairwise alignment of sequences and/or multiple sequence alignments.
#'
#' \code{alignpair} uses the Viterbi algorithm to find the optimal alignment
#' between two sequences, a sequence and an alignment, or two alignments.
#'
#' @param x,y character vectors or character matrices of aligned sequences.
#' @param d,e gap opening and gap extension penalties for pairwise sequence
#' alignment.
#' @param S an optional substitution matrix with \code{dimnames} attributes
#' corresponding to the residue alphabet. If NULL matches are
#' scored as 1 and mismatches as -1.
#' @param type a character string specifying whether the alignment should be
#' 'global' (penalized end gaps), 'semiglobal' (default; free end gaps) or
#' local (highest scoring subalignment).
#' @param residues either NULL (default; emitted residues are automatically
#' detected from the input sequences), or a case sensitive character vector specifying the
#' residue alphabet (e.g. c(A, C, G, T) for DNA).
#' Note that the former option can be slow for large character vectors;
#' therefore specifying the residue alphabet can increase speed in these cases.
#' Also note that the default setting \code{residues = NULL} will not
#' detect rare residues that are not present in the input sequences, and thus will
#' not assign them emission probabilities.
#' @return a character matrix of aligned sequences.
#' @references Soding...
#' @examples
#' x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#' y <- c("P", "A", "W", "H", "E", "A", "E")
#' z <- align(x, y)
#' alignpair(x, z)
#'
alignpair <- function(x, y, d = 8, e = 2, S = NULL, qe = NULL,
                      type = "semiglobal", offset = 0, windowspace = "WilburLipman",
                      pseudocounts = "background",
                      residues = NULL, gapchar = "-", cpp = TRUE){
  DNA <- is.DNA(x)
  AA <- is.AA(x)
  if(DNA){
    if(!is.DNA(y)) stop("class(x) and class(y) must match")
    gapchar <- as.raw(4)
    # changes here need also apply in Viterbi.default and Viterbi.PHMM
    if(is.list(x)){
      if(length(x) == 1){
        if(!is.matrix(x)) x <- matrix(x[[1]], nrow = 1, dimnames = list(names(x), NULL))
        class(x) <- "DNAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(x)) x <- matrix(x, nrow = 1, dimnames = list(deparse(substitute(x)), NULL))
      class(x) <- "DNAbin"
    }
    if(is.list(y)){
      if(length(y) == 1){
        if(!is.matrix(y)) y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "DNAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(y)) y <- matrix(y, nrow = 1, dimnames = list(deparse(substitute(y)), NULL))
      class(y) <- "DNAbin"
    }
  }else if(AA){
    if(!is.AA(y)) stop("class(x) and class(y) must match")
    gapchar <- as.raw(45)
    # changes here need also apply in Viterbi.default and Viterbi.PHMM
    if(is.list(x)){
      if(length(x) == 1){
        if(!is.matrix(x)) x <- matrix(x[[1]], nrow = 1, dimnames = list(names(x), NULL))
        class(x) <- "AAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(x)) x <- matrix(x, nrow = 1, dimnames = list(deparse(substitute(x)), NULL))
      class(x) <- "AAbin"
    }
    if(is.list(y)){
      if(length(y) == 1){
        if(!is.matrix(y)) y <- matrix(y[[1]], nrow = 1, dimnames = list(names(y), NULL))
        class(y) <- "AAbin"
      }else stop("Invalid input: multi-sequence list")
    }else{
      if(!is.matrix(y)) y <- matrix(y, nrow = 1, dimnames = list(deparse(substitute(y)), NULL))
      class(y) <- "AAbin"
    }
  }else{
    if(!is.matrix(x)){
      x <- matrix(x, nrow = 1, dimnames = list(deparse(substitute(x)), NULL))
    }
    if(!is.matrix(y)){
      y <- matrix(y, nrow = 1, dimnames = list(deparse(substitute(y)), NULL))
    }
  }
  if(nrow(x) == 1 & nrow(y) == 1){
    alig <- Viterbi(x, y, d = d, e = e, S = S, type = type,
                    windowspace = windowspace, cpp = cpp)
    #, offset = offset) ###not necessary for vec vs vec
    xind <- yind <- alig$path
    xind[alig$path != 2] <- 1:length(x)
    xind[alig$path == 2] <- 0
    newx <- c(gapchar, as.vector(x))[xind + 1]
    yind[alig$path != 0] <- 1:length(y)
    yind[alig$path == 0] <- 0
    newy <- c(gapchar, as.vector(y))[yind + 1]
    res <- rbind(newx, newy)
    rownames(res) <- c(rownames(x), rownames(y))
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    return(res)
  }else if(sum(c(nrow(x) == 1, nrow(y) == 1)) == 1){
    if(nrow(x) == 1){
      tmp <- x
      x <- y
      y <- tmp
      rm(tmp) # the old switcharoo
    }
    #residues <- alphadetect(x, residues = residues, gapchar = gapchar)
    n <- nrow(x)
    z <- derive.PHMM(x, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    l <- z$size
    alig <- Viterbi(z, y, qe = qe, logspace = TRUE, type = type,
                    offset = offset, windowspace = windowspace, cpp = cpp)
    yind <- alig$path
    yind[alig$path != 0] <- 1:length(y)
    yind[alig$path == 0] <- 0 ### necessary?
    tmp <- matrix(gapchar) #single element matrix
    rownames(tmp) <- rownames(y)
    newy <- cbind(tmp, y)[,yind + 1]
    #also need to account for inserts in x
    ynotinsert <- alig$path != 2 #logical vector
    yinsertlengths <- insertlengths(ynotinsert) #tabulate insert lengths
    xinsertlengths <- z$insertlengths
    #reconcile x and y insert lengths
    # gls = gap lengths, gps = gap positions
    diffsx <- diffsy <- yinsertlengths - xinsertlengths
    newx <- x
    if(any(diffsx > 0)){
      diffsx[diffsx < 0] <- 0
      xgls <- diffsx[diffsx > 0]
      xgps <- c(z$alignment, ncol(x) + 1)[which(diffsx > 0)] - 1
      newx <- insertgaps(newx, xgps, xgls, gapchar = gapchar)
    }
    if(any(diffsy < 0)){
      diffsy[diffsy > 0] <- 0
      #align y to model
      #yprog <- alig$progression
      yprog <- progression(alig$path, alig$start)[1, ]
      ygls <- -1 * diffsy[diffsy < 0]
      ygps <- sapply(which(diffsy < 0), match, c(yprog, l + 1)) - 1
      newy <- insertgaps(newy, ygps, ygls, gapchar = gapchar)
    }
    res <- rbind(newx, newy)
    #rownames(res)[n + 1] <- tmp1
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    return(res)
  }else if(nrow(x) > 1 & nrow(y) > 1){
    #residues <- alphadetect(x, residues = residues, gapchar = gapchar)
    nx <- nrow(x)
    ny <- nrow(y)
    zx <- derive.PHMM(x, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    zy <- derive.PHMM(y, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    lx <- zx$size
    ly <- zy$size
    alig <- Viterbi(zx, zy, qe = qe, logspace = TRUE, type = type,
                    offset = offset, windowspace = windowspace, cpp = cpp)
    #vectors same length as model (+ 1 for begin state) with counts of
    #no of gaps to insert after each position
    xinsertlengths <- insertlengths(alig$path < 3) - zx$insertlengths
    yinsertlengths <- insertlengths(alig$path > 1) - zy$insertlengths
    #these vecs correspond to orig alig cols,
    #keeps track of how many gaps to insert after each
    resx <- rep(0, ncol(x) + 1)
    names(resx) <- 0:ncol(x)
    resy <- rep(0, ncol(y) + 1)
    names(resy) <- 0:ncol(y)
    zxali <- c(0, zx$alignment, ncol(x) + 1)
    zyali <- c(0, zy$alignment, ncol(y) + 1)
    names(zxali) <- 0:(lx + 1)
    names(zyali) <- 0:(ly + 1) # includes begin and end states
    #prog <- alig$progression
    prog <- progression2(alig$path, alig$start)
    if(!(all(prog[, 1] == 0))) prog <- cbind(c(0, 0), prog)
    ypve <- which(yinsertlengths > 0)
    ynve <- which(yinsertlengths < 0)
    xpve <- which(xinsertlengths > 0)
    xnve <- which(xinsertlengths < 0)
    #ptiga = positions to insert gaps after
    #ctiga = columns (of original alignment) to insert gaps after
    if(length(ypve) > 0){
      yctiga <- zyali[ypve + 1] - 1
      ygls <- yinsertlengths[ypve]
      resy[yctiga + 1] <- ygls
    }
    if(length(ynve) > 0){
      xptiga <- prog[1, sapply(ynve - 1, match, prog[2,])] ## #may need checking
      xctiga <- zxali[xptiga + 2] - 1
      xgls <- -1 * yinsertlengths[ynve]
      resx[xctiga + 1] <- xgls
    }
    if(length(xpve) > 0){
      xctiga <- zxali[xpve + 1] - 1
      xgls <- xinsertlengths[xpve]
      resx[xctiga + 1] <- resx[xctiga + 1] + xgls
    }
    if(length(xnve) > 0){
      yptiga <- prog[2, sapply(xnve - 1, match, prog[1,])]
      yctiga <- zyali[yptiga + 2] - 1
      ygls <- -1 * xinsertlengths[xnve]
      resy[yctiga + 1] <- resy[yctiga + 1] + ygls
    }
    newx <- insertgaps(x, which(resx > 0) - 1, resx[resx > 0], gapchar = gapchar)
    newy <- insertgaps(y, which(resy > 0) - 1, resy[resy > 0], gapchar = gapchar)
    res <- rbind(newx, newy)
    res <- res[, apply(res, 2, function(v) !all(v == gapchar))]
    class(res) <- if(DNA) "DNAbin" else if(AA) "AAbin" else NULL
    return(res)
  }else{
    stop("invalid arguments provided for x and or y")
  }
}

#' Align sequences to a profile HMM.
#'
#' Blah blah
#'
#' @param sequences a list of ...
#' @param model an object of class \code{"PHMM"}
#' @param ... further arguments to be passed to \code{"Viterbi"}.
#'
align2phmm <- function(sequences, model, gapchar = "-", ...){
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
