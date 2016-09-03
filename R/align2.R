#' Multiple sequence alignment.
#'
#' \code{align} finds the optimal alignment for a list of sequences using a hybrid
#' algorithm that involves a progressive alignment step followed by an iterative refinement
#' stage.
#'
#' @param sequences a list of character vectors consisting of symbols from
#' the residue alphabet
#' @param gapchar the character used to represent gaps in the alignment matrix.
#' @param residues one of the character strings "autodetect" (default), "aminos" or
#' "bases", or a case sensitive character vector matching the residue
#' alphabet. Specifying the symbol type ('aminos' or 'bases') can increase speed
#' for larger alignments. Note that setting \code{residues = "autodetect"} will not
#' detect rare residues that are not present in the alignment and thus will
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
#'
#'
align <- function(sequences, type = "global", residues = "autodetect",
                  gapchar = "-", DI = FALSE, ID = FALSE, refine = "Viterbi",
                  quiet = TRUE, ...){
  if(!(is.list(sequences))) stop("invalid 'sequences' agrument")
  nsq <- length(sequences)
  if(is.null(attr(sequences, "names"))) names(sequences) <- paste0("SEQ", 1:nsq)
  # seqlengths <- sapply(sequences, length)
  # nmodules <- round(mean(seqlengths))
  residues <- alphabet(sequences, residues = residues, gapchar = gapchar)
  fun <- Vectorize(function(i, j) ktup(sequences[[i]], sequences[[j]]))
  qds <- outer(1:nsq, 1:nsq, fun)
  ### need to speed this up in C++
  dimnames(qds) <- list(names(sequences),names(sequences))
  guidetree <- as.dendrogram(hclust(as.dist(qds), method = "average"))
  newick <- write.dendrogram(guidetree, strip.edges = TRUE)
  newick <- gsub(";", "", newick)
  newick <- gsub("\\(", "alignpair\\(", newick)
  if(type == 'semiglobal') newick <- gsub("\\)", ", type = 'semiglobal'\\)", newick)
  msa1 <- with(sequences, eval(parse(text = newick)))
  omniphmm <- derivePHMM(msa1, DI = DI, ID = ID, pseudocounts = "background")
  if(refine == "Viterbi"){
    finalphmm <- train(omniphmm, sequences, maxiter = 300, DI = DI, ID = ID,
                       inserts = "map", quiet = quiet, ... = ...)
  }else if(refine == "BaumWelch"){
    finalphmm <- BaumWelch(omniphmm, sequences, maxiter = 300, DI = DI, ID = ID,
                           quiet = quiet, ... = ...)
  } else if (refine == "none"){
    return(msa1)
  } else stop("argiument 'refine' must be either 'Viterbi', 'BaumWelch' or 'none'")
  align2phmm(sequences, model = finalphmm)
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
#' @return a character matrix of aligned sequences.
#' @references Soding...
#' @examples
#' x <- c("H", "E", "A", "G", "A", "W", "G", "H", "E", "E")
#' y <- c("P", "A", "W", "H", "E", "A", "E")
#' z <- align(x, y)
#' alignpair(x, z)
#'
alignpair <- function(x, y, d = 8, e = 2, S = NULL, qe = NULL,
                      type = "global",
                      offset = -0.1, itertab = NULL, pseudocounts = "background",
                      residues = "autodetect", gapchar = "-"){
  if(is.vector(x) & is.vector(y)){
    alig <- Viterbi(x, y, d = d, e = e, S = S, type = type,
                    itertab = itertab)#, offset = offset) ###not necessary for vec vs vec
    xind <- yind <- alig$path
    xind[alig$path != 3] <- 1:length(x)
    xind[alig$path == 3] <- 0
    newx <- c(gapchar, x)[xind + 1]
    yind[alig$path != 1] <- 1:length(y)
    yind[alig$path == 1] <- 0
    newy <- c(gapchar, y)[yind + 1]
    res <- rbind(newx, newy)
    rownames(res) <- c(deparse(substitute(x)), deparse(substitute(y)))
    return(res)
  }else if((is.matrix(x) & is.vector(y)) | (is.vector(x) & is.matrix(y))){ ###also need option to flip
    vm <- is.vector(x) & is.matrix(y)
    if(vm){
      tmp1 <- deparse(substitute(x))
      tmp2 <- x
      x <- y
      y <- tmp2
    }else{
      tmp1 <- deparse(substitute(y))
    }
    if(identical(residues, "autodetect")) residues <- sort(unique(c(as.vector(x), y)))
    n <- nrow(x)
    z <- derivePHMM(x, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    l <- z$size
    alig <- Viterbi(z, y, qe = qe, logspace = TRUE, type = type,
                    offset = offset, itertab = itertab)
    yind <- alig$path
    yind[alig$path != 1] <- 1:length(y)
    yind[alig$path == 1] <- 0
    newy <- c(gapchar, y)[yind + 1]
    #also need to account for inserts in x
    ynotinsert <- alig$path != 3 #logical vector
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
      yprog <- alig$progression
      ygls <- -1 * diffsy[diffsy < 0]
      ygps <- sapply(which(diffsy < 0), match, c(yprog, l + 1)) - 1
      newy <- insertgaps(newy, ygps, ygls, gapchar = gapchar)
    }
    res <- rbind(newx, newy)
    rownames(res)[n + 1] <- tmp1
    return(res)
  }else if(is.matrix(x) & is.matrix(y)){
    if(identical(residues, "autodetect")) residues <- sort(unique(c(as.vector(x), as.vector(y))))
    nx <- nrow(x)
    ny <- nrow(y)
    zx <- derivePHMM(x, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    zy <- derivePHMM(y, pseudocounts = pseudocounts, residues = residues, logspace = TRUE)
    lx <- zx$size
    ly <- zy$size
    alig <- Viterbi(zx, zy, qe = qe, logspace = TRUE, type = type,
                    offset = offset, itertab = itertab)
    #vectors same length as model (+ 1 for begin state) with counts of
    #no of gaps to insert after each position
    xinsertlengths <- insertlengths(alig$path < 4) - zx$insertlengths
    yinsertlengths <- insertlengths(alig$path > 2) - zy$insertlengths
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
    prog <- alig$progression
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
    return(res)
  }else{
    stop("invalid arguments provided for x and or y")
  }
}

#' Align sequences to a profile HMM.
#'
#' @param ... further arguments to be passed to \code{Viterbi}.
#'
align2phmm <- function(sequences, model, gapchar = "-", ...){
  if(is.list(sequences)){
  }else if(is.vector(sequences, mode = "character")){
    sequences <- list(sequences)
  }else stop("invalid 'sequences' argument")
  out <- matrix(nrow = length(sequences), ncol = 2 * model$size + 1)
  rownames(out) <- attr(sequences, "names")
  colnames(out) <- rep("I", 2 * model$size + 1)
  colnames(out)[seq(2, 2 * model$size, by = 2)] <- 1:model$size
  for(i in 1:length(sequences)){
    alignment <- Viterbi(model, sequences[[i]], type = 'global', ...)
    newrow <- rep(gapchar, length(alignment$path))
    newrow[alignment$path > 1] <- sequences[[i]]
    inserts <- alignment$path == 3
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
      tuples <- tuples[, decs != 2] #remove true->falses
      tuples <- tuples[,-ncol(tuples)]
      tuples[1, -1] <- TRUE
    }else{
      tuples <- rbind(c(FALSE, rep(TRUE, model$size)), rep(FALSE, model$size + 1))
    }
    indices <- as.vector(tuples)[-1]
    newrow2 <- rep(gapchar, 2 * model$size + 1)
    names(newrow2) <-  c(rep(c("I", "M"), model$size), "I") #seq(0.5, x$size + 0.5, by = 0.5)
    newrow2[indices] <- newrow
    out[i,] <- newrow2
  }
  discardcols <- apply(out, 2 ,function(v) all(v == gapchar))
  matchcols <- c(FALSE, rep(c(TRUE, FALSE), model$size))
  out <- out[, matchcols | !discardcols]
  out <- as.list(as.data.frame(out, stringsAsFactors = F))
  fun <- function(e){
    ee <- strsplit(e, split = "")
    elengths <- sapply(ee, length)
    maxlen <- max(elengths)
    no.gapstoappend <- maxlen - elengths
    appges <- lapply(no.gapstoappend, function(v) rep(gapchar, v))
    res <- t(mapply(c, ee, appges))
    res
  }
  out[names(out) == "I"] <- lapply(out[names(out) == "I"], fun)
  out <- as.matrix(as.data.frame(out, stringsAsFactors = F, optional = T))
  inserts <- sapply(colnames(out), function(v) grepl("I", v))
  colnames(out)[inserts] <- "I"
  rownames(out) <- names(sequences)
  return(out)
}
