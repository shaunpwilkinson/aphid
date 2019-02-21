#' Plot profile hidden Markov models.
#'
#' \code{plot.PHMM} provides a visual representation of a profile hidden
#'   Markov model.
#'
#' @param x an object of class \code{"PHMM"}.
#' @param from an integer giving the module number to start the plot
#'   sequence from. Also accepts the chracter string "start" (module 0; default).
#' @param to an integer giving the module number to terminate the plot sequence.
#'   Also accepts the chracter string "end" (default).
#' @param just a character string giving the justfication of the plot relative
#'   to the device. Accepted values are "left", "center" and "right".
#' @param arrexp the expansion factor to be applied to the arrows in the plot.
#' @param textexp the expansion factor to be applied to the text in the plot.
#' @param ... additional arguments to be passed to \code{\link{plot}}.
#' @return NULL (invisibly).
#' @details \code{"plot.PHMM"} Plots a \code{"PHMM"} object as a directed graph
#'   with sequential modules consisting of squares, diamonds and circles
#'   representing match, insert and delete states, respectively.
#'   Modules are interconnected by directed
#'   lines with line-weights proportional to the transition probabilities between
#'   the states. Since the plotted models are generally much longer than they are
#'   high, it is usually better to output the plot to a PDF file as demonstrated
#'   in the example below.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @seealso \code{\link{plot.HMM}}
#' @examples
#'   ## Small globin alignment example from Durbin et al (1998) Figure 5.3
#'   data(globins)
#'   ## derive a profile hidden Markov model from the alignment
#'   globins.PHMM <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
#'   ## plot the PHMM
#'   plot(globins.PHMM, main = "Profile hidden Markov model for globins")
#'   ##
#'   ## derive a profile hidden Markov model from the woodmouse dataset in the
#'   ## ape package
#'   library(ape)
#'   data(woodmouse)
#'   woodmouse.PHMM <- derivePHMM(woodmouse)
#'   ## plot partial model to viewer device
#'   plot(woodmouse.PHMM, from = 0, to = 5)
#'   ## plot the entire model to a PDF in the current working directory
#'   \donttest{
#'   tmpf <- tempfile(fileext = ".pdf")
#'   nr <- ceiling((woodmouse.PHMM$size + 2)/10)
#'   pdf(file = tmpf, width = 8.27, height = nr * 2)
#'   par(mfrow = c(nr, 1), mar = c(0, 0, 0, 0) + 0.1)
#'   from <- 0
#'   to <- 10
#'   for(i in 1:nr){
#'     plot(woodmouse.PHMM, from = from, to = to, just = "left")
#'     from <- from + 10
#'     to <- min(to + 10, woodmouse.PHMM$size + 1)
#'   }
#'   dev.off()
#'   }
################################################################################
plot.PHMM <- function(x, from = "start", to = "end", just = "center",
                      arrexp = 1, textexp = 1, ...){
  logspace <- .logdetect(x)
  if(logspace){
    x$A <- exp(x$A)
    x$E <- exp(x$E)
  }
  plot(0:1, 0:1, type = 'n', axes = FALSE, ylab = "", xlab = "", ... = ...)
  parusr <- par()$usr
  residues <- rownames(x$E)
  symblen <- length(residues)
  pHMMlength <- ncol(x$E)
  #maxemiss <- if(ncol(x$E) > 0) max(x$E) else 1
  maxemiss <- 1
  if(from == "start") from <- 0
  if(to == "end") to <- pHMMlength + 1
  no.boxes <- to  + 1 - from
  boxhgt <- symblen/4 # in yunits
  plotwid <- 2 * no.boxes - 1 # in xunits
  plothgt <- 4 + boxhgt # in yunits
  plotratio <- plothgt/plotwid
  canvasdims <- par()$pin #inches
  canwid <- canvasdims[1]/(parusr[2] - parusr[1]) #inches
  canhgt <- canvasdims[2]/(parusr[4] - parusr[3]) #inches
  canratio <- canhgt/canwid
  lim <- if(canratio > plotratio) 'wid' else 'hgt' # is plot width or height limited
  scalefac <- if(lim == 'wid') canwid/plotwid else canhgt/plothgt
  # one box is *scalefac* inches wide (and high)
  xunit <- scalefac/canwid
  yunit <- scalefac/canhgt
  boxhgt <- boxhgt * yunit
  plotwid <- plotwid * xunit
  plothgt <- plothgt * yunit
  coords <- matrix(nrow = 3 * no.boxes, ncol = 2)
  coords[, 1] <- rep(seq(xunit/2, xunit/2 + (no.boxes * 2 - 1) * xunit,
                         by = 2 * xunit), each = 3)
  if(just == 'center') coords[, 1] <- coords[, 1] + (1 - plotwid)/2
  if(just == 'right') coords[, 1] <- coords[, 1] + (1 - plotwid)
  ccy <- 0.5 + plothgt/2 - yunit/2 #circle centre y coord
  dcy <- ccy - 2 * yunit # diamond centre y coord
  bcy <- dcy - 1.5 * yunit - boxhgt/2 # box centre y coord
  coords[, 2] <- rep(c(ccy, dcy, bcy), no.boxes)# - (1 - plothgt)/2l
  fromto <- cbind(rep(1:(3 * (pHMMlength + 1)), each = 3),
                  rep(rep(c(2, 4, 6), 3), pHMMlength + 1) +
                    rep(seq(0, 3 * pHMMlength, 3), each = 9))
  nr <- nrow(fromto)
  transitions <- rep(c("DI", "DD", "DM", "II", "ID", "IM", "MI", "MD", "MM"), pHMMlength + 1)
  modstates <- paste(rep(0:pHMMlength, each = 9))
  nullarrows <- rep(FALSE, nr)
  nullarrows[c(1:3, nr - c(1, 4, 7))] <- TRUE
  if(from != 0) nullarrows[1:(from * 9)] <- TRUE #9 diff arrows from each pos
  if(to != pHMMlength + 1) nullarrows[((to * 9 + 2):nr)[-c(3, 6)]] <- TRUE
  for(i in (1:nr)[!nullarrows]){
    arrwgt <- x$A[transitions[i], modstates[i]] * 4 * arrexp # factor 4 looks abt right
    fromi <- fromto[i - 9 * from, 1]
    toi <- fromto[i - 9 * from, 2]
    coordsi <- coords[c(fromi, toi), ]
    ## coordsi[1,1] is x0, [1,2] is y0, [2,1] is x1, [2,2] is y1
    if(fromi %% 3 == 1 & toi %% 3 == 2){
      arrows(x0 = coordsi[1, 1], y0 = coordsi[1, 2],
             x1 = coordsi[2, 1], y1 =coordsi[2, 2] + sqrt((yunit^2)/2),
             length = scalefac/10, lwd = arrwgt)# fix length
    }else if(fromi %% 3 == 0 & toi %% 3 == 2){
      arrows(x0 = coordsi[1, 1], y0 = coordsi[1, 2],
             x1 = coordsi[2, 1], y1 =coordsi[2, 2] - sqrt((yunit^2)/2),
             length = scalefac/10, lwd = arrwgt)# fix length
    }else{
      lines(coordsi, lwd = arrwgt)
    }
  }
  pos <- from
  symbcoords <- matrix(nrow = symblen, ncol = 2)
  symbcoords[, 1] <- if(pos == 0) coords[6, 1] else coords[3, 1]
  symbcoords[, 1] <- symbcoords[, 1] - xunit * 0.4
  interval <- boxhgt/(symblen + 1)
  tmpseq <- seq(from = coords[3, 2] + boxhgt/2,
                to = coords[3, 2] - boxhgt/2,
                by = -interval)
  symbcoords[,2] <- tmpseq[-c(1, length(tmpseq))]
  for(i in seq(from = 1, to = 3 * no.boxes, by = 3)){
    if(pos != pHMMlength + 1){
      .diamond(coords[i + 1, 1], coords[i + 1, 2],
              radx = sqrt((xunit^2)/2),
              rady = sqrt((yunit^2)/2),
              col = 'white')
      text(x = coords[i + 1, 1], y = coords[i + 1, 2],
           labels = paste(round(x$A["II", pos + 1] * 100, 0)),
           cex = textexp)
      if(pos != 0){
        .ellipse(coords[i, 1], coords[i, 2],
                radx = xunit/2,
                rady = yunit/2,
                col = 'white')
        text(x = coords[i, 1], y = coords[i, 2], labels = paste(pos),
             cex = textexp)
        rect(xleft = coords[i + 2, 1] - xunit/2,
             ybottom = coords[i + 2, 2] - boxhgt/2,
             xright = coords[i + 2, 1] + xunit/2,
             ytop = coords[i + 2, 2] + boxhgt/2,
             col = 'white')
        text(symbcoords, labels = residues, adj = 0, cex = 0.4 * textexp)
        for(j in 1:symblen){
          rect(xleft = symbcoords[j, 1] + xunit/3,
               ybottom = symbcoords[j, 2] - 0.4 * interval,
               xright = symbcoords[j, 1] + (x$E[j, pos] *
                                              xunit/2)/maxemiss + xunit/3,
               ytop = symbcoords[j, 2] + 0.4 * interval,
               lwd = 0,
               col = 'black')
        }
        symbcoords[, 1] <- symbcoords[, 1] + (2 * xunit)
      }
    }
    if(pos == 0 | pos == pHMMlength + 1){
      rect(xleft = coords[i + 2, 1] - xunit/2,
           ybottom = coords[i + 2, 2] - yunit/2,
           xright = coords[i + 2, 1] + xunit/2,
           ytop = coords[i + 2, 2] + yunit/2,
           col = 'white')
      text(x = coords[i + 2, 1], y = coords[i + 2, 2],
           labels = if(pos == 0) "B" else "E", cex = textexp)
    }
    pos <- pos + 1
  }
  invisible()
}
################################################################################
#' Plot standard hidden Markov models.
#'
#' \code{plot.HMM} provides a visual representation of a standard hidden Markov
#'   model.
#'
#' @param x an object of class \code{"HMM"}.
#' @param begin logical indicating whether the begin/end state should be plotted.
#'   Defaults to FALSE.
#' @inheritParams plot.PHMM
#' @return NULL (invisibly).
#' @details \code{"plot.HMM"} Plots a \code{"HMM"} object as a directed graph.
#'   States (rectangles) are interconnected by directed
#'   lines with line-weights proportional to the transition probabilities between
#'   the states.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @seealso \code{\link{plot.PHMM}}
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @examples
#'   ## the dishonest casino example from Durbin et al (1998)
#'   states <- c("Begin", "Fair", "Loaded")
#'   residues = paste(1:6)
#'   A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9), nrow = 3)
#'   dimnames(A) <- list(from = states, to = states)
#'   E <- matrix(c(rep(1/6, 6), rep(1/10, 5), 1/2), nrow = 2, byrow = TRUE)
#'   dimnames(E) <- list(states = states[-1], residues = residues)
#'   x <- structure(list(A = A, E = E), class = "HMM")
#'   plot(x, main = "Dishonest casino hidden Markov model")
################################################################################
plot.HMM <- function(x, just = "center", arrexp = 1, textexp = 1,
                     begin = FALSE, ...){
  logspace <- .logdetect(x)
  if(logspace){
    x$A <- exp(x$A)
    x$E <- exp(x$E)
  }
  plot(0:1, 0:1, type = 'n', axes = FALSE, ylab = "", xlab = "", ... = ...)
  parusr <- par()$usr
  states <- rownames(x$E)
  residues <- colnames(x$E)
  statelen <- length(states)
  symblen <- length(residues)
  #maxemiss <- max(x$E)
  maxemiss <- 1
  no.boxes <- if(begin) statelen + 1 else statelen
  boxhgt <- symblen/4 # in yunits
  plotwid <- 2 * no.boxes # in xunits
  plothgt <- (no.boxes - 1) * boxhgt # in yunits
  plotratio <- plothgt/plotwid
  canvasdims <- par()$pin #inches
  canwid <- canvasdims[1]/(parusr[2] - parusr[1]) #inches
  canhgt <- canvasdims[2]/(parusr[4] - parusr[3]) #inches
  canratio <- canhgt/canwid
  lim <- if(canratio > plotratio) 'wid' else 'hgt' # is plot width or height limited
  scalefac <- if(lim == 'wid') canwid/plotwid else canhgt/plothgt
  xunit <- scalefac/canwid
  yunit <- scalefac/canhgt
  boxhgt <- boxhgt * yunit
  plotwid <- plotwid * xunit
  plothgt <- plothgt * yunit
  coords <- matrix(nrow = no.boxes, ncol = 2)
  coords[, 1] <- seq(xunit, (no.boxes * 2 - 1) * xunit,
                     by = 2 * xunit)
  if(just == 'center') coords[, 1] <- coords[, 1] + (1 - plotwid)/2
  if(just == 'right') coords[, 1] <- coords[, 1] + (1 - plotwid)
  coords[, 2] <- 0.5
  if(begin){
    coords1 <- coords[1,]
    coords <- coords[-1,]
    #self arrow
    arrwgt0 <- x$A[1, 1] * 4 * arrexp
    .arc(x = coords1[1] - xunit/2, y = coords1[2],
        from = pi/2, to = 5 * pi/2,
        radx = xunit/4, rady = yunit/4,
        arrow = 0.5, lwd = arrwgt0)
    for(i in 1:statelen){
      arrwgt1 <- x$A[1, i + 1] * 4 * arrexp
      arrwgt2 <- x$A[i + 1, 1] * 4 * arrexp
      .arc(x = coords1[1] + i * xunit, y = coords1[2],
          from = -pi/2, to = pi/2,
          radx = i * xunit, rady = i * boxhgt/2,
          arrow = 0.5, lwd = arrwgt1)
      .arc(x = coords1[1] + i * xunit, y = coords1[2],
          from = pi/2, to = 3*pi/2,
          radx = i * xunit, rady = i * boxhgt/2,
          arrow = 0.5, lwd = arrwgt2)
    }
    rect(xleft = coords1 [1] - xunit/2,
         ybottom = coords1[2] - boxhgt/2,
         xright = coords1[1] + xunit/2,
         ytop = coords1[2] + boxhgt/2,
         col = 'white')
    text(x = coords1[1], y = coords1[2], labels = c("BEGIN\nEND"),
         cex = 0.6 * textexp)
  }
  for(i in 1:statelen){
    for(j in 1:statelen){
      arrwgt <- x$A[i + 1, j + 1] * 4 * arrexp
      if(i == j){
        .arc(x = coords[i, 1] - xunit/2, y = coords[i, 2],
            from = pi/2, to = 5 * pi/2,
            radx = xunit/4, rady = yunit/4,
            arrow = 0.5, lwd = arrwgt)
      }else{
        dif <- abs(diff(c(i, j)))
        .arc(x = mean(coords[c(i, j),1]), y = coords[i, 2],
            from = if(j > i) -pi/2 else pi/2,
            to = if(j > i) pi/2 else 3 * pi/2,
            radx = dif * xunit, rady = dif * boxhgt/2,
            arrow = 0.5, lwd = arrwgt)
      }
    }
  }
  interval <- boxhgt/(symblen + 2)
  symbcoords <- matrix(nrow = symblen + 1, ncol = 2)
  symbcoords[, 1] <- coords[1, 1]
  symbcoords[, 1] <- symbcoords[, 1] - xunit * 0.4
  tmpseq <- seq(from = coords[1, 2] + boxhgt/2,
                to = coords[1, 2] - boxhgt/2,
                by = -interval)
  symbcoords[,2] <- tmpseq[-c(1, length(tmpseq))]
  for(i in 1:statelen){
    rect(xleft = coords[i, 1] - xunit/2,
         ybottom = coords[i, 2] - boxhgt/2,
         xright = coords[i, 1] + xunit/2,
         ytop = coords[i, 2] + boxhgt/2,
         col = 'white')
    text(symbcoords, labels = c(states[i], residues), adj = 0,
         cex = 0.6 * textexp)
    for(j in 1:symblen){
      rect(xleft = symbcoords[j + 1, 1] + xunit/4,
           ybottom = symbcoords[j + 1, 2] - 0.4 * interval,
           xright = symbcoords[j + 1, 1] +
             (x$E[i, j] * xunit/2)/maxemiss + xunit/4,
           ytop = symbcoords[j + 1, 2] + 0.4 * interval,
           lwd = 0,
           col = 'black')
    }
    symbcoords[, 1] <- symbcoords[, 1] + (2 * xunit)
  }
  invisible()
}
################################################################################
#################### Internal geometric functions ##############################
################################################################################
.arc <- function(x, y, radx, rady = radx, from = 0, to = 2 * pi,
                no.points = 100, fill = NULL,
                arrow = NULL, arrowsize = 0.08, code = 2, ...){
  piseq <- seq(from, to, by = (to - from)/no.points)
  coords <- matrix(nrow = no.points + 1, ncol = 2)
  coords[, 1] <- x + radx * sin(piseq)
  coords[, 2] <- y + rady * cos(piseq)
  if(from == 0 & to == 2 * pi){
    polygon(coords, col = fill, ... = ...)
  }else{
    lines(coords, ... = ...)
  }
  if(!(is.null(arrow))){
    if(arrow < 0 | arrow > 1) stop ("arrow argument must be between 0 and 1")
    if(arrow == 0){
      arrows(x0 = coords[1, 1], y0 = coords[1, 2], x1 = coords[2, 1],
             y1 = coords[2, 2], code = 1, length = arrowsize, ... = ...)
    }else{
      l <- ceiling(arrow * no.points)
      arrows(x0 = coords[l, 1], y0 = coords[l, 2], x1 = coords[l + 1, 1],
             y1 = coords[l + 1, 2], code = code, length = arrowsize, ... = ...)
    }
  }
  invisible()
}
################################################################################
.chord <- function(x, y, rad, type = 'outer', no.points = 100,
                  arrow = TRUE, arrowlength = 0.08, reversearrow = FALSE,
                  ...){#x and y are vectors of from, to
  distance <- sqrt(diff(x)^2 + diff(y)^2)
  tmp <- sqrt(rad^2 - (distance/2)^2) #distance from midpoint to centre
  midpoint <- c(mean(x), mean(y)) #coordinates
  angle <- atan2(diff(x), diff(y))
  opp <- sin(angle) * tmp
  adj <- cos(angle) * tmp
  centcoords <- midpoint + (c(adj, -opp))
  segangle <- 2 * acos(tmp/rad)
  startang <- atan2(x[1]- centcoords[1], y[1] - centcoords[2])
  if(startang < 0) startang <- -startang + (2 * (pi + startang))
  if(type == 'outer'){
    piseq <- seq(startang, startang + segangle, by = segangle/no.points)
  } else if (type == 'inner'){
    piseq <- seq(startang, startang - (2 * pi - segangle),
                 by = -(2 * pi - segangle)/no.points)
  } else{
    stop("type argument must be set to either 'outer' or 'inner'")
  }
  coords <- matrix(nrow = no.points + 1, ncol = 2)
  coords[, 1] <- centcoords[1] + rad * sin(piseq)
  coords[, 2] <- centcoords[2] + rad * cos(piseq)
  lines(coords, ... = ...)
  if(arrow){
    arrowdir <- if(reversearrow) -1 else 1
    arrows(x0 = coords[no.points/2 , 1],
           y0 = coords[no.points/2, 2],
           x1 = coords[no.points/2 + arrowdir, 1],
           y1 = coords[no.points/2 + arrowdir, 2],
           length = arrowlength,
           ... = ...)
  }
  invisible()
}
################################################################################
.ellipse <- function(x, y, radx, rady = radx,
                    no.points = 100, col = NULL, lwd = 1){
  piseq <- seq(0, 2 * pi, by = (2 * pi)/no.points)
  coords <- matrix(nrow = no.points + 1, ncol = 2)
  coords[, 1] <- x + radx * sin(piseq)
  coords[, 2] <- y - rady * cos(piseq)
  polygon(coords, col = col, lwd = lwd)
  invisible()
}
################################################################################
.diamond <- function(x, y, radx, rady = radx, col = NULL, lwd = 1){
  xcoords <- c(x + c(0, radx, 0, -radx))
  ycoords <- c(y + c(rady, 0, -rady, 0))
  polygon(xcoords, ycoords, col = col, lwd = lwd)
  invisible()
}
################################################################################
