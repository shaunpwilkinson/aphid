plot.PHMM <- function(x, from = "start", to = "end", just = "center",
                      arrexp = 1, textexp = 1, ...){
  plot(0:1, 0:1, type = 'n', axes = FALSE, ylab = "", xlab = "", ... = ...)
  parusr <- par()$usr
  residues <- rownames(x$E)
  symblen <- length(residues)
  pHMMlength <- ncol(x$E)
  maxemiss <- max(x$E)
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
  lim <- if(canratio > plotratio) 'wid' else 'hgt' #is plot width or height limited?
  scalefac <- if(lim == 'wid') canwid/plotwid else canhgt/plothgt
  # so one box is *scalefac* inches wide (and high)
  xunit <- scalefac/canwid
  yunit <- scalefac/canhgt
  boxhgt <- boxhgt * yunit
  plotwid <- plotwid * xunit
  plothgt <- plothgt * yunit
  ### so how many inches is one yunit?

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
  fromstates <- rep(rep(c("D", "I", "M"), each = 3), times = pHMMlength + 1)
  tostates <- rep(c("I", "D", "M"), times = 3 * (pHMMlength + 1))
  modstates <- paste(rep(0:pHMMlength, each = 9))
  nullarrows <- rep(FALSE, nr)
  nullarrows[c(1:3, nr - c(1, 4, 7))] <- TRUE
  if(from != 0) nullarrows[1:(from * 9)] <- TRUE #9 diff arrows from each pos
  if(to != pHMMlength + 1) nullarrows[((to * 9 + 2):nr)[-c(3, 6)]] <- TRUE
  for(i in (1:nr)[!nullarrows]){
    arrwgt <- x$A[fromstates[i], modstates[i], tostates[i]] * 4 * arrexp
    ## 4 just looks about right
    # lines(rbind(coords[fromto[i - 9 * from, 1], ],
    #             coords[fromto[i - 9 * from, 2], ]), lwd = arrwgt)
    fromi <- fromto[i - 9 * from, 1]
    toi <- fromto[i - 9 * from, 2]
    coordsi <- coords[c(fromi, toi), ]
    #coordsi[1,1] is x0, [1,2] is y0, [2,1] is x1, [2,2] is y1
    if(fromi %% 3 == 1 & toi %% 3 == 2){
      arrows(x0 = coordsi[1, 1], y0 = coordsi[1, 2],
             x1 = coordsi[2, 1], y1 =coordsi[2, 2] + sqrt((yunit^2)/2),
             length = scalefac/10, lwd = arrwgt)####### fix length
    }else if(fromi %% 3 == 0 & toi %% 3 == 2){
      arrows(x0 = coordsi[1, 1], y0 = coordsi[1, 2],
             x1 = coordsi[2, 1], y1 =coordsi[2, 2] - sqrt((yunit^2)/2),
             length = scalefac/10, lwd = arrwgt)####### fix length
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
      diamond(coords[i + 1, 1], coords[i + 1, 2],
              radx = sqrt((xunit^2)/2),
              rady = sqrt((yunit^2)/2),
              col = 'white')
      text(x = coords[i + 1, 1], y = coords[i + 1, 2],
           labels = paste(round(x$A["I", pos + 1, "I"] * 100, 0)),
           cex = textexp)
      if(pos != 0){
        ellipse(coords[i, 1], coords[i, 2],
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
               col = 'darkgrey')
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
}

plot.HMM <- function(x, just = "center", arrexp = 1, textexp = 1,
                     includebegin = FALSE, ...){
  plot(0:1, 0:1, type = 'n', axes = FALSE, ylab = "", xlab = "")#, ... = ...)
  parusr <- par()$usr
  states <- rownames(x$E)
  residues <- colnames(x$E)
  statelen <- length(states)
  symblen <- length(residues)
  maxemiss <- max(x$E)
  no.boxes <- if(includebegin) statelen + 1 else statelen
  boxhgt <- symblen/4 # in yunits
  plotwid <- 2 * no.boxes # in xunits
  plothgt <- (no.boxes - 1) * boxhgt # in yunits
  plotratio <- plothgt/plotwid
  canvasdims <- par()$pin #inches
  canwid <- canvasdims[1]/(parusr[2] - parusr[1]) #inches
  canhgt <- canvasdims[2]/(parusr[4] - parusr[3]) #inches
  canratio <- canhgt/canwid
  lim <- if(canratio > plotratio) 'wid' else 'hgt' #is plot width or height limited?
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
  if(includebegin){
    coords1 <- coords[1,]
    coords <- coords[-1,]
    #self arrow
    arrwgt0 <- x$A[1, 1] * 4 * arrexp
    arc(x = coords1[1] - xunit/2, y = coords1[2],
        from = pi/2, to = 5 * pi/2,
        radx = xunit/4, rady = yunit/4,
        arrow = 0.5, lwd = arrwgt0)
    for(i in 1:statelen){
      arrwgt1 <- x$A[1, i + 1] * 4 * arrexp
      arrwgt2 <- x$A[i + 1, 1] * 4 * arrexp
      arc(x = coords1[1] + i * xunit, y = coords1[2],
          from = -pi/2, to = pi/2,
          radx = i * xunit, rady = i * boxhgt/2,
          arrow = 0.5, lwd = arrwgt1)
      arc(x = coords1[1] + i * xunit, y = coords1[2],
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
        arc(x = coords[i, 1] - xunit/2, y = coords[i, 2],
            from = pi/2, to = 5 * pi/2,
            radx = xunit/4, rady = yunit/4,
            arrow = 0.5, lwd = arrwgt)
      }else{
        dif <- abs(diff(c(i, j)))
        arc(x = mean(coords[c(i, j),1]), y = coords[i, 2],
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
           col = 'darkgrey')
    }
    symbcoords[, 1] <- symbcoords[, 1] + (2 * xunit)
  }
}

