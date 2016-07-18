# #old matrix vs matrix code


# function progress takes an integer vector and returns a sequence starting from 0
# where a zero progresses 1, a negative number progresses by its value
# and a positive removes its value of following integers 
progress <- function(v){
  vl <- as.list(v)
  vl[which(v < 0)] <- lapply(vl[which(v < 0)], function(x) c(0, rep(x, -1 * x)))
  vl[which(v > 0)] <- lapply(vl[which(v > 0)], function(x) c(0, rep(NA, x)))
  vlv <- unlist(vl)
  res <- c(vlv, 0)
  res[!is.na(vlv)] <- 0:length(vlv[!is.na(vlv)])
  res <- res[vlv >= 0]
  res
}
# v <- c(0, 0, -1, 0, 2, 0, -2, 0)
# progress(v) ### 0  1  2  4  5 NA NA  6  7 10




ofset <- -1 * min(c(xinsertlengths, yinsertlengths))
xil2 <- xinsertlengths + ofset
yil2 <- yinsertlengths + ofset
insertzeros <- function(x, y) c(rep(1, x), rep(0, y))

x2 <- unlist(mapply(insertzeros, zx$insertlengths + 1, xil2, SIMPLIFY = FALSE))
y2 <- unlist(mapply(insertzeros, zy$insertlengths + 1, yil2, SIMPLIFY = FALSE))

#xinsertlengths <- c(xinsertlengths, 0) # include end state as well
#yinsertlengths <- c(yinsertlengths, 0)

#progression along alignment (not model) 
#px and py should be of equal length (same as alig$path ?)
#first element is 0, last is ncol of alignment +  1 (includes imaginary begin and end cols)
px <- progress(xinsertlengths)
py <- progress(yinsertlengths)
#replace NAs with current state
px <- unlist(mapply(rep, px[!is.na(px)], insertlengths(!is.na(px))[-1] + 1, SIMPLIFY = FALSE))
py <- unlist(mapply(rep, py[!is.na(py)], insertlengths(!is.na(py))[-1] + 1, SIMPLIFY = FALSE)) 
vx <- px[-1] - px[-length(px)] - 1
vy <- py[-1] - py[-length(py)] - 1

# vx[vx < 0] <- 0
# vy[vy < 0] <- 0
# tmpx <- vx - vy
# vx[tmpx >= 0] <- tmpx[tmpx >= 0]
# vy[tmpx >= 0] <- 0
# tmpy <- vy - vx
# vy[tmpy >= 0] <- tmpy[tmpy >= 0]
# vx[tmpy >= 0] <- 0

pxl <- as.list(px)
pyl <- as.list(py)
func2 <- function(x) if(x > 0) rep(NA, x) else NULL
pxl2 <- lapply(vy, func2) 
pyl2 <- lapply(vx, func2)
func3 <- function(x, y) if(y > 0) seq(x + 1, x + y) else NULL
pxl3 <- mapply(func3, px[-length(px)], vx, SIMPLIFY = FALSE)
pyl3 <- mapply(func3, py[-length(py)], vy, SIMPLIFY = FALSE)
px <- unlist(mapply(c, pxl[-length(pxl)], pxl2, pxl3, SIMPLIFY = FALSE))
py <- unlist(mapply(c, pyl[-length(pyl)], pyl2, pyl3, SIMPLIFY = FALSE))
px <- px[-1]
py <- py[-1]
px[is.na(px)] <- 1
py[is.na(py)] <- 1
#func <- function(x, y) c(x, rep(NA, y))
#pxl[which(vy > 0)] <- mapply(func, pxl[which(vy > 0)], vy[vy > 0], SIMPLIFY = FALSE) 
#pyl[which(vx > 0)] <- mapply(func, pyl[which(vx > 0)], vx[vx > 0], SIMPLIFY = FALSE) 
newx <- cbind(gapchar, x)[, px]
newy <- cbind(gapchar, y)[, py]
res <- rbind(newx, newy)
res <- res[, apply(res, 2, function(v) !all(v == gapchar))] # can fix this earlier on
return(res)