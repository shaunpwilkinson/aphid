library(aphid)
context("build, train and apply profile HMMs")

# simulate a DNA sequence dataset
set.seed(999)
bases <- c("A", "C", "G", "T")
x <- list(sample(bases, replace = TRUE, size = 100))
evolve <- function(a) if(runif(1) > 0.95) sample(bases, 1) else a
for(i in 2:10) x[[i]] <- unname(sapply(x[[i - 1]], evolve))
names(x) <- paste("Sequence", 1:10)
# convert to DNAbin object
rawbases <- as.raw(c(136, 40, 72, 24))
xDNA <- lapply(x, function(s) rawbases[match(s, bases)])
class(xDNA) <- "DNAbin"

# simulate an AA sequence dataset, this time an alignment
set.seed(999)
aminos <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
y <- matrix(sample(aminos, replace = TRUE, size = 100), nrow = 1)
evolve <- function(a) if(runif(1) > 0.95) sample(aminos, 1) else a
for(i in 2:10) y <- rbind(y, sapply(y[i - 1,], evolve))
rownames(y) <- paste("Sequence", 1:10)
# convert to AAbin object
rawaminos <- as.raw((65:89)[-c(2, 10, 15, 21, 24, 26)])
yAA <- apply(y, c(1, 2), function(s) rawaminos[match(s, aminos)])
class(yAA) <- "AAbin"

set.seed(999)
x.PHMM <- derivePHMM(x, residues = "DNA")
set.seed(999)
xDNA.PHMM <- derivePHMM(xDNA)
#plot(x.PHMM, from = 0, to = 10)

set.seed(999)
y.PHMM <- derivePHMM(y, residues = "AMINO", k = 2)
set.seed(999)
yAA.PHMM <- derivePHMM(yAA, k = 2)
# plot(y.PHMM, from = 0, to = 10)

test_that("objects have correct classes", {
  expect_is(x.PHMM, "PHMM")
  expect_is(xDNA.PHMM, "PHMM")
  expect_is(y.PHMM, "PHMM")
  expect_is(yAA.PHMM, "PHMM")
})
