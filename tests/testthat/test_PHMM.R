library(aphid)
context("build, train and apply profile HMMs")

# simulate a DNA sequence dataset
set.seed(999)
bases <- c("A", "C", "G", "T")
x <- list(sample(bases, replace = TRUE, size = 100))
evolve <- function(a) if(runif(1) > 0.95) sample(bases, 1) else a
for(i in 2:5) x[[i]] <- unname(sapply(x[[i - 1]], evolve))
names(x) <- paste("Sequence", 1:5)
# convert to DNAbin object
rawbases <- as.raw(c(136, 40, 72, 24))
xDNA <- lapply(x, function(s) rawbases[match(s, bases)])
class(xDNA) <- "DNAbin"

# simulate an AA sequence dataset, this time an alignment
set.seed(999)
aminos <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
y <- matrix(sample(aminos, replace = TRUE, size = 100), nrow = 1)
evolve <- function(a) if(runif(1) > 0.95) sample(aminos, 1) else a
for(i in 2:5) y <- rbind(y, sapply(y[i - 1,], evolve))
rownames(y) <- paste("Sequence", 1:5)
# convert to AAbin object
rawaminos <- as.raw((65:89)[-c(2, 10, 15, 21, 24, 26)])
yAA <- apply(y, c(1, 2), function(s) rawaminos[match(s, aminos)])
class(yAA) <- "AAbin"

# derive PHMMs for DNA sequences
set.seed(999)
x.PHMM <- derivePHMM(x, residues = "DNA", quiet = TRUE)
set.seed(999)
xDNA.PHMM <- derivePHMM(xDNA, quiet = TRUE)
#plot(x.PHMM, from = 0, to = 10)

# derive PHMMs for AA sequences
set.seed(999)
y.PHMM <- derivePHMM(y, residues = "AMINO", k = 2)
set.seed(999)
yAA.PHMM <- derivePHMM(yAA, k = 2)
# plot(y.PHMM, from = 0, to = 10)

# test forward and backward and Viterbi algorithms
x.for <- forward(x.PHMM, x[[1]], cpp = FALSE, residues = "DNA")
xDNA.for <- forward(xDNA.PHMM, xDNA[[1]])
x.bck <- backward(x.PHMM, x[[1]], cpp = FALSE, residues = "DNA")
xDNA.bck <- backward(xDNA.PHMM, xDNA[[1]])
x.vit <- Viterbi(x.PHMM, x[[1]], cpp = FALSE, residues = "DNA", S = substitution$NUC.4.4)
xDNA.vit <- Viterbi(xDNA.PHMM, xDNA[[1]], S = substitution$NUC.4.4)

y.for <- forward(y.PHMM, y[1, ], cpp = FALSE, residues = "AMINO")
yAA.for <- forward(yAA.PHMM, yAA[1, ])
y.bck <- backward(y.PHMM, y[1, ], cpp = FALSE, residues = "AMINO")
yAA.bck <- backward(yAA.PHMM, yAA[1, ])
y.vit <- Viterbi(y.PHMM, y[1, ], cpp = FALSE, residues = "AMINO", S = substitution$MATCH)
yAA.vit <- Viterbi(yAA.PHMM, yAA[1, ], S = substitution$MATCH)

# Ensure cpp and R evaluations give same results
x2.vit <- Viterbi(x.PHMM, x[[1]], cpp = FALSE, type = "semiglobal")
xDNA2.vit <- Viterbi(xDNA.PHMM, xDNA[[1]], type = "semiglobal")
x3.vit <- Viterbi(x.PHMM, x[[1]], cpp = FALSE, type = "local")
xDNA3.vit <- Viterbi(xDNA.PHMM, xDNA[[1]], type = "local")


# Generate random sequences
x.sim <-generate(x.PHMM, size = 200)
xDNA.sim <-generate(xDNA.PHMM, size = 200, DNA = TRUE)
y.sim <-generate(y.PHMM, size = 200)
yAA.sim <-generate(yAA.PHMM, size = 200, AA = TRUE)

# Baum Welch training
xbw.PHMM <- train(x.PHMM, c(x, x.sim), method = "BaumWelch", deltaLL = 0.1, quiet = TRUE)
xDNAbw.PHMM <- train(xDNA.PHMM, c(xDNA, xDNA.sim), method = "BaumWelch", deltaLL = 0.1, quiet = TRUE)
ybw.PHMM <- train(y.PHMM, c(y, y.sim), method = "BaumWelch", deltaLL = 0.1, quiet = TRUE)
yAAbw.PHMM <- train(yAA.PHMM, c(yAA, yAA.sim), method = "BaumWelch", deltaLL = 0.1, quiet = TRUE)

# Read and write HMMER files
fl <- tempfile()
x.HMMER <- writePHMM(x.PHMM, file = fl)
x2.PHMM <- readPHMM(fl)
fl <- tempfile()
y.HMMER <- writePHMM(y.PHMM, file = fl)
y2.PHMM <- readPHMM(fl)

# Multiple sequence alignments
set.seed(999)
x.alig <- align(x, quiet = TRUE)
set.seed(999)
x2.alig <- align(x, progressive = TRUE, quiet = TRUE)



test_that("objects have correct classes", {
  expect_is(x.PHMM, "PHMM")
  expect_is(xDNA.PHMM, "PHMM")
  expect_is(y.PHMM, "PHMM")
  expect_is(yAA.PHMM, "PHMM")
  expect_is(x2.PHMM, "PHMM")
  expect_is(y2.PHMM, "PHMM")
  expect_is(x.for, "fullprob")
  expect_is(x.bck, "fullprob")
  expect_is(x.vit, "Viterbi")
})


test_that("Character and raw inputs give same output", {
  expect_equal(round(unname(x.PHMM$A), 2), round(unname(xDNA.PHMM$A), 2))
  expect_equal(round(unname(x.PHMM$E), 2), round(unname(xDNA.PHMM$E), 2))
  expect_equal(round(unname(y.PHMM$A), 2), round(unname(yAA.PHMM$A), 2))
  expect_equal(round(unname(y.PHMM$E), 2), round(unname(yAA.PHMM$E), 2))
  expect_equal(round(unname(y.PHMM$E), 2), round(unname(y2.PHMM$E), 2))
  expect_equal(x.PHMM$size, xDNA.PHMM$size)
  expect_equal(x.PHMM$size, x2.PHMM$size)
  expect_equal(y.PHMM$size, yAA.PHMM$size)
  expect_equal(y.PHMM$size, y2.PHMM$size)
  expect_equal(x.alig, x2.alig)
  expect_equal(x.for$score, x.bck$score)
  expect_equal(xDNA.for$score, xDNA.bck$score)
  expect_equal(y.for$score, y.bck$score)
  expect_equal(yAA.for$score, yAA.bck$score)
  expect_equal(round(y.for$score, 2), round(yAA.bck$score, 2))
  expect_equal(x.vit$path, xDNA.vit$path)
  expect_equal(round(x.vit$score), round(xDNA.vit$score))
  expect_equal(y.vit$path, yAA.vit$path)
  expect_equal(y.vit$score, yAA.vit$score)
})

test_that("Plotting command functions as expected", {
  expect_identical(plot(x.PHMM, from = 0, to = 10), NULL)
})
