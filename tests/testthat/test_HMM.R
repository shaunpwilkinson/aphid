library(aphid)
context("build, train and apply standard HMMs")

### Dishonest casino example
states <- c("Begin", "Fair", "Loaded")
residues <- paste(1:6)
### Define transition probability matrix A
A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9), nrow = 3)
dimnames(A) <- list(from = states, to = states)
### Define emission probability matrix E
E <- matrix(c(rep(1/6, 6), rep(1/10, 5), 1/2), nrow = 2, byrow = TRUE)
dimnames(E) <- list(states = states[-1], residues = residues)
### Create the HMM object
x <- structure(list(A = A, E = E), class = "HMM")
### Load dice-roll data
data(casino)
actual <- c("F", "L")[match(names(casino), c("Fair", "Loaded"))]
### Actual path is stored in the names attribute of the sequence
### Calculate optimal path of hidden states
vit1 <- Viterbi(x, casino)
predicted <- c("F", "L")[vit1$path + 1]
### Implement forward and backward algorithms
for1 <- forward(x, casino)
bck1 <- backward(x, casino)

test_that("Dynammic programming found correct path and probs", {
  expect_equal(sum(predicted == "F"), 216)
  expect_equal(sum(predicted == "L"), 84)
  expect_equal(round(vit1$score, 2), -538.81)
  expect_equal(round(for1$score, 2), -516.45)
  expect_equal(for1$score, bck1$score)
})
