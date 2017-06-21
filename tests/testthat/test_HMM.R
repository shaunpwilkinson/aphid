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
vit2 <- Viterbi(x, casino, cpp = FALSE)
predicted <- c("F", "L")[vit1$path + 1]

### Implement forward, backward algorithms
for1 <- forward(x, casino)
for2 <- forward(x, casino, cpp = FALSE)
bck1 <- backward(x, casino)
bck2 <- backward(x, casino, cpp = FALSE)
post1 <- posterior(x, casino)
post2 <- posterior(x, casino, cpp = FALSE)

### Simulate five random sequences
sim <- vector(mode = "list", length = 5)
for(i in 1:5) sim[[i]] <- generate(x, size = 300)

### Train the models with simulated data
x2 <- train(x, sim, method = "BaumWelch", deltaLL = 0.01, quiet = TRUE)
x3 <- train(x, sim, method = "Viterbi", quiet = TRUE)

### Plotting
fl <- tempfile(fileext=".pdf")
pdf(file = fl)
x.plot <- plot(x, begin = TRUE)
dev.off()


test_that("Dynammic programming found correct path and probs", {
  expect_equal(sum(predicted == "F"), 216)
  expect_equal(sum(predicted == "L"), 84)
  expect_equal(round(vit1$score, 2), -538.81)
  expect_equal(round(for1$score, 2), -516.45)
  expect_equal(for1$score, bck1$score)
  expect_equal(sum(apply(post1, 2, sum)), 300)
})

test_that("Cpp and R gave congruent results", {
  expect_equal(vit1$score, vit2$score)
  expect_equal(vit1$path, vit2$path)
  expect_equal(for1$score, for2$score)
  expect_equal(bck1$score, bck2$score)
  expect_equal(post1, post2)
})

test_that("Model training yielded acceptable results", {
  expect_is(x2, "HMM")
  expect_is(x3, "HMM")
})

test_that("Plotting command functions as expected", {
  expect_identical(x.plot, NULL)
})
