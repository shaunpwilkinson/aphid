## ---- echo = FALSE, message = FALSE, warning = FALSE---------------------
#knitr::opts_chunk$set(out.width='750px', dpi=200)
#knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- fig.width=7.15, fig.height=4---------------------------------------
library("aphid")
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
### Plot the model
plot(x, textexp = 1.5)
### Optionally add the transition probabilities as text
text(x = 0.02, y = 0.5, labels = "0.95")
text(x = 0.51, y = 0.5, labels = "0.90")
text(x = 0.5, y = 0.9, labels = "0.05")
text(x = 0.5, y = 0.1, labels = "0.10")

## ---- echo = FALSE-------------------------------------------------------
data(casino)
cat("", 
    paste0(casino[1:50], collapse = ""), "\n", 
    paste0(casino[51:100], collapse = ""), "\n", 
    paste0(casino[101:150], collapse = ""), "\n", 
    paste0(casino[151:200], collapse = ""), "\n", 
    paste0(casino[201:250], collapse = ""), "\n", 
    paste0(casino[251:300], collapse = ""), "\n")

## ------------------------------------------------------------------------
data(casino)
### The actual path is stored in the names attribute of the sequence
actual <- c("F", "L")[match(names(casino), c("Fair", "Loaded"))]
### Find the predicted path
vit1 <- Viterbi(x, casino)
predicted <- c("F", "L")[vit1$path + 1]
### Note the path element of the output Viterbi object is an integer vector
### the addition of 1 to the path converts from C/C++ to R's indexing style

## ---- echo = FALSE-------------------------------------------------------

cat("", 
    "Actual    ", paste0(actual[1:50], collapse = ""), "\n", 
    "Predicted ", paste0(predicted[1:50], collapse = ""), "\n\n", 
    "Actual    ", paste0(actual[51:100], collapse = ""), "\n",
    "Predicted ", paste0(predicted[51:100], collapse = ""), "\n\n", 
    "Actual    ", paste0(actual[101:150], collapse = ""), "\n", 
    "Predicted ", paste0(predicted[101:150], collapse = ""), "\n\n", 
    "Actual    ", paste0(actual[151:200], collapse = ""), "\n", 
    "Predicted ", paste0(predicted[151:200], collapse = ""), "\n\n", 
    "Actual    ", paste0(actual[201:250], collapse = ""), "\n", 
    "Predicted ", paste0(predicted[201:250], collapse = ""), "\n\n", 
    "Actual    ", paste0(actual[251:300], collapse = ""), "\n",
    "Predicted ", paste0(predicted[251:300], collapse = ""), "\n")

## ---- fig.width=7.15, fig.height=4---------------------------------------
casino.post <- posterior(x, casino)
plot(1:300, seq(0, 1, length.out = 300), type = "n", xlab = "Roll number",
     ylab = "Posterior probability of dice being fair")
starts <- which(c("L", actual) == "F" & c(actual, "F") == "L")
ends <- which(c("F", actual) == "L" & c(actual, "L") == "F") - 1
for(i in 1:6) rect(starts[i], 0, ends[i], 1, col = "grey", border = NA)
lines(1:300, casino.post[1, ])

## ---- fig.width=7.15, fig.height=4---------------------------------------
y <- deriveHMM(list(casino), logspace = FALSE)
plot(y, textexp = 1.5)

### Optionally add the transition probabilities as text
text(x = 0.02, y = 0.5, labels = round(y$A["Fair", "Fair"], 2))
text(x = 0.51, y = 0.5, labels = round(y$A["Loaded", "Loaded"], 2))
text(x = 0.5, y = 0.9, labels = round(y$A["Fair", "Loaded"], 2))
text(x = 0.5, y = 0.1, labels = round(y$A["Loaded", "Fair"], 2))

## ------------------------------------------------------------------------
data(globins)
globins

## ---- fig.width=7.15, fig.height=5---------------------------------------
globins.PHMM <- derivePHMM(globins, residues = "AMINO", seqweights = NULL)
plot(globins.PHMM)

## ------------------------------------------------------------------------
path <- Viterbi(globins.PHMM, globins["GLB1_GLYDI", ])$path
path

## ------------------------------------------------------------------------
c("D", "M", "I")[path + 1]

## ------------------------------------------------------------------------
sim <- list(length = 10)
set.seed(9999)
for(i in 1:10) sim[[i]] <- generate(globins.PHMM, size = 20)
sim

## ------------------------------------------------------------------------
sim <- lapply(sim, function(s) s[names(s) != "D"])

## ------------------------------------------------------------------------
globins2.PHMM <- train(globins.PHMM, sim, method = "BaumWelch", 
                       deltaLL = 0.01, seqweights = NULL)

## ------------------------------------------------------------------------
globins <- unalign(globins)
align(globins, model = globins.PHMM, seqweights = NULL, residues = "AMINO")

