### to do list ####
# multi sequence training for Baum Welch
# backward.PHMM function (and posterior)
# fix print.forward, etc
# read and write HMMs and PHMMs
# sequence simulation from phmm
# build trees
# multiple algnments (finish by aligning seqs to model)
# option to output derive in logspace
# fix includebegin arg in plot.HMM
# BLOSUM and PAM matrices
# integrate DNAbin and AAbin compatibility
# Viterbi output key as matrix
# end gaps in derive
# speed up Viterbi with Rcpp
# help documentation
# write vignette
# write paper

### less urgent ####
# bayesian derive methods
# insert-state allocation (see pg 122)


##### pairwise alignment example ####

# method = 'background'
# itertab = NULL
# offset = -0.1
# type = 'global'
# logspace = FALSE
# S = NULL
# qe = NULL
# residues = 'aminos'
# d = 8
# e = 2

x <- unlist(strsplit("VGATEVENVDHIKAEV", split = ""))
y <- unlist(strsplit("VGASTEVENVAEV", split = ""))
a <- align(x, y, S = BLOSUM50)
a
z <- unlist(strsplit("VGASTEVEVDHIKAEV", split = ""))
b <- align(z, a, S = BLOSUM50, residues = "aminos")
b

y <- unlist(strsplit("VGAAAASTEVENVDIKAEV", split = ""))
x <- align(x, y, S = BLOSUM50, residues = "aminos")
x
y <- unlist(strsplit("VGAAAASTEVENVDIKAEV", split = ""))

x <- align(x, y, S = BLOSUM50, residues = "aminos", type = 'global')
x

y <- unlist(strsplit("KPSVGAAADSTEVENVDIKAEV", split = ""))
x <- align(x, y, S = BLOSUM50, residues = "aminos", type = 'semiglobal')
x

hello <- align(x, y, S = BLOSUM50, residues = "aminos", type = 'global')
hello <- derive(hello, residues = 'aminos')
plot(hello, textexp = 0.5)
y <- x[,-(c(1, 5, 6, 7))]
x
align(x, y)

##############


# hello <- x
# for(i in 1:100) hello <- cbind(hello, x)
# system.time(derive(hello))
# Rprof(tmp <- tempfile())
# hello2 <- derive(hello)
# Rprof()
# summaryRprof(tmp)

a <- unlist(strsplit("VGATEVENVTAEV", split = ""))
b <- unlist(strsplit("VGATEVENVSAEV", split = ""))
a <- align(a, b, S = BLOSUM50)
a
b <- unlist(strsplit("VGATEENVPAEV", split = ""))
a <- align(a, b, S = BLOSUM50, residues = "aminos")
a

d <- unlist(strsplit("VGPTEENVTAEV", split = ""))
e <- unlist(strsplit("VGPTEENVSAEV", split = ""))
d <- align(d, e, S = BLOSUM50)
d
e <- unlist(strsplit("VGPTEENVPAEV", split = ""))
d <- align(d, e, S = BLOSUM50, residues = "aminos")
d
align(a, d)







x <- unlist(strsplit("VGATEVENVDHIKAEV", split = ""))
y <- unlist(strsplit("KPSVGAAADSTEVENVDIKAEV", split = ""))
hello <- Viterbi(x, y, type = 'global')
hello$path
hello$score
hello$progression
hello


hello <- Viterbi(x, y, type = 'semiglobal', offset = 0.2)
hello$path
hello$score
hello$progression

hello <- Viterbi(x, y, type = 'semiglobal', d = 2, e = 1)
hello$path
hello$score
hello$progression

hello <- Viterbi(x, y, type = 'local')
hello$path
hello$score
hello$progression


##### small alignment example ######

HBA_HUMAN  <- unlist(strsplit("VGA--HAGEY", split = ""))
HBB_HUMAN  <- unlist(strsplit("V----NVDEV", split = ""))
MYG_PHYCA  <- unlist(strsplit("VEA--DVAGH", split = ""))
GLB3_CHITP <- unlist(strsplit("VKG------D", split = ""))
GLB5_PETMA <- unlist(strsplit("VYS--TYETS", split = ""))
LGB2_LUPLU <- unlist(strsplit("FNA--NIPKH", split = ""))
GLB1_GLYDI <- unlist(strsplit("IAGADNGAGV", split = ""))

alig <- rbind(HBA_HUMAN, HBB_HUMAN, MYG_PHYCA, GLB3_CHITP,
              GLB5_PETMA, LGB2_LUPLU, GLB1_GLYDI)

x <- derive(alig)
x
plot(x)

y <- unlist(strsplit("VGASTEVENVEV", split = ""))
##  true alignment "VGA--NVAEV"
hello <- Viterbi(x, y)
hello$path
hello$progression
hello <- Viterbi(x, y, offset = 1)
hello$path








#### aligning two PHMMS

x1 <- unlist(strsplit("VGA--NAGEH", split = ""))
x2 <- unlist(strsplit("VGA--NVDEH", split = ""))
x3 <- unlist(strsplit("VGA--NVAEH", split = ""))
x4 <- unlist(strsplit("VGG---VAEH", split = ""))
x5 <- unlist(strsplit("VGA--NVEEH", split = ""))
x6 <- unlist(strsplit("FGA--NVPEH", split = ""))
x7 <- unlist(strsplit("VNADYNVAEH", split = ""))

x <- rbind(x1, x2, x3, x4, x5, x6, x7)
x <- derive(x, method = 'background', residues = 'aminos')
x

y1 <- unlist(strsplit("VGA--AGEH", split = ""))
y2 <- unlist(strsplit("VGA--VDEH", split = ""))
y3 <- unlist(strsplit("VGA--VAEH", split = ""))
y4 <- unlist(strsplit("VGG--VAEH", split = ""))
y5 <- unlist(strsplit("VGA--VEEH", split = ""))
y6 <- unlist(strsplit("FGA--SPEH", split = ""))
y7 <- unlist(strsplit("VNADYVAEH", split = ""))


y <- rbind(y1, y2, y3, y4, y5, y6, y7)
#align(x, y, S = BLOSUM50, residues = "aminos", type = 'global')
y <- derive(y, method = 'background', residues = 'aminos')
y

hello <- Viterbi(x, y, type = 'global')
hello$path
hello$score
hello$progression


hello <- Viterbi(x, y, type = 'semiglobal', offset = 0.5)
hello$path
hello$score
hello$progression


hello <- Viterbi(x, y, type = 'local')
hello$path
hello$score
hello$progression


plot(x)
plot(y)

qe <- rep(0.05, 20)
names(qe) <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
hello <- Viterbi(x, y, offset = -0.1, qe = qe)
hello
hello$path
hello$score
hello$progression





#### dishonest casino example ####
s <- c(0.99, 0.01)
names(s) <- c('Fair', 'Loaded')
A <- matrix(c(0.95, 0.1, 0.05, 0.9), nrow = 2)
dimnames(A) <- list(from = c('Fair', 'Loaded'),
                              to = c('Fair', 'Loaded'))
E <- matrix(c((1/6), (1/6), (1/6), (1/6), (1/6), (1/6),
                      (1/10),(1/10),(1/10),(1/10),(1/10),(1/2)),
                    nrow = 2, byrow = TRUE)
dimnames(E) <- list(states = c('Fair', 'Loaded'),
                            residues = paste(1:6))

#x <- buildHMM(s = s, A = A, E = E)
x <- structure(list(s = s, A = A, E = E), class = "HMM")

y = paste(c(3,1,5,1,1,6,2,4,6,4,4,6,6,4,4,2,4,5,3,1,1,3,2,1,6,3,1,1,6,4,1,5,2,1,3,3,
              6,2,5,1,4,4,5,4,3,6,3,1,6,5,6,6,2,6,5,6,6,6,6,6,6,5,1,1,6,6,4,5,3,1,3,2,
              6,5,1,2,4,5,6,3,6,6,6,4,6,3,1,6,3,6,6,6,3,1,6,2,3,2,6,4,5,5,2,3,6,2,6,6,
              6,6,6,6,2,5,1,5,1,6,3,1,2,2,2,5,5,5,4,4,1,6,6,6,5,6,6,5,6,3,5,6,4,3,2,4,
              3,6,4,1,3,1,5,1,3,4,6,5,1,4,6,3,5,3,4,1,1,1,2,6,4,1,4,6,2,6,2,5,3,3,5,6,
              3,6,6,1,6,3,6,6,6,4,6,6,2,3,2,5,3,4,4,1,3,6,6,1,6,6,1,1,6,3,2,5,2,5,6,2,
              4,6,2,2,5,5,2,6,5,2,5,2,2,6,6,4,3,5,3,5,3,3,3,6,2,3,3,1,2,1,6,2,5,3,6,4,
              4,1,4,4,3,2,3,3,5,1,6,3,2,4,3,6,3,3,6,6,5,5,6,2,5,6,6,6,6,2,6,3,2,6,6,6,
              6,1,2,3,5,5,2,4,5,2,4,2))


hello <- Viterbi(x, y)
hello$path
forward(x, obs)
backward(x, obs)
posterior(x, obs)
hello <- BaumWelch(x, obs)
hello
plot(hello, includebegin = T)
plot(x)
box()


#### multistate HMM ####
states <- c('A', 'B', 'C', 'D', 'E')
s <- c(0.5, 0.1, 0.1, 0.1, 0.1)
names(s) <- states

set.seed(123)
A <- matrix(rnorm(25, 100, 20), nrow = 5)
A <- A/apply(A, 1, sum)
dimnames(A) <- list(from = states, to = states)

E <- matrix(c((1/6), (1/6), (1/6), (1/6), (1/6), (1/6),
                      (1/10),(1/10),(1/10),(1/10),(1/10),(1/2),
                      (1/6), (1/6), (1/6), (1/6), (1/6), (1/6),
                      (1/6), (1/6), (1/6), (1/6), (1/6), (1/6),
                      (1/2),(1/10),(1/10),(1/10),(1/10),(1/10)),
                    nrow = 5, byrow = TRUE)
dimnames(E) <- list(states = states,
                            residues = paste(1:6))

x <- buildHMM(s = s, A = A, E = E)
plot(x, textexp = 0.9, includebegin = F)
box()






## Clustal omega multiple sequence alignmment
x <- list(
  A1 = unlist(strsplit('tgcggtatctacctacgacggtcgta', split = "")),
  A2 = unlist(strsplit('tgcggtatttacctacgacggtcgta', split = "")),
  A3 = unlist(strsplit('tggggtatctacctacgacggtcgta', split = "")),
  A4 = unlist(strsplit('tgcggtatctacccacgaccgtcgta', split = "")),
  A5 = unlist(strsplit('tgcggtatctacctacgacggttgta', split = "")),
  A6 = unlist(strsplit('tgcgccatctacctacgacggtcgta', split = "")),
  A7 = unlist(strsplit('agcggtatctacctacgacggtcgta', split = "")),
  A8 = unlist(strsplit('tgcggtatctaccaacgacggtcgta', split = "")),
  A9 = unlist(strsplit('tgcggtatctactacgacggtcgta', split = "")),
  A10 = unlist(strsplit('tgcggtatctaccgacggtcgta', split = "")),
  A11 = unlist(strsplit('tgcggtactacctacgacggcgta', split = "")),
  A12 = unlist(strsplit('tgcggtctacctaccggtcgta', split = "")),
  A13 = unlist(strsplit('tgcggtatctttacctacgacggtcgta', split = "")),
  A14 = unlist(strsplit('tgcggtatcttacctacgacggtcgta', split = "")),
  A15 = unlist(strsplit('tgcggtatctttgacctacgacggtcgta', split = ""))
)
res <- clustal(x)
n <- nrow(res)
hello3 <- matrix(nrow = n, ncol = n)
dimnames(hello3) <- list(rownames(res), rownames(res))
for(i in 1:n){
  for(j in 1:n){
    hello3[i, j] <- JC69(res[c(i, j),])$K
  }
}
newtree <- as.dendrogram(hclust(dist(hello3), method = "average"))
plot(newtree)

### some real data - coral ITS2 barcodes
# x <- c(paste("AY320", 289:352, sep=""),
#        paste("AY322", 575:612, sep=""),
#        paste("AY4580", 21:63, sep=""))
# library(ape)
# host.seq<-read.GenBank(x)
# write.dna(host.seq,"/home/shaun/Desktop/host.fas", format ='fasta')
# corals <- read.dna("/home/shaun/Desktop/host.fas", format ='fasta')
# corals <- lapply(corals, ape::as.character.DNAbin)
# save(corals, file = "data/corals.RData")


x <- lapply(corals, function(y) y[1:50]) #truncate to speed up comp
x <- x[1:20]
hello <- profile::clustal(x)
hello


# Rprof(tmp <- tempfile())
# hello <- clustalo(x[1:10])
# Rprof()
# summaryRprof(tmp)

#a0 m a or c
#c0 r a or g
#90 w a or t
#d0 d a g or t
#...

myseq <- randomSequence(2000)
mynewick <- "(A:1,(B:0.5,C:0.2):0.6);"
mydendrogram <- read.dendrogram(text = mynewick)
plot(mydendrogram)
mysimulation <- treesim(myseq, K = 0.01, n = 100, d = mydendrogram)
x <- treesim.as.matrix(mysimulation)
x <- rbind(x, x, x, x, x, x)
system.time(derive(x, residues = 'bases', method = 'background'))

