# aphid

[![Build Status](https://travis-ci.org/shaunpwilkinson/aphid.svg?branch=master)](https://travis-ci.org/shaunpwilkinson/aphid)
[![codecov](https://codecov.io/github/shaunpwilkinson/aphid/branch/master/graphs/badge.svg)](https://codecov.io/github/shaunpwilkinson/aphid)
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--7332--7931-brightgreen.svg)](http://orcid.org/0000-0002-7332-7931) 
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

Analysis with profile hidden Markov models in R

--------------------------------------------------------------------------------

`aphid` is an R package containing functions for building and using 
profile hidden Markov models for biological sequence analysis. 
Functions are included for multiple and pairwise sequence alignment, 
model construction and parameter optimization, calculation of conditional 
probabilities (using the forward, backward and Viterbi algorithms),
tree-based sequence weighting, sequence simulation, and file import/export 
compatible with the [HMMER](http://www.hmmer.org/) software package. 
`aphid` also includes functions for developing and working with 
standard hidden Markov models.

This package was written based on the algorithms described in the book 
[Biological Sequence Analysis](
https://www.amazon.com/Biological-Sequence-Analysis-Probabilistic-Proteins/dp/0521629713)
by Richard Durbin, Sean Eddy, Anders Krogh and Graeme Mitchison. 
This book offers an in depth explanation of hidden Markov models and 
profile HMMs for users of all levels of familiarity. 
Many of the examples and datasets in the package are directly derived from the 
text, which serves as a useful primer for this package.

### Installation
To download the development version of the package from 
GitHub you will first need to ensure you have a C/C++ compliler and the 
[devtools](https://github.com/hadley/devtools) R package installed. 
Linux users will generally have a compiler such as `gcc` installed by default; 
however Windows users will need to download 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac 
OSX users will need [Xcode](https://developer.apple.com/xcode) 
(note that Rtools and Xcode are not R packages). To download and install 
devtools, run 
```R
install.packages("devtools")
``` 
then install and load `aphid` by running 
```R
devtools::install_github("shaunpwilkinson/aphid") 
library("aphid")
```

## Use and Examples

### Help
An overview of the package and it's functions can be found by running
```R
?aphid
```
To build the vignette users will need to have LaTeX installed. RStudio recommends 
[MiKTeX Complete](http://miktex.org/2.9/setup) for Windows and
[TexLive 2013 Full](http://tug.org/) for Mac OS X and Linux.

If you experience a problem using this package please feel free to
raise it as an issue on [GitHub](http://github.com/shaunpwilkinson/aphid/issues).


### Acknowledgements
This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.

