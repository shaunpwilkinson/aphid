# aphid

Analysis with profile hidden Markov models

[![DOI](https://zenodo.org/badge/63536088.svg)](https://zenodo.org/badge/latestdoi/63536088)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/aphid)](https://cran.r-project.org/package=aphid)
[![](http://cranlogs.r-pkg.org/badges/grand-total/aphid)](https://cran.r-project.org/package=aphid)
[![Build Status](https://travis-ci.org/shaunpwilkinson/aphid.svg?branch=master)](https://travis-ci.org/shaunpwilkinson/aphid)
[![codecov](https://codecov.io/github/shaunpwilkinson/aphid/branch/master/graphs/badge.svg)](https://codecov.io/github/shaunpwilkinson/aphid)
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--7332--7931-brightgreen.svg)](http://orcid.org/0000-0002-7332-7931)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)


--------------------------------------------------------------------------------

`aphid` is an R package for the development and application of
hidden Markov models and profile HMMs for biological sequence analysis.
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
To download `aphid` from CRAN and load the package, run

```R
install.packages("aphid")
library("aphid")
```

To download the development version from 
GitHub, first ensure a C/C++ compliler is available and the 
[devtools](https://github.com/hadley/devtools) R package is installed. 
Linux users will generally have a compiler installed by default; 
however Windows users may need to download 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac 
OSX users will need Xcode (note that these are not R packages). 
Install and load the package by running 

```R
devtools::install_github("shaunpwilkinson/aphid", build_vignettes = TRUE) 
library("aphid")
```

### Use and Examples
An overview of the package and its functions can be found by running

```R
?aphid
```

To view the tutorial, run

```R
vignette("aphid-vignette")
```

### Issues
If you experience a problem using this package please feel free to
raise it as an issue on [GitHub](http://github.com/shaunpwilkinson/aphid/issues).


### Acknowledgements
This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.

