% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.PHMM}
\alias{plot.PHMM}
\title{Plot profile hidden Markov models.}
\usage{
\method{plot}{PHMM}(x, from = "start", to = "end", just = "center",
  arrexp = 1, textexp = 1, ...)
}
\arguments{
\item{x}{an object of class \code{"PHMM"}}

\item{from}{an integer giving the module number from where to begin the plot sequence.
Defaults to "start" which is converted to 0. Only applicable when plotting \code{"PHMM"} objects.}

\item{to}{an integer giving the module number at which to to finish the plot sequence.
Defaults to "start" which is converted to the total number of internal modules in the
model plus two (to account for the begin and end states).
Only applicable when plotting \code{"PHMM"} objects.}

\item{just}{a character string giving the justfication of the plot relative to the device.
Accepted values are "left", "center" and "right".}

\item{arrexp}{the expansion factor to be applied to the arrows in the plot.}

\item{textexp}{the expansion factor to be applied to the text in the plot.}
}
\description{
\code{plot.PHMM} provides a visual representation of a profile hidden Markov model.
}
\details{
Since the plotted models are generally much longer than they are
high, it is usually recommended to output the plot to a PDF file as demonstrated
in the example below.
}
\examples{
library(ape)
data(woodmouse)
woodmouse.PHMM <- derive.PHMM(woodmouse)

## plot partial model to viewer device
plot(woodmouse.PHMM, from = 0, to = 5)

## plot the entire model to a PDF
nr <- ceiling((woodmouse.PHMM$size + 2)/10)
pdf(file = "woodmouse.pdf",
   width = 8.27, height = nr * 2)
par(mfrow = c(nr, 1), mar = c(0, 0, 0, 0) + 0.1)
from <- 0
to <- 10
for(i in 1:nr){
  plot(woodmouse.PHMM, from = from, to = to, just = "left")
  from <- from + 10
  to <- min(to + 10, woodmouse.PHMM$size + 1)
}
dev.off()

}
