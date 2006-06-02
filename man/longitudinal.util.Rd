\name{longitudinal.util}
\alias{longitudinal.util}

\alias{get.time.repeats}
\alias{is.equally.spaced}
\alias{is.regularly.sampled}
\alias{has.repeated.measurements}
\alias{combine.longitudinal}
\alias{condense.longitudinal}


\title{Utility Functions for the "Longitudinal" Data Structure}
\usage{
get.time.repeats(x)
is.equally.spaced(x)
is.regularly.sampled(x)
has.repeated.measurements(x)
combine.longitudinal(x1, x2)
condense.longitudinal(x, s, func=median)
}
\arguments{
  \item{x, x1, x2}{\code{\link{longitudinal}} time series objects}
  
  \item{s}{An integer, or a vector of integers, that designate the set of
           time series (variables) to condense.}
  \item{func}{Univariate function used to summarize the multiple measurements
              per time point.}

}
\description{

  The above functions are all utility functions for  \code{\link{longitudinal}} objects.
  
  \code{get.time.repeats} returns the measurement design, i.e. the time points
  and the number of repeats per time point. 
  
  \code{is.equally.spaced} checks whether the distances between subsequent time points
  are all equal.   
  
  \code{is.regularly.sampled} checks whether the number of measurements are identical
  across time points.   
 
  \code{has.repeated.measurements} checks whether any time point as been measured more
  than once.
  
  
  \code{combine.longitudinal} combines the measurements of two \code{longitudinal} objects. 
  These objects must have the same (number of) variables. 
  
  \code{condense.longitudinal} condenses the multiple measurements per time point using
  an arbitrary function (e.g., mean, median, var).
}

\value{

  \code{get.time.repeats} returns a list containing two vectors (\code{time} and \code{repeats}).
  
  \code{is.equally.spaced}, \code{is.regularly.sampled}, and \code{has.repeated.measurements}
  return either TRUE or FALSE.
  
  \code{combine.longitudinal} returns a \code{\link{longitudinal}} object.
  
  \code{condense.longitudinal} returns a matrix.
}


\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{longitudinal}}, \code{\link{tcell}}.}

\examples{
# load "longitudinal" library
library("longitudinal")


# load tcell data set
data(tcell)
dim(tcell.34)
is.longitudinal(tcell.34)
summary(tcell.34)

# information
get.time.repeats(tcell.34)
is.equally.spaced(tcell.34)
is.regularly.sampled(tcell.34)
has.repeated.measurements(tcell.34)

# compute the mean value at each time point for the first two gene
condense.longitudinal(tcell.34, 1:2, mean)


# combine two time series
m1 <- matrix(rnorm(100), 50, 2)
m2 <- matrix(rnorm(100), 50, 2)
z1 <- as.longitudinal(m1, repeats=c(10,5,5,10,20), time=c(2,8,9,15,16))
z2 <- as.longitudinal(m2, repeats=c(10,5,5,10,20), time=c(1,8,9,15,20))

z3 <- combine.longitudinal(z1,z2)
summary(z3)
get.time.repeats(z3)  # compare with z1 and z2
}
\keyword{ts}
