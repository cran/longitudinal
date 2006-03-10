\name{longitudinal}
\alias{longitudinal}
\alias{as.longitudinal}
\alias{is.longitudinal}
\alias{summary.longitudinal}
\alias{print.longitudinal}
\alias{plot.longitudinal}
  
  
\title{Data Structure for Longitudinal Data}
\usage{
as.longitudinal(x, repeats=1, time)
is.longitudinal(x)
summary.longitudinal(object, \ldots)
print.longitudinal(x, \ldots)
plot.longitudinal(x, series=1, type=c("median", "mean"), autolayout=TRUE, \ldots)
}
\arguments{
  \item{x, object}{Time series data, contained in a \code{longitudinal} object
           or in matrix form (\code{as.longitudinal}).} 
	   
  
  \item{repeats}{Integer, or a vector of integers, that specifies the number
                 of available measurements per time point.  If only one number is given
		 then it is assumed the time series is regularily sampled. If instead a vector
		 is specified, then each time point may have a different number of samples.}  

  \item{time}{A vector with the dates for the time points. If not specified,
               equally spaced time points  1, 2, 3, \ldots are assumed. }	   

  
  \item{series}{Number, or a vector of numbers, that indicates which time series
           (=variables and columns of x) are plotted.}
  
  \item{type}{Determines whether the plotted line corresponds to the mean or median value
               of the samples per time point (default: "median").}
  
  \item{autolayout}{determine the number of
                    columns and rows in the plot automatically in the display of multiple time series  (default: TRUE). }	         	   
  
  \item{\ldots}{Additional optional parameters}
}
\description{

  The data type \code{longitudinal} stores multiple time series data.  It allows repeated
  measurements, irregular sampling, and unequal temporal spacing of the time points. 
  
  \code{as.longitudinal} converts a matrix into a \code{longitudinal} object.  The columns of the
  input matrix are considered as individual variables (time series).  The rows contain the
  measurements in temporal order (for instance, rows 1-10 could contain 10 repeated measurements taken at
  time point 1, rows 11-20 further 10 measurements taken at time point 2, and so on).
  The dates for the time points can be specified with the argument \code{times} and need not
  be equally spaced.  With the argument \code{repeats} it is possible to specify the number
  of measurements per time point (this may be different from time point to time point).  In the resulting
  \code{longitudinal} matrix object the rownames will indicate both the time points and the 
  repetition number (e.g., "10-1", "10-2", "10-3", \ldots,  "20-1", "20-2", "20-3", etc.).
  
  \code{is.longitudinal} checks whether a matrix has the \code{longitudinal} attributes.
  
  The functions \code{summary}, \code{print}, \code{plot} are the standard generic functions
  adapted to \code{longitudinal} objects.
}

\value{

  \code{as.longitudinal} returns a \code{longitudinal} object.
  
  \code{is.longitudinal} returns TRUE or false.
 
}


\author{
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}
\seealso{\code{\link{longitudinal.util}}, \code{\link{tcell}}, \code{\link{ts}}.}

\examples{
# load "longitudinal" library
library("longitudinal")

# load data set
data(tcell)
is.longitudinal(tcell.34)
attributes(tcell.34)
tcell.34[,1:3]

# how many samples and how many genes?
dim(tcell.34)
summary(tcell.34)

# plot first nine time series
plot(tcell.34, 1:9)

#####

# an artificial example with repeated measurements, irregular sampling, and unequal spacing 
m <- matrix(rnorm(200), 50, 4)
z <- as.longitudinal(m, repeats=c(10,5,5,10,20), time=c(2,8,9,15,16))
plot(z, 1:4)
}
\keyword{ts}
