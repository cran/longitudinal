\name{tcell}
\alias{tcell}
\alias{tcell.10}
\alias{tcell.34}
\alias{tcell.gene.description}


\title{Microarray Time Series Data for T-Cell Activation }

\description{
  The data result from two experiments investigating
  the expression response of human T cells to PMA and ionomicin treatment.
 
  The first data set (\code{tcell.34}) contains the temporal expression levels of 58 genes
  for 10 unequally spaced time points.  At each time point there are 34 
  separate measurements. The second data set (\code{tcell.10}) stems from
  a related experiment considering 
  the same genes and identical time points, and contains 
  10 further measurements per time point. See Rangel et al. (2004) for more details.
}

\usage{
data(tcell)
}

\format{
  \code{tcell.10} and \code{tcell.34} are \code{\link{longitudinal}} objects,
  i.e. matrices with 58 colums each and a number of extra attributes
  (see \code{\link{longitudinal}} and \code{\link{longitudinal.util}}).
  
  The vector \code{tcell.gene.descriptions} contains the 
  description of the functions of the 58 investigated genes.
}

\source{
  This data is described in Rangel et al. (2004).
}


\references{
  Rangel, C., Angus, J., Ghahramani, Z., Lioumi, M., Sotheran, E., Gaiba, A.,
  Wild, D. L., and Falciani, F.
  (2004) Modeling T-cell activation using gene expression profiling and
  state-space models. \emph{Bioinformatics}, \bold{20}, 1361--1372.
}

\examples{
# load "longitudinal" library
library("longitudinal")

# load data sets
data(tcell)

# data set with 10 repeats 
dim(tcell.10)
summary(tcell.10)
is.longitudinal(tcell.10)
is.regularly.sampled(tcell.10)
is.equally.spaced(tcell.10)
get.time.repeats(tcell.10)

# data set with 34 repeats 
dim(tcell.34)
summary(tcell.34)
is.longitudinal(tcell.34)
is.regularly.sampled(tcell.34)
is.equally.spaced(tcell.34)
get.time.repeats(tcell.34)

# descriptions of the first nine genes
tcell.gene.description[1:9]

# plot the first nine time series
plot(tcell.10, 1:9)
plot(tcell.34, 1:9)

# Rangel et al. use the combined data set
tcell.44 <- combine.longitudinal(tcell.34, tcell.10)
plot(tcell.44, 1:9)
}
\keyword{datasets}

