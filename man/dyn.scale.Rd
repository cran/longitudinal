\name{dyn.scale}
\alias{dyn.scale}
\alias{dyn.moments}
\alias{dyn.weights}
\alias{time2weights}


\title{Dynamical Scale, Moments, and Weights}

\description{
  \code{time2weights} computes weights corresponding to time points

  \code{dyn.weights} computes these weights for a given \code{\link{longitudinal}} matrix.
  
  \code{dyn.moments} computes means and variances for the variables in
        a \code{\link{longitudinal}} object.
  
  \code{dyn.scale} centers and standardizes a \code{\link{longitudinal}} matrix.
  
 }

\usage{
  time2weights(t)
  dyn.weights(x)
  dyn.moments(x)
  dyn.scale(x, center=TRUE, scale=TRUE)
}

\arguments{
  \item{t}{a vector with time points}
  \item{x}{a \code{\link{longitudinal}} object, or a matrix}
  \item{center}{logical value}
  \item{scale}{logical value}
}


\value{
  A vector with weights (\code{time2weights} and \code{dyn.weights}), 
  a list containing the column means and variances (\code{dyn.moments}), 
  or a rescaled longitudinal matrix (\code{dyn.scale}).
}

\details{
  The dynamical weights are computed assuming a linear spline - see Opgen-Rhein
  and Strimmer (2006a,b).  The dynamical mean and variance etc. are then simply 
  weighted versions of the usual empirical estimators.
} 


\author{
  Rainer Opgen-Rhein and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}

\references{

 Opgen-Rhein, R., and K. Strimmer. 2006a. Inferring gene dependency networks from 
 genomic longitudinal data: a functional data approach.
 REVSTAT \bold{4}:53-65. (\url{http://www.ine.pt/revstat/pdf/rs060103.pdf})  

 Opgen-Rhein, R., and K. Strimmer. 2006b. Using regularized dynamic correlation 
 to infer gene dependency networks from time-series microarray data. 
 The 4th International Workshop on Computational Systems Biology,
 WCSB 2006 (June 12-13, 2006, Tampere, Finland). 
 (\url{http://www.cs.tut.fi/wcsb06/})  
}


\seealso{\code{\link[corpcor]{wt.scale}}.}

\examples{
# load "longitudinal" library
library("longitudinal")

# weights of for the data points in tcell data
data(tcell)
dyn.weights(tcell.34)

# dynamical moments
dyn.moments(tcell.34)

}
\keyword{multivariate}
