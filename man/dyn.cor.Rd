\name{dyn.cor}
\alias{dyn.cor}
\alias{dyn.pcor}
\alias{dyn.invcor}
\alias{dyn.cov}
\alias{dyn.pcov}
\alias{dyn.invcov}


\title{Dynamical Correlation and Covariance}
\description{
  The functions estimate the dynamical correlation and covariance, and
  related quantities.  
}
\usage{
dyn.cor(x, lambda, verbose=TRUE)
dyn.cov(x, lambda, verbose=TRUE)
dyn.invcor(x, lambda, verbose=TRUE)
dyn.invcov(x, lambda, verbose=TRUE)
dyn.pcor(x, lambda, verbose=TRUE)
dyn.pcov(x, lambda, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the shrinkage intensity (range 0-1) employed for computing the underlying
  correlation matrix. If \eqn{\lambda} is is not specified (the default)
        a suitable value is automatically chosen - see \code{\link[corpcor]{cor.shrink}}.  
	 For \eqn{\lambda=0} no shrinkage occurs
	 and the empirical estimates are computed.}
  \item{verbose}{report progress while computing (default: TRUE)}
}
\details{
   The dynamical correlation and related quantities implemented here
   follow the definition of Opgen-Rhein and K. Strimmer (2006a,b).
   This approach derives from a FDA perspective. Essentially, it takes account 
   of the distances between the various time points by introducing weights.
   If these weights are all equal the usual estimators are obtained.
  
   For details about the analytic shrinkage procedure consult 
   Opgen-Rhein and K. Strimmer (2006b) and  Schaefer and Strimmer (2005)
   as well as the help page of \code{\link[corpcor]{cov.shrink}}.
   
}
\value{

  \code{dyn.cor} returns the dynamical correlation matrix.
  \code{dyn.cov} returns the dynamical covariance matrix. 

  \code{dyn.invcor} returns the inverse dynamical correlation matrix.
  \code{dyn.invcov} returns the inverse dynamical covariance matrix. 

  \code{dyn.pcor} returns the partial dynamical correlation matrix.
  \code{dyn.pcov} returns the partial dynamical covariance matrix.  
}
\author{
  Rainer Opgen-Rhein (\url{http://www.statistik.lmu.de/~opgen/}) and
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}
\references{

 Opgen-Rhein, R., and K. Strimmer. 2006a. Inferring gene dependency networks from 
 genomic longitudinal data: a functional data approach.
 REVSTAT \bold{4}:53-65. (\url{http://http://www.ine.pt/revstat/}) 

 Opgen-Rhein, R., and K. Strimmer. 2006b. Using regularized dynamic correlation 
 to infer gene dependency networks from time-series microarray data. 
 The 4th International Workshop on Computational Systems Biology,
 WCSB 2006 (June 12-13, 2006, Tampere, Finland). 
 (\url{http://www.cs.tut.fi/wcsb06/})  
  
 Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
 covariance estimation and implications for functional genomics. 
 Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
 (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})  
}
\seealso{\code{{dyn.weights}}, \code{\link[corpcor]{cov.shrink}}}

\examples{
# load "longitudinal" library
library("longitudinal")

# load tcell data
data(tcell)
get.time.repeats(tcell.34)


# dynamical partial correlation 
# (this takes into account the unequally distant time points) 
dynpc <- dyn.pcor(tcell.34, lambda=0)

# static partial correlation
statpc <- pcor.shrink(tcell.34, lambda=0)

# this is NOT the same
sum((dynpc - statpc)^2)

}
\keyword{multivariate}
