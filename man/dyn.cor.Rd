\name{dyn.cor}
\alias{dyn.cor}
\alias{dyn.pcor}
\alias{dyn.invcor}
\alias{dyn.cov}
\alias{dyn.pvar}
\alias{dyn.var}
\alias{dyn.invcov}


\title{Dynamical Correlation and Covariance}
\description{
  The functions estimate dynamical correlation and covariance, and
  related quantities.  
}
\usage{
dyn.cor(x, lambda, verbose=TRUE)
dyn.var(x, lambda.var, verbose=TRUE)
dyn.cov(x, lambda, lambda.var, verbose=TRUE)
dyn.invcor(x, lambda, verbose=TRUE)
dyn.invcov(x, lambda, lambda.var, verbose=TRUE)
dyn.pvar(x, lambda, lambda.var, verbose=TRUE)
dyn.pcor(x, lambda, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the correlation shrinkage intensity (range 0-1). If \code{lambda} is not specified 
	 (the default) it is estimated using an analytic formula from Sch\"afer 
	 and Strimmer (2005) and Opgen-Rhein and Strimmer (2006a,b).
	  For \code{lambda=0} the empirical correlations are recovered.
	  See also \code{\link[corpcor]{cor.shrink}}.}
  \item{lambda.var}{the variance shrinkage intensity (range 0-1). If \code{lambda.var} is not specified 
	 (the default) it is estimated using an analytic formula from Sch\"afer 
	 and Strimmer (2005) and Opgen-Rhein and Strimmer (2006a,b).
	  For \code{lambda.var=0} the empirical variances are recovered.
	  See also \code{\link[corpcor]{var.shrink}}.}  
  \item{verbose}{report progress while computing (default: TRUE)}
}
\details{
   The dynamical correlation and related quantities implemented here
   follow the definition of Opgen-Rhein and Strimmer (2006a,b).
   This approach is derived from a FDA perspective. Essentially, it takes account 
   of the distances between the various time points by assigning weights to samples.
   If these weights are all equal the usual iid estimators are obtained.
  
   For details about the analytic shrinkage procedure consult 
   Opgen-Rhein and Strimmer (2006b) and  Sch\"afer and Strimmer (2005)
   as well as the help page of \code{\link[corpcor]{cov.shrink}}.
   
}
\value{

  \code{dyn.cor} returns the dynamical correlation matrix.
  \code{dyn.var} returns the vector of dynamical variances.
  \code{dyn.cov} returns the dynamical covariance matrix. 

  \code{dyn.invcor} returns the inverse dynamical correlation matrix.
  \code{dyn.invcov} returns the inverse dynamical covariance matrix. 

  \code{dyn.pvar} returns the vector of partial dynamical variances.
  \code{dyn.pcor} returns the partial dynamical correlation matrix.  
}
\author{
  Rainer Opgen-Rhein and
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\references{

 Opgen-Rhein, R., and K. Strimmer. 2006a. Inferring gene dependency networks from 
 genomic longitudinal data: a functional data approach.
 REVSTAT \bold{4}:53-65.

 Opgen-Rhein, R., and K. Strimmer. 2006b. Using regularized dynamic correlation 
 to infer gene dependency networks from time-series microarray data. 
 The 4th International Workshop on Computational Systems Biology,
 WCSB 2006 (June 12-13, 2006, Tampere, Finland). 
  
 Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
 covariance estimation and implications for functional genomics. 
 Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
 <DOI:10.2202/1544-6115.1175>


}
\seealso{\code{\link{dyn.weights}}, \code{\link[corpcor]{cov.shrink}}, \code{\link[corpcor]{pcor.shrink}}}

\examples{
# load "longitudinal" library
library("longitudinal")

# load tcell data
data(tcell)
get.time.repeats(tcell.34)


# dynamical partial correlation 
# (this takes into account of the unequal spacings between time points) 
dynpc <- dyn.pcor(tcell.34, lambda=0)

# static partial correlation
statpc <- pcor.shrink(tcell.34, lambda=0)

# this is NOT the same
sum((dynpc - statpc)^2)

}
\keyword{multivariate}
