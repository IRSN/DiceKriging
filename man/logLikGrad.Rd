\name{logLikGrad}
\alias{logLikGrad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Concentrated log-Likelihood of a km object - Analytical gradient }
\description{
   Returns the exact gradient vector of the concentrated log-likelihood.
}

\usage{
logLikGrad(param, model, envir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{ a vector containing the optimization variables. }
  \item{model}{ an object of class \code{km}. }
  \item{envir}{ an environment specifying where to get intermediate values calculated in \code{logLikFun}. }
}

\value{
  the gradient of the concentrated log-likelihood.
  }

\references{

J.-S. Park and J. Baek (2001), Efficient computation of maximum likelihood estimators in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences}, \bold{27} no. 1, 1-7.
	
}

\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne}

\seealso{ \code{\link{km}}, \code{\link{logLikFun}} }

\keyword{models}
\keyword{htest}
\keyword{internal}
