\name{logLikFun}
\alias{logLikFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Concentrated log-likelihood of a km object }
\description{
Returns the concentrated log-likelihood, obtained from the likelihood by plugging in the estimators of the parameters that can be expressed in function of the other ones.
}
\usage{
logLikFun(param, model, envir=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{ a vector containing the optimization variables. }
  \item{model}{ an object of class \code{km}. }
  \item{envir}{ an optional environment specifying where to assign intermediate values for future gradient calculations. Default is NULL.}
}

\value{
The concentrated log-likelihood value.
}

\details{
When there is no nugget effect nor observation noise, the concentrated log-likelihood is obtained by plugging in the variance and the trend MLE. Maximizing the likelihood is then equivalent to maximizing the concentrated log-likelihood with respect to the covariance parameters. In the other cases, the maximization of the concentrated log-likelihood also involves other parameters (the variance explained by the stationary part of the process for noisy observations, and this variance divided by the total variance if there is an unknown homogeneous nugget effect).
}

\references{ 

J.-S. Park and J. Baek (2001), Efficient computation of maximum likelihood estimators in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences}, \bold{27} no. 1, 1-7.
	
}

\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne}


\seealso{ \code{\link{logLik,km-method}}, \code{\link{km}}, \code{\link{logLikGrad}} }

\keyword{models}
\keyword{htest}

