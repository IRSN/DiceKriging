\name{computeAuxVariables} 
\alias{computeAuxVariables}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Auxiliary variables for kriging }

\description{
Computes or updates some auxiliary variables used for kriging (see below). This is useful in several situations : when all parameters are known (as for one basic step in Bayesian analysis), or when some new data is added but one does not want to re-estimate the model coefficients. On the other hand, \code{computeAuxVariables} is not used during the estimation of covariance parameters, since this function requires to compute the trend coefficients at each optimization step; the alternative given by (Park, Baek, 2001) is preferred.
} 

\usage{
computeAuxVariables(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ an object of class \code{km} with missing (or non updated) items. }
}

\value{
   An updated \code{km} objet, where the changes concern the following items:
  \item{T}{ a matrix equal to the upper triangular factor of the Choleski decomposition of \code{C}, such that \code{t(T)*T = C} (where C is the covariance matrix). }
  \item{z}{ a vector equal to \code{inv(t(T))*(y - F*beta)}, with \code{y}, \code{F}, \code{beta} are respectively the response, the experimental matrix and the trend coefficients specified in \code{model@trend.coef}. If \code{model@trend.coef} is empty, \code{z} is not computed.}
  \item{M}{ a matrix equal to \code{inv(t(T))*F}. }
  }
  
\note{
\code{T} is computed with the base function \code{chol}. \code{z} and \code{M} are computed by solving triangular linear systems with \code{backsolve}. \code{z} is not computed if \code{model@trend.coef} is empty.
}

\references{ 
J.-S. Park and J. Baek (2001), Efficient computation of maximum likelihood estimators in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences}, \bold{27} no. 1, 1-7. 
}
                      
\author{O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne}

\seealso{ \code{\link{covMatrix}}, \code{\link[base]{chol}}, \code{\link[base]{backsolve}}. }

\keyword{models}
