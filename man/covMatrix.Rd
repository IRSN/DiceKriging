\name{covMatrix} 
\alias{covMatrix}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Covariance matrix}
\description{
  Computes the covariance matrix at a set of locations for a spatial random process with a given covariance structure.
}

\usage{
covMatrix(object, X, noise.var = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{an object specifying the covariance structure.}
 \item{X}{ a matrix whose columns represent locations. }
 \item{noise.var}{for noisy observations : an optional vector containing the noise variance at each observation}
}

\value{
  a list with the following items :
  \item{C}{ a matrix representing the covariance matrix for the locations specified in the X argument, including a possible nugget effect or observation noise.}
  \item{vn}{ a vector of length \code{n} (\code{X} size) containing a replication of the nugget effet or the observation noise (so that \code{C-diag(vn)} contains the covariance matrix when there is no nugget effect nor observation noise)}
  }


\author{ Olivier Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{  \code{\link{covMatrixDerivative}} }

\keyword{models}