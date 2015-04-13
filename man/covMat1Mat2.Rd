\name{covMat1Mat2}
\alias{covMat1Mat2} 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cross covariance matrix }
\description{
  Computes the cross covariance matrix between two sets of locations for a spatial random process with a given covariance structure. Typically the two sets are a learning set and a test set.
}
\usage{
covMat1Mat2(object, X1, X2, nugget.flag=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object specifying the covariance structure.}
  \item{X1}{ a matrix whose rows represent the locations of a first set (for instance a set of learning points).}
  \item{X2}{ a matrix whose rows represent the locations of a second set (for instance a set of test points).}
  \item{nugget.flag}{an optional boolean. If \code{TRUE}, the covariance between 2 equal locations takes into account the nugget effect (if any). Locations are considered equal if their euclidian distance is inferior to \code{1e-15}. Default is \code{FALSE}.}
}

\value{
  a matrix of size \code{(nb of rows of X1 * nb of rows of X2)} whose element \code{(i1,i2)} is equal to the covariance between the locations specified by row \code{i1} of \code{X1} and row \code{i2} of \code{X2}.
  }


\author{ Olivier Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{  \code{\link{covMatrix}} }

\keyword{models}