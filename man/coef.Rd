\name{coef}
\alias{coef}
\alias{coef<-}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Get coefficients values }
\description{ Get or set coefficients values.}
\usage{
  coef(object, ...)
  coef(object, ...) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{object}{an object specifying a covariance structure or a km object.}
%   \item{type}{a character string indicating what to get. Default is 'all'.}
%   \item{as.list}{a boolean indicating the output type. Default is \code{TRUE}.}
   \item{...}{other arguments (undocumented at this stage).}
   \item{value}{value for replacement. At this stage, only for the class \code{covTensorProduct}.}
}

%\value{ A list or vector containing the desired coefficients values.}
\note{The replacement method \code{coef<-} is in construction.}

\author{Y. Deville, O. Roustant}

\keyword{models}
