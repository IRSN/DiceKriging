\name{covParametersBounds}
\alias{covParametersBounds}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Boundaries for covariance parameters}
\description{
  Default boundaries for covariance parameters. }
\usage{
covParametersBounds(object, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object specifying the covariance structure.}
  \item{X}{a matrix representing the design of experiments.}
}

\value{
  a list with items \code{lower, upper} containing default boundaries for the covariance parameters.
}

\details{The default values are chosen as follows : 
\tabular{ll}{
Range parameters (all covariances) \tab \code{lower=1e-10}, \code{upper}=2 times the difference between \cr
\tab the max. and min. values of \code{X} for each coordinate \cr
Shape parameters (\code{powexp} covariance) \tab  \code{lower=1e-10}, \code{upper=2} for each coordinate
}
}


\author{ Olivier Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{  \code{\link{km}} }

\keyword{models}