\name{SCAD.derivative}
\alias{SCAD.derivative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Penalty function derivative}
\description{
  Derivative of SCAD function.
}
\usage{SCAD.derivative(x, lambda)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector where the function is to be evaluated.}
  \item{lambda}{ a number representing a tuning parameter.}
}

\value{
  A vector containing the SCAD derivative values at \code{x}.
}


\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\references{ 
R. Li and A. Sudjianto (2005), Analysis of Computer Experiments Using Penalized Likelihood in Gaussian Kriging Models, \emph{Technometrics}, \bold{47} no. 2, 111-120.
 }

\seealso{ \code{\link{SCAD}}, \code{\link{km}} }

\keyword{models}
\keyword{htest}
\keyword{internal}

