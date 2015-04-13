\name{SCAD}
\alias{SCAD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Penalty function}
\description{
  Smoothly Clipped Absolute Deviation function.
}
\usage{SCAD(x, lambda)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector where the function is to be evaluated.}
  \item{lambda}{ a number representing a tuning parameter.}
}

\value{
  A vector containing the SCAD values at \code{x}.
}

\details{
 SCAD  is an even continuous function equal to 0 at \code{x=0}, and defined piecewise with derivative \code{lambda} in \code{[0, lambda]}, \code{(a*lambda - x)/(a-1)} in \code{[lambda, a*lambda]}, and \code{0} for \code{x} larger than  \code{a*lambda}. As suggested by (Li, Sudjianto, 2005), we set \code{a=3.7}.
}

\note{
 In MLE problems, the penalty value \code{lambda} should tend to 0 when the sample size tends to infinity to insure that the asymptotic properties of Penalized-MLE and MLE are the same (see Li, Sudjianto, 2005).
}

\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\references{ 
R. Li and A. Sudjianto (2005), Analysis of Computer Experiments Using Penalized Likelihood in Gaussian Kriging Models, \emph{Technometrics}, \bold{47} no. 2, 111-120.
}

\seealso{ \code{\link{SCAD.derivative}} and \code{\link{km}} for a famous example }

\keyword{models}
\keyword{htest}

\examples{ 
x <- seq(-8,8, length=200)
a <- 3.7

lambda <- 1.5
y <- SCAD(x, lambda)
plot(x, y, type="l", ylim=c(0,6))
x.knots <- c(-a*lambda, -lambda, 0, lambda, a*lambda)
points(x.knots, SCAD(x.knots, lambda), pch=19, cex=0.5)
text(6, SCAD(6, lambda)+0.3, paste("lambda =", lambda))

for (i in 1:2) {
   lambda <- lambda - 0.5
   y <- SCAD(x, lambda)
   lines(x, y, type="l")
   x.knots <- c(-a*lambda, -lambda, 0, lambda, a*lambda)
   points(x.knots, SCAD(x.knots, lambda), pch=19, cex=0.5)
   text(6, SCAD(6, lambda)+0.3, paste("lambda =", lambda))
}

abline(v=0, h=0, lty="dotted")
title("SCAD function")
}

