\name{kmNoNugget.init}
\alias{kmNoNugget.init}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models}
\description{
  \code{kmNoNugget.init} is used to give initial values to fit kriging models when there is no nugget effect nor noisy observations.
}

\usage{
kmNoNugget.init(model, fn, fnscale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ an object of class \code{km}. }
  \item{fn}{ the function considered: \code{\link{logLikFun}} or  \code{\link{leaveOneOutFun}}.}
  \item{fnscale}{ a real number which sign determines the direction for optimization: <0 for \code{logLikFun}, >0 for \code{leaveOneOutFun}.}
}

\value{
  \item{par}{ a matrix whose rows contain initial vectors of parameters.}
  \item{value}{ a vector containing the function values corresponding to \code{par}.}
  \item{cov}{ a list containing the covariance objects corresponding to \code{par}.}
  \item{lower}{ , }
  \item{upper}{ vectors containing lower and upper bounds for parameters.}
}

\details{
  The procedure can be summarized in 2 stages:
  \tabular{ll}{ 
   1) \tab If no initial value is provided by the user for the covariance parameters, simulate them uniformly inside the domain delimited by \code{model@lower} and \code{model@upper}. The number of simulations is the one given in \code{model@control$pop.size}. \cr
   2) \tab Compute the likelihood for each parameters set, and select the one(s) that gives the highest value(s). The number of values considered can be set by the argument \code{multistart} in \code{\link{km}}.
   	}
}


\author{ O. Roustant, David Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{km}}, \code{\link{kmEstimate}} }

\keyword{models}
\keyword{htest}
\keyword{internal}




