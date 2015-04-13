\name{kmEstimate}
\alias{kmEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting Kriging Models }
\description{
  \code{kmEstimate} is used to fit kriging models. This function should not be called directly, due to the environments defined in \code{km} to avoid computing twice \code{nxn} matrices. Call \code{km} instead.
}
\usage{
kmEstimate(model, envir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  an object of class \code{km}. }
  \item{envir}{  an environment specifying where to assign intermediate values for future gradient calculations. }
}

\value{
  An object of class \code{km}.
}


\references{ Park J-S, Baek J. (2001), Efficient computation of maximum likelihood estimators
                      in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences},
                      \bold{27}, 1-7. }
\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{kmNoNugget.init}}, \code{\link{km}} }

\keyword{models}
\keyword{htest}
\keyword{internal}



