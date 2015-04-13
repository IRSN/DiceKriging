\name{plot}

%- Also NEED an '\alias' for EACH other topic documented here.
\alias{plot}
\alias{plot.km}
\alias{plot,km-method}
\title{Diagnostic plot for the validation of a km object}

\description{
Three plots are currently available, based on the \code{leaveOneOut.km} results: one plot of fitted values against response values, one plot of standardized residuals, and one qqplot of standardized residuals.
}
\usage{
\S4method{plot}{km}(x, y, kriging.type = "UK", trend.reestim = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class "km" without noisy observations.}
  \item{y}{ not used.  }
  \item{kriging.type}{ an optional character string corresponding to the kriging family, to be chosen between simple kriging ("SK") or universal kriging ("UK").}
  \item{trend.reestim}{ should the trend be reestimated when removing an observation? Default to FALSE.}
  \item{...}{ no other argument for this method.}
}

\details{
The diagnostic plot has not been implemented yet for noisy observations. The standardized residuals are defined by \code{( y(xi) - yhat_{-i}(xi) ) / sigmahat_{-i}(xi)}, where \code{y(xi)} is the response at the point \code{xi}, \code{yhat_{-i}(xi)} is the fitted value when removing the observation \code{xi} (see \code{\link{leaveOneOut.km}}), and \code{sigmahat_{-i}(xi)} is the corresponding kriging standard deviation. }

\section{Warning}{
 Kriging parameters are not re-estimated when removing one observation. With few points, the re-estimated values can be far from those obtained with the entire learning set. One option is to reestimate the trend coefficients, by setting \code{trend.reestim=TRUE}.
}

\value{ A list composed of:
  \item{mean }{ a vector of length \emph{n}. The ith coordinate is equal to the kriging mean (including the trend) at the ith observation number when removing it from the learning set, }
  \item{sd }{ a vector of length \emph{n}. The ith coordinate is equal to the kriging standard deviation at the ith observation number when removing it from the learning set,}
where \emph{n} is the total number of observations.
}


\references{ 

N.A.C. Cressie (1993), \emph{Statistics for spatial data}, Wiley series in probability and mathematical statistics.

J.D. Martin and T.W. Simpson (2005), Use of kriging models to approximate deterministic computer models, \emph{AIAA Journal}, \bold{43} no. 4, 853-863.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.
	
}

\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{  \code{\link{predict,km-method}},  \code{\link{leaveOneOut.km}} }

\keyword{models}
\keyword{methods}

\examples{
# A 2D example - Branin-Hoo function

# a 16-points factorial design, and the corresponding response
d <- 2; n <- 16
fact.design <- expand.grid(seq(0,1,length=4), seq(0,1,length=4))
fact.design <- data.frame(fact.design); names(fact.design)<-c("x1", "x2")
branin.resp <- data.frame(branin(fact.design)); names(branin.resp) <- "y" 

# kriging model 1 : gaussian covariance structure, no trend, 
#                   no nugget effect
m1 <- km(~.^2, design=fact.design, response=branin.resp, covtype="gauss")
plot(m1)  # LOO without parameter reestimation
plot(m1, trend.reestim=TRUE)  # LOO with trend parameters reestimation 
                              # (gives nearly the same result here)
}

