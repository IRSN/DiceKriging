\name{kmData}
\alias{kmData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit and/or create kriging models }
\description{ \code{kmData} is equivalent to \code{km}, except for the interface with the data. In \code{kmData}, the user must supply both the design and the response within a single data.frame \code{data}. To supply them separately, use \code{km}.
}
\usage{
kmData(formula, data, inputnames = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class "formula" specifying the linear trend of the kriging model (see \code{\link{lm}}). At this stage, transformations of the response are not taken into account.}
  \item{data}{ a data.frame containing both the design (input variables) and the response (1-dimensional output given by the objective function at the design points).}
  \item{inputnames}{ an optional vector of character containing the names of variables in \code{data} to be considered as input variables. By default, all variables but the response are input variables.}
  \item{\dots}{ other arguments for creating or fitting Kriging models, to be taken among the arguments of \code{km} function apart from \code{design} and \code{response}.}
}
\value{
An object of class \code{km} (see \code{\link{km-class}}).
}
\author{
O. Roustant
}

\seealso{
\code{\link{km}}
}

\examples{
# a 16-points factorial design, and the corresponding response
d <- 2; n <- 16
design.fact <- expand.grid(x1=seq(0,1,length=4), x2=seq(0,1,length=4))
y <- apply(design.fact, 1, branin)
data <- cbind(design.fact, y=y)

# kriging model 1 : matern5_2 covariance structure, no trend, no nugget effect
m1 <- kmData(y~1, data=data)
# this is equivalent to: m1 <- km(design=design.fact, response=y)

# now, add a second response to data:
data2 <- cbind(data, y2=-y)
# the previous model is now obtained with:
m1_2 <- kmData(y~1, data=data2, inputnames=c("x1", "x2"))

}
\keyword{models}
\keyword{htest}
