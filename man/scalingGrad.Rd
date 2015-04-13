\name{scalingGrad}
\alias{scalingGrad}

\title{
Gradient of the dimensional Scaling function
}
\description{
Gradient of the Scaling function (marginal in dimension k) of Xiong et al. with respect to eta
}
\usage{
scalingGrad(X, knots, k)
}

\arguments{
  \item{X}{ an n*d matrix standing for a design of n experiments in d-dimensional space. }
  \item{knots}{ a list of knots parametrizing the transformation.}
  \item{k}{ dimension of the input variables for which the gradient is calculated. }
}

\value{
Gradient of the Scaling function of Xiong et al. with respect to eta
}

\references{
Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design. 
}

\seealso{ \code{\link{scalingFun}} }

\keyword{models}
