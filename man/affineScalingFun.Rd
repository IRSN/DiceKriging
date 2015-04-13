\name{affineScalingFun}
\alias{affineScalingFun}

\title{
Scaling function (affine case)
}
\description{
Parametric transformation of the input space variables. The transformation is obtained coordinatewise by integrating piecewise affine marginal "densities" parametrized by a vector of knots and a matrix of density values at the knots. See references for more detail.
}
\usage{
affineScalingFun(X, knots, eta)
}

\arguments{
  \item{X}{ an n*d matrix standing for a design of n experiments in d-dimensional space }
  \item{knots}{ a (K+1) vector of knots parametrizing the transformation. The knots are here the same in all dimensions. }
  \item{eta}{ a d*(K+1) matrix of coefficients parametrizing the d marginal transformations. Each line stands for a set of (K+1) marginal density values at the knots defined above.}
}

\value{
The image of X by a scaling transformation of parameters knots and eta
}

\seealso{ \code{\link{scalingFun}} }

\references{
Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design. 
}

\keyword{models}
