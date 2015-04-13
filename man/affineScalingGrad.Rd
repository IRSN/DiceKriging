\name{affineScalingGrad}
\alias{affineScalingGrad}

\title{
Gradient of the Scaling function (affine case)
}
\description{
Gradient of the Scaling function (marginal in dimension k) of Xiong et al. with respect to eta
}
\usage{
affineScalingGrad(X, knots, k)
}

\arguments{
  \item{X}{ an n*d matrix standing for a design of n experiments in d-dimensional space }
  \item{knots}{ a (K+1) vector of knots parametrizing the transformation. The knots are here the same in all dimensions. }
  \item{k}{ dimension of the input variables for which the gradient is calculated }
}

\value{
Gradient at X of the scaling transformation of Xiong et al. with respect to eta. In the present version of the function, the gradient calculation is restricted to the case of 2 parameters by dimension (corresponding to the CovAffineScaling class).
}

\seealso{ \code{\link{scalingGrad}} }

\references{
Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design.
}

\keyword{models}
