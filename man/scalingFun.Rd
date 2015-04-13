\name{scalingFun}
\alias{scalingFun}

\title{
Scaling function 
}
\description{
Parametric transformation of the input space variables. The transformation is obtained coordinatewise by integrating piecewise affine marginal "densities" parametrized by a vector of knots and a matrix of density values at the knots. See references for more detail.
}
\usage{
scalingFun(X, knots, eta, plot=FALSE)
}

\arguments{
  \item{X}{ an n*d matrix standing for a design of n experiments in d-dimensional space }
  \item{knots}{ a list of knots parametrizing the transformation. }
  \item{eta}{ a list of coefficients parametrizing the d marginal transformations. Each element stands for a set of marginal density values at the knots defined above.}
  \item{plot}{ if TRUE plots the image of the columns of X according to the corresponding marginal transformations.}
}

\value{
The image of X by a scaling transformation of parameters knots and eta
}

\references{
Y. Xiong, W. Chen, D. Apley, and X. Ding (2007), \emph{Int. J. Numer. Meth. Engng}, A non-stationary covariance-based Kriging method for metamodelling in engineering design. 
}

\seealso{ \code{\link{scalingGrad}} }

\examples{
## 1D Transform of Xiong et al.
knots <- c(0, 0.3, 0.8, 1); eta <- c(2, 0.4, 1.4, 1.1)
nk <- length(knots)
t <- seq(from = 0, to = 1, length = 200)
f <- scalingFun(X = matrix(t), knots = list(knots), eta = list(eta))

## for text positions only
itext <- round(length(t) * 0.7)
xtext <- t[itext]; ftext <- f[itext] / 2; etamax <- max(eta)

## plot the transform function
opar <- par(mfrow = c(2, 1))
par(mar = c(0, 4, 5, 4))
plot(x = t, y = f, type = "l", lwd = 2, col = "orangered",
     main = "scaling transform f(x) and density g(x)",
     xlab = "", ylab = "", xaxt = "n", yaxt = "n") 
axis(side = 4)
abline(v = knots, lty = "dotted"); abline(h = 0)
text(x = xtext, y = ftext, cex = 1.4,
     labels = expression(f(x) == integral(g(t)*dt, 0, x)))

## plot the density function, which is piecewise linear
scalingDens1d <- approxfun(x = knots, y = eta)
g <- scalingDens1d(t)
gtext <- 0.5 * g[itext] + 0.6 * etamax
par(mar = c(5, 4, 0, 4))
plot(t, g, type = "l", lwd = 2, ylim = c(0, etamax * 1.2),
     col = "SpringGreen4", xlab = expression(x), ylab ="")
abline(v = knots, lty = "dotted")
lines(x = knots, y = eta, lty = 1, lwd = 2, type = "h", col = "SpringGreen4")
abline(h = 0)
text(x = 0.7, y = gtext, cex = 1.4, labels = expression(g(x)))

## show knots with math symbols eta, zeta
for (i in 1:nk) {
  text(x = knots[i], y = eta[i] + 0.12 * etamax, cex = 1.4,
       labels = substitute(eta[i], list(i = i)))
  mtext(side = 1, cex = 1.4, at = knots[i], line = 2.4,
        text = substitute(zeta[i], list(i = i)))
}
polygon(x = c(knots, knots[nk], knots[1]),  y = c(eta, 0, 0),
        density = 15, angle = 45, col = "SpringGreen", border = NA)
par(opar)
}

\keyword{models}
