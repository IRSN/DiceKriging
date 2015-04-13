\name{goldsteinPrice}
\alias{goldsteinPrice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2D test function}
\description{
  Goldstein price 2-dimensional test function.
}
\usage{goldsteinPrice(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 2-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the Goldstein price function values at \code{x}
}

\details{
The Goldstein price function is usually defined over the domain [-2,-2] x [2, 2]. 
Here, the function is adapted to the domain [0,1] x [0,1]. It has 1 global minimum : 
x1 = c(0.5, 0.25)
}

\author{ D. Ginsbourger, Ecole des Mines de St-Etienne. }

\keyword{optimize}
\keyword{internal}


\examples{ 
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, goldsteinPrice)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid, y.grid, z.grid, 40)
x1 = c(0.5, 0.25)
points(t(x1), pch=19, col="red")
title("Fonction de Goldstein price")
}


