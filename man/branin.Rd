\name{branin} 
\alias{branin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2D test function}
\description{
  Branin-Hoo 2-dimensional test function.
}
\usage{branin(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 2-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the Branin-Hoo function values at \code{x}
}

\details{
The Branin-Hoo function is defined here over [0,1] x [0,1], instead of [-5,0] x [10,15] as usual. 
It has 3 global minima : 
x1 = c(0.9616520, 0.15); x2 = c(0.1238946, 0.8166644); x3 = c(0.5427730, 0.15)
}

\author{ D. Ginsbourger, Ecole des Mines de St-Etienne. }

\examples{ 
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, branin)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,40)
x1 = c(0.9616520, 0.15); x2 = c(0.1238946, 0.8166644); x3 = c(0.5427730, 0.15)
points(rbind(t(x1), t(x2), t(x3)), pch=19, col="red")
title("Fonction de Branin")
}

\keyword{optimize}
\keyword{internal}