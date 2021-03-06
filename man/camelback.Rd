\name{camelback} 
\alias{camelback}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2D test function}
\description{
  Camelback 2-dimensional test function.
}
\usage{camelback(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 2-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the Camelback function values at \code{x}
}

\details{
The Camelback function is usually defined over the domain [-3,-2] x [3, 2]. 
Here, the function is adapted to the domain [0,1] x [0,1]. It has 2 global minima : 
x1 = c(0.5149730,0.3218374); x2 = c(0.4850263,0.6781641)
}

\author{ D. Ginsbourger, Ecole des Mines de St-Etienne. }

\keyword{optimize}
\keyword{internal}

\examples{ 
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, camelback)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,20)
x1 = c(0.5149730,0.3218374); x2 = c(0.4850263,0.6781641)
points(rbind(t(x1), t(x2)), pch=19, col="red")
title("Fonction Camelback")
}


