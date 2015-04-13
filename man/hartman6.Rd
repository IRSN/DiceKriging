\name{hartman6}
\alias{hartman6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 6D test function}
\description{
  Hartman 6-dimensional test function.
}
\usage{hartman6(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 6-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the hartman6 function values at \code{x}
}

\details{
The hartman6 function is defined over the domain \code{[0,1]^6}. It has 1 global minimum : 
x1 = c(0.20168952, 0.15001069, 0.47687398, 0.27533243, 0.31165162, 0.65730054)
}

\author{ D. Ginsbourger, Ecole des Mines de St-Etienne. }

\keyword{optimize}
\keyword{internal}


\examples{ 
design <- matrix(runif(600), 100, 6)
response <- apply(design, 1, hartman6)
}


