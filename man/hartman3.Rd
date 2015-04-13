\name{hartman3}
\alias{hartman3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 3D test function}
\description{
  Hartman 3-dimensional test function.
}
\usage{hartman3(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 3-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the hartman3 function values at \code{x}
}

\details{
The hartman3 function is defined over the domain \code{[0,1]^3}. It has 1 global minimum : 
x1 = c(0.1, 0.55592003, 0.85218259)
}

\author{ D. Ginsbourger, Ecole des Mines de St-Etienne. }

\keyword{optimize}
\keyword{internal}


\examples{ 
design <- matrix(runif(300), 100, 3)
response <- apply(design, 1, hartman3)
}


