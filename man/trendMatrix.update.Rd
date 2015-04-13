\name{trendMatrix.update}
\alias{trendMatrix.update}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trend model matrix operation }
\description{  Updates the trend linear model matrix when new data are given. }
\usage{
trendMatrix.update(model, Xnew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  an object of class \code{km} corresponding to the non-updated kriging model. }
  \item{Xnew}{ a data frame containing the new data points.   }
}

\value{ The model matrix corresponding to known design points and new data points. }

\note{ The model design and model response are not updated. This should be done before using \code{trendMatrix.update}.}

\author{ O. Roustant, D. Ginsbourger }


\keyword{models}
\keyword{regression}
\keyword{internal}
