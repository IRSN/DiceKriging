\name{checkNames}
\alias{checkNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Consistency test between the column names of two matrices
}
\description{
Tests if the names of a second matrix are equal to a given matrix up to a permutation, and permute its columns accordingly. When the second one has no column names, the names of the first one are used in the same order.
}
\usage{
checkNames(X1, X2, X1.name = "X1", X2.name = "X2")
}
\arguments{
  \item{X1}{a matrix containing column names.}
  \item{X2}{a matrix containing the same number of columns.}
  \item{X1.name}{, }
  \item{X2.name}{optional names for the matrix X1 and X2 theirselves (useful for error messages).}
}

\details{
If \code{X2} does not contain variable names, then the names of \code{X1} are used in the same order, and \code{X2} is returned with these names. Otherwise, if the column names of \code{X1} and \code{X2} are equal up to a permutation, the column of \code{X2} are permuted according to the order of \code{X1}' names.
}

\value{
The matrix \code{X2}, with columns possibly permuted. See details.
}

\author{O. Roustant}

\seealso{
\code{\link{predict,km-method}}, \code{\link{simulate,km-method}}
}

\examples{
X1 <- matrix(1, 2, 3)
X2 <- matrix(1:6, 2, 3)

colnames(X1) <- c("x1", "x2", "x3")
checkNames(X1, X2)
# attributes the same names for X2, and returns X2

colnames(X2) <- c("x1", "x2", "x5")
\dontrun{checkNames(X1, X2)}
# returns an error since the names of X1 and X2 are different

colnames(X2) <- c("x2", "x1", "x3")
checkNames(X1, X2)
# returns the matrix X2, but with permuted columns
}

\keyword{models}
