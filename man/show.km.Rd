\name{show}
\alias{show,km-method}

\title{ Print values of a km object}
\description{
  Show method for \code{km} object. Printing the main features of a kriging model.
}

\usage{
   \S4method{show}{km}(object)
}

\arguments{
  \item{object}{ an object of class \code{km}.}
}

\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{ \code{\link{km}} }


\examples{ 
# A 2D example - Branin-Hoo function

# a 16-points factorial design, and the corresponding response
d <- 2; n <- 16
fact.design <- expand.grid(seq(0,1,length=4), seq(0,1,length=4))
fact.design <- data.frame(fact.design); names(fact.design)<-c("x1", "x2")
branin.resp <- data.frame(branin(fact.design)); names(branin.resp) <- "y" 

# kriging model 1 : power-exponential covariance structure, no trend, 
#                   no nugget effect
m1 <- km(y~1, design=fact.design, response=branin.resp, covtype="powexp")
m1    # equivalently : show(m1)
}


\keyword{methods}

