require("DiceKriging")
require("testthat")
source("test-km.R")

set.seed(1)

# a 16-points factorial design, and the corresponding response
d <- 2; n <- 16
design.fact <- expand.grid(x1=seq(0,1,length=4), x2=seq(0,1,length=4))
y <- apply(design.fact, 1, branin)

context("Checking km scaling: 2D example - Branin-Hoo function")

# kriging model 1 : matern5_2 covariance structure, no trend, no nugget effect
m1 <- km(design=design.fact, response=y,control=list(trace=FALSE),scaling=TRUE)

test_that.km(m1,trend.coef = 279.2024,covariance.sd2 = 99777.04,covariance.eta = matrix(c(2.1286391 , 0.5116505 , 0.5, 0.5),ncol=2))


context("Checking km scaling: different number of nodes per dimension (not passing with DiceKriging 1.5-4)")

design.fact <- expand.grid(aa=seq(0,1,length=4), ba=seq(0,1,length=4),ac=seq(0,1,length=4))
y <- apply(design.fact, 1, function(x)branin(x[1:2])*x[3])

set.seed(12)
# kriging model 1 : matern5_2 covariance structure, no trend, no nugget effect
m2 <- km(design=design.fact, response=y,control=list(trace=FALSE),scaling=TRUE, knots = list(c(0,1),c(0,1),c(0,.5,1)))

test_that.km(m2,trend.coef = 127.2252,covariance.sd2 = 85424.41)
