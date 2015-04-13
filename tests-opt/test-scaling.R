# To test that scaling=TRUE increases kriging predicting power (in the sense of MSE)
# @see http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

require("testthat")

scaling_surpass.mse <- function(d,npoints,f,control=NULL,plot=FALSE,msg=FALSE,...){
    mse_noscaling = mse(d=d,npoints=npoints,f=f,
                            scaling=FALSE,
                            control=control,...)
    
    mse_scaling = mse(d=d,npoints=npoints,f=f,
                          scaling=TRUE,
                          control=control,...)
    
    if (msg){
    cat(paste(sep="","Kriging ",d,"D function, on ",npoints," conditional points, ",list(...),"\n"))
    cat("NO scaling -> ")
    cat(summary(mse_noscaling,na.rm=TRUE))
    cat("\nscaling -> ")
    cat(summary(mse_scaling,na.rm=TRUE))
    cat("\n")
    }
    
    if(plot) boxplot(cbind(mse_noscaling,mse_scaling),names=c("MSE without scaling","MSE with scaling"))
    
    return ( (mean(mse_noscaling) - mean(mse_scaling))/mean(mse_noscaling) )
}

context("Testing kriging MSE with/out scaling")

test_that(desc="scaling:2D,25 cond. points",expect_that(sign(scaling_surpass.mse(d=2, npoints=5, f = trignD) + .1),equals(1)))

test_that(desc="scaling:2D,25 cond. points,x<-x^5",expect_that(sign(scaling_surpass.mse(d=2, npoints=5, f = function(x){trignD(x^5)})),equals(1)))

test_that(desc="scaling:2D,100 cond. points",expect_that(sign(scaling_surpass.mse(d=2, npoints=10, f = trignD) + .1),equals(1)))

test_that(desc="scaling:2D,100 cond. points,x<-x^5",expect_that(sign(scaling_surpass.mse(d=2, npoints=10, f = function(x){trignD(x^5)})),equals(1)))
