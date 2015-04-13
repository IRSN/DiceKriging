# To test that low dimensions kriging are well working
# @see http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

require("testthat")

success.rate <- function(d,npoints,plot=FALSE,msg=FALSE,...){
    theta.estim = estim.theta(d,npoints,...)
    
    ok = length(which(apply(FUN=min,X=theta.estim,MARGIN=1)>0.1))/nrow(theta.estim)
    
    if(msg){
        cat(paste(sep="","Kriging ",d,"D function, on ",npoints," conditional points", list(...),"\n"))
        cat(paste(100-100*ok,"% failure\n"))
    }
    
    if(plot) hist(apply(FUN=min,X=theta.estim,MARGIN=1),xlab="min(theta[1,d])")
    
    return(ok)
}

context("Testing kriging success")

test_that(desc="Kriging:1D,10 cond. points",expect_that(success.rate(d=1, npoints=10),equals(1)))

test_that(desc="Kriging:2D,100 cond. points",expect_that(success.rate(d=1, npoints=10),equals(1)))

test_that(desc="Kriging:3D,343 cond. points",expect_that(success.rate(d=1, npoints=7),equals(1)))


