# To test that logLikFailover is increasing kriging success (success=theta is not estimated to NaN or 0+)
# @see http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

require("testthat")

failover_surpass.rate <- function(d,npoints,plot=FALSE,msg=FALSE,...){
    theta.estim_nofailover = estim.theta(d=d,npoints=npoints,
                                         control=list(logLikFailOver=FALSE),
                                         ...)
    
    theta.estim_failover = estim.theta(d=d,npoints=npoints,
                                       control=list(logLikFailOver=TRUE),
                                       ...)
    
    if (plot) boxplot(cbind(apply(FUN=min,X=theta.estim_nofailover,MARGIN=1),apply(FUN=min,X=theta.estim_failover,MARGIN=1)),names=c("Theta without fail-over","Theta with fail-over"))
    
    ok_nofailover = length(which(apply(FUN=min,X=theta.estim_nofailover,MARGIN=1)>0.1))/nrow(theta.estim_nofailover)
    ok_failover = length(which(apply(FUN=min,X=theta.estim_failover,MARGIN=1)>0.1))/nrow(theta.estim_failover)
    
    if (msg){
        cat(paste(sep="","Kriging ",d,"D function, on ",npoints," conditional points, ",list(...),"\n"))
        cat(paste("NO loglikFailover -> ",100-100*ok_nofailover,"% failure\n"))
        cat(paste("loglikFailover -> ",100-100*ok_failover,"% failure\n"))
    }
    
    return(ok_failover-ok_nofailover)
}

context("Testing kriging success with/out logLik failover")

test_that(desc="loglikFailOver:1D,10 cond. points,Gauss cov",expect_that(sign(failover_surpass.rate(d=1, npoints=10,covtype='gauss')),equals(1)))

test_that(desc="loglikFailOver:2D,100 cond. points,Gauss cov",expect_that(sign(failover_surpass.rate(d=2, npoints=10,covtype='gauss')),equals(1)))

