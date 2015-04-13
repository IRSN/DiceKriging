# For each covtype, test that 
# * theta estimation is well performed (in the sense that a given theta is not to badly estimated: theta/10 < mean(theta_hat) < 10*theta ) 
# * MSE is not too hight
# @see http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

require("testthat")

build.fun <- function(d, covtype, coef.trend=0, coef.cov = .1, coef.var=1, control=NULL,plot=FALSE, n.grid=10,...){
    X <-data.frame(make.grid(c(0,1),d))
    names(X) <- paste("x",1:d,sep="")
    
    k <- km(design=X, response=(rep(0,2^d)), 
            covtype=covtype, coef.trend=coef.trend, coef.cov = rep(coef.cov,d), coef.var=rep(coef.var,d), 
            control=control,...)
        
    if (plot) sectionview.km(k,center=rep(0.5,d))
    
    # create an uncond. process with these properties
    n.grid = min(n.grid,floor(40^(1/d))) # To limit simulate points to 40 (numerical issues)
    .k.X = as.matrix(make.grid(seq(0,1,length=n.grid),d))
    .k.Y = simulate(k,cond=FALSE,newdata=.k.X,checkNames=FALSE)
    
    # create a basic interpolation function from the previous uncond. process
    f <- function(x) {
        x=array(x)
        dist= (rowSums((as.matrix(.k.X,ncol=d)-t(matrix(as.matrix(x,ncol=d),ncol=nrow(.k.X),nrow=d)))^2))^.5
        
        if(min(dist)<1E-6) return(.k.Y[which.min(dist)])
        
        return(sum(.k.Y/dist)/sum(1/dist))
                   
#         sort.dist=sort(dist, index.return = TRUE)
#                    
#         if (sort.dist$x[1]<1E-6) return(.k.Y[sort.dist$ix[1]])
#                    
#         x1=.k.X[sort.dist$ix[1],]
#         x2=.k.X[sort.dist$ix[2],]
#         y1=.k.Y[sort.dist$ix[1]]
#         y2=.k.Y[sort.dist$ix[2]]
#         return( (y1*sort.dist$x[2]+y2*sort.dist$x[1]) / (sort.dist$x[1]+sort.dist$x[2]) )
    }
    
    if (plot) {
        sectionview.fun(fun=f,center=rep(0.5,d),add=TRUE,npoints=10*n.grid,col='red')
        if (d==1) points(.k.X,.k.Y,col='red')
    }
    
    return(f)
}
# plot(Vectorize(build.fun(d=1,covtype="gauss")))
# plot(Vectorize(build.fun(d=1,covtype="matern3_2")))
# plot(Vectorize(build.fun(d=1,covtype="exp")))
# sectionview3d.fun(build.fun(d=2,covtype="gauss",coef.cov=.9))

context("Testing kriging range estimation")
theta_estim <- function(d,npoints,covtype,theta,control=NULL,plot=FALSE,msg=FALSE,...) {
    
    f = build.fun(d=d,
                  covtype=covtype,coef.cov=theta,
                  control=control,plot=plot,n.grid=npoints,...)
    
    estim = estim.theta(d=d,npoints=npoints,f=f,
                        covtype=covtype,
                        control=control,...)
    
    # to clean bad kriging (= theta << 1)
    estim = estim[which(rowSums(estim)>0.0001),]
    
    if (msg) {
        cat(paste(sep="","Kriging (cov ",covtype,") ",d,"D function, on ",npoints," conditional points, ",list(...),"\n"))
        print(summary(estim,na.rm=TRUE))
    }
    
    if (plot) hist(estim)
    
    return(estim)
}

covtypes=c("gauss","exp","matern3_2","matern5_2")
for (covtype in covtypes) {
    theta=0.1
    
    theta_est = theta_estim(d=1, npoints=20, covtype = covtype, theta = theta, control=list(logLikFailOver=TRUE))
    test_that(desc=paste("CovEstimTheta:1D,covtype=",covtype,",theta=",theta),expect_that( sign(mean(theta_est) - 0.1*theta),equals(1)))
    test_that(desc=paste("CovEstimTheta:1D,covtype=",covtype,",theta=",theta),expect_that( sign(10*theta - mean(theta_est)),equals(1)))
    
    theta_est = theta_estim(d=2, npoints=10, covtype = covtype, theta = theta, control=list(logLikFailOver=TRUE))
    test_that(desc=paste("CovEstimTheta:2D,covtype=",covtype,",theta=",theta),expect_that( sign(mean(theta_est) - 0.1*theta),equals(1)))
    test_that(desc=paste("CovEstimTheta:2D,covtype=",covtype,",theta=",theta),expect_that( sign(10*theta - mean(theta_est)),equals(1)))
}

context("Testing kriging MSE")
mse_estim <- function(d,npoints,covtype,theta,control=NULL,plot=FALSE,msg=FALSE,...){
    
    f = build.fun(d=d,
                  covtype=covtype,
                  control=control,plot=plot,...)
    
    estim = mse(d=d,npoints=npoints,f=f,
                covtype=covtype,
                control=control,...)
    
    if (msg) {
        cat(paste(sep="","Kriging (cov ",covtype,") ",d,"D function, on ",npoints," conditional points, ",list(...),"\n"))
        print(summary(estim,na.rm=TRUE))
    }
    
    if (plot) hist(estim)    
    
    return(estim)
}

covtypes=c("gauss","exp","matern3_2","matern5_2")
for (covtype in covtypes) {
    theta=0.1
    mse_est = mse_estim(d=1, npoints=20, covtype = covtype, theta = theta, control=list(logLikFailOver=TRUE))
    test_that(desc=paste("CovEstimMSE:1D,covtype=",covtype,",theta=",theta),expect_that( sign(mean(mse_est,na.rm=TRUE) - 10),equals(-1)))
}
covtypes=c("gauss","exp","matern3_2","matern5_2")
for (covtype in covtypes) {
    theta=0.1
    mse_est = mse_estim(d=2, npoints=10, covtype = covtype, theta = theta, control=list(logLikFailOver=TRUE))
    test_that(desc=paste("CovEstimMSE:2D,covtype=",covtype,",theta=",theta),expect_that( sign(mean(mse_est,na.rm=TRUE) - 10),equals(-1)))
}


