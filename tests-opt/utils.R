require(DiceKriging)
require(DiceView)

trignD <- function(x) {
    prod(sin(2*pi*x))
}

merge.list <- function(...){
    args <- list(...)
    l = args[[1]]
    for (i in 2:length(args)){
        for (arg in names(args[[i]])){
            l[[arg]]=args[[i]][[arg]]
        }
    }
    return(l)
}

# expand.grid.rows <- function(add.seq,prev.seq) {
#     add.seq = as.matrix(add.seq)
#     m = matrix(prev.seq,nrow=nrow(prev.seq)*nrow(add.seq),ncol=ncol(prev.seq)+1)
#     for (i in 1:nrow(add.seq)) {
#         m[((i-1)*nrow(prev.seq)+1):(i*nrow(prev.seq)), ncol(prev.seq)+1] = add.seq[i,1]  
#     }
#     return(m)
# }
make.grid <- function(s=c(0,.5,1),d=1) {
    l=list(s)
    if (d>1) for (i in 2:d) {
        l=c(l,list(s))
    }
    return(do.call("expand.grid",l))
}

build.km <- function(d=3, npoints=5, f=trignD, seed=1, plot=FALSE,control=NULL,...) {
    
    # function to create a kriging based on a deterministic function (f) and a factorial sampling of npoints
    km.nD <- function(d, f, npoints, control=NULL,...) {
        control = merge.list(list(trace=FALSE),control)
        
        design.fact <- matrix(runif(d*npoints^d),ncol=d)#make.grid(s=seq(0,1,length=npoints),d=d)
        design.fact <- data.frame(design.fact); 
        names(design.fact) <- paste("x",1:d,sep="")
        
        model.nD = km(design=design.fact, response=apply(FUN=f,X=design.fact,MARGIN=1),control=control,...)
        
        return(model.nD)
    }
    
    set.seed(seed)
    
    model <- NULL
    try(model <- km.nD(d=d,f=f,npoints=npoints,control=control,...),silent=TRUE)
    if (!is.null(model)) {
        if(plot) sectionview.km(model,center=rep(0.5,d))
    }
    
    return(model)
}

# model = build.km(d=2, f=trignD, npoints=7, seed=1, control=list(trace=FALSE), plot=TRUE)
# if ( is.null(model) || any(model@covariance@range.val<0.1) ) {
#     stop("[FAILED] utils.R: failed to build kriging model.") 
# } else cat("[PASSED] utils.R\n")


estim.theta <- function(d,npoints,f=trignD,n=10,control=NULL,...){
    control = merge.list(list(trace=FALSE),control)
    
    theta.estim = matrix(NaN,n,d)
    pb <- txtProgressBar(min = 1, max = n, style = 3)
    for (i in 1:nrow(theta.estim)) {
        setTxtProgressBar(pb, i)
        
        model = NULL        
        model = build.km(d=d,npoints=npoints,f=f,seed=i,plot=FALSE,control=control,...)
        if(!is.null(model))
            theta.estim[i,] = model@covariance@range.val
        else cat("!")
    }
    cat("\n")
    
    return(theta.estim)
}

# theta.estim = estim.theta(d=1,npoints=10,n=10)
# ok = length(which(apply(FUN=min,X=theta.estim,MARGIN=1)>0.1))/nrow(theta.estim)
# if ( ok != 1 ) {
#     stop("[FAILED] utils.R: failed to estimate theta.") 
# } else cat("[PASSED] utils.R\n")

mse.model <- function(model,d,f,npoints=50) {
    #n=1000
    #X = seq(0,1,length=n) # should be Monte Carlo / qMC ...
    
    X <- make.grid(seq(0,1,length=npoints),d=d)
    X <- data.frame(X); 
    names(X) <- paste("x",1:d,sep="")
    
    pred = predict.km(object=model,X,type="SK")
    m = pred$mean
    sd = pred$sd
    Y=apply(FUN=f,X=X,MARGIN=1)
    mse = (sum((Y-m)^2) + sum(sd^2)) / npoints^d
    return(mse)
}

mse <- function(d,npoints,f=trignD,n=10,control=NULL, ...){
    mse = matrix(NaN,n,1)
    pb <- txtProgressBar(min = 1, max = n, style = 3)
    for (i in 1:nrow(mse)) {
        setTxtProgressBar(pb, i)
        
        model = NULL
        model = build.km(d=d,npoints=npoints,f=f,seed=i,plot=FALSE, 
                         control=control,...)
        
        if(!is.null(model))        
            mse[i,1] = mse.model(model,d=d,f=f)
        else cat("!")
        
    }
    cat("\n")
    
    return(mse)
}

# if ( max(mse(d=1,npoints=10,n=10)) > .1 ) {
#     stop("[FAILED] utils.R: failed to build good (low IMSE) kriging model.") 
# } else cat("[PASSED] utils.R\n")
