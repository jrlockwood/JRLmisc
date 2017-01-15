mdiaR <- function(x, target, sweights=rep(1,nrow(x)), maxdelta=1, maxiter=1000, tol=1e-08){
    ############################################################################
    ## 4/27 - dan's function for MDIA with optional base weights
    ############################################################################
    ## x - a dataframe or matrix 
    ## target - a vector of length(ncol(x))
    ## sweights -- optional sampling weights, 
    ## maxdelta -- maximum step size for Newton-Raphson
    ## maxiter -- max iteration
    ## tol -- tolerance for convergence and checking that solution yields weights that achieve balance
    ##
    ## Value:
    ## beta -- coefficients of MDIA weight function
    ## se_beta -- Standard error of beta
    ## covbeta -- variance-covariance matrix for beta
    ## weights -- MDIA weights
    ## exitcode -- 0: solution found and it yields balance; -1: no solution, probably not in convex hull; 
    ##             -2: solution that didn't yield balance
    ############################################################################
    stopifnot(is.numeric(target), ncol(x)==length(target), nrow(x)==length(sweights), maxiter > 0, maxdelta > 0, tol>0)
    .mdia <- function(x, target, sweights=rep(1,nrow(x)), maxdelta=1, maxiter=1000, tol=1e-08){
        x <- as.matrix(x)
        n <- nrow(x)
        m <- ncol(x)
        stopifnot(n > 1, m == length(target))
        
        ## set initial values for weights ##
        weight <- sweights/sum(sweights)
        
        ## set initial values for coefficents ##
        beta <- rep(0, m)
        
        z <- target
        
        delta <- rep(99, m)
        
        it <- 1
        while(max(abs(delta)) >= tol & it < maxiter){ 
            zz <- apply(x * weight, 2, sum)
            xx <- x - matrix(zz, nrow=n, ncol=m, byrow=T)
            xxw <- xx * weight
            sigma <- t(xxw) %*% xx
            delta <- qr.solve(sigma,z-zz, tol=1e-10)
            w2 <- max(abs(delta))
            if(w2 > maxdelta){ delta <- maxdelta*delta/w2}
            beta <- delta+beta
            
            weight <- as.numeric(exp(x %*% beta))*sweights;
            weight <- weight/sum(weight);
            
            it <- it+1
        }
        
        invsigma <- solve(sigma)
        covbeta <- invsigma %*% (t(xxw) %*% xxw) %*% invsigma
        stdbeta <- sqrt(diag(covbeta))
        
        weight <- as.numeric(exp(x %*% beta))
        weight <- weight/sum(weight*sweights)
        
        res <- list(beta=beta, se_beta=stdbeta, covbeta=covbeta, weights=weight)
    }
    
    res <- try(.mdia(x=x, target=target, sweights=sweights, maxdelta=maxdelta, maxiter=maxiter, tol=tol) )
    if(is.null(attr(res, "condition"))){
        if(any(abs(apply(x*res$weights*sweights, 2, sum)-target) > tol)){res$exitcode <- -2}else{
                                                                                               res$exitcode <- 0 }
    }else{
        res <- list(beta=NULL, se_beta=NULL, covbeta=NULL, weights=NULL, exitcode=-1)
    }
    
    return(res)
}
