winbugs.wishart.params <- function(target, nsim, reltol = 0.05){
    ## "target" is a diagonal covariance matrix passed as a vector of the
    ## diagonals.  winbugs.df is set to one greater than the dimension of target.
    ## the goal is then to find "winbugs.scalemat" such that samples of precision
    ## matrices from this wishart have prior MEDIAN of diagonals of the covariance
    ## matrix (its inverse) very close to "target".  "nsim" is number of monte
    ## carlo simulations used to approximate medians.  "reltol" is relative
    ## tolerance to determine we are close enough.
    ##
    ## NOTE: use median rather than mean because of extreme skew in distribution
    ##
    ## NOTE: winbugs parameterization of wishart has prior mean equal to
    ## winbugs.df * solve(winbugs.scalemat), see
    ## http://statacumen.com/2009/07/02/wishart-distribution-in-winbugs-nonstandard-parameterization/
    ## which means that our final result needs to be the inverse of whatever scalemat we locate
    stopifnot( all(target > 0) && (reltol > 0) )
    winbugs.df <- length(target) + 1
    
    ## rwish() is parameterized so that E(W) = vS and is supposed to a prior for
    ## precision matrix with mean of its inverse = diag(target), so start scalemat
    ## at a bad moment estimate
    scalemat <- (1 / winbugs.df) * diag(1 / target)
    ## note that solve(E(W)) = solve(winbugs.df * scalemat) = diag(target)
    
    simV <- vector(nsim, mode="list")  
    done <- FALSE
    while(!done){
        .scalemat <- scalemat
        print(diag(.scalemat))
        
        ## generate "nsim" *covariance* matrices from the wishart
        for(i in 1:nsim){
            simV[[i]] <- solve(rWishart(1, df = winbugs.df,  Sigma = scalemat)[,,1])
        }
        
        ## get the medians of the diagonals and compare to target
        d <- apply(do.call("rbind",lapply(simV, diag)), 2, median) - target
        absdif <- abs( d / target )
        done <- all( absdif < reltol )
        ## if we are too big, need to increase diagonals, if we are too small, need to reduce
        diag(scalemat) <- diag(scalemat) * ifelse(d > 0, (1 + absdif), 1 / (1 + absdif))
    }
    
    ## print summaries of variances and correlations
    .V <- do.call("rbind",lapply(simV, diag))
    .R <- do.call("rbind",lapply(simV, function(x){ .tmp <- diag(1/sqrt(diag(x))) %*% x %*% diag(1/sqrt(diag(x))); .tmp[lower.tri(.tmp)] }))
    colnames(.V) <- paste("var",1:ncol(.V), sep="")
    colnames(.R) <- paste("cor",1:ncol(.R), sep="")
    print(apply(.V, 2, quantile, prob = c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975)), digits=4)
    print(apply(.R, 2, quantile, prob = c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975)), digits=4)
    
    return(list(winbugs.df = winbugs.df, winbugs.scalemat = solve(.scalemat)))
}
## test
## winbugs.wishart.params(target = c(10,4,4,8), nsim = 20000, reltol = 0.02)
## winbugs.wishart.params(target = c(1,2,10,3), nsim = 30000, reltol = 0.02)
