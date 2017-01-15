gchol2 <- function(x, stol=1e-10, ltol=1e-8, rtol=1e-4, check_target_rank = 0){
    ## lower triangular cholesky for (d x d) positive semi-definite symmetric matrix.
    ## alternative to gchol(bdsmatrix) that seems to function better
    ##
    ## stol = tolerance for determining symmetry
    ## ltol = tolerance for setting column of L to zero and also for eigenvalues of x
    ## rtol = tolerance for difference between L %*% t(L) and x
    ## if check_target_rank > 0, checks whether evident rank equals check_target_rank
    
    ## check symmetry
    d <- nrow(x)
    stopifnot( (ncol(x) == d) & (ma(x - t(x)) < stol))
    
    ## inefficient but want to be sure...
    ##
    ## evals <- eigen(x)$values
    ## evals <- evals[which(abs(evals) > ltol)]
    ## stopifnot(all(evals > 0))
    ## .rank <- length(evals)
    
    ## calculate cholesky
    ok <- c(1,rep(0,d-1))
    L <- matrix(0.0, ncol=d, nrow=d)
    
    L[1,1]     <- sqrt(x[1,1])
    L[(2:d),1] <- x[(2:d),1] / L[1,1]
    
    for(j in 2:d){
        .diag <- x[j,j] - sum(L[j,(1:(j-1))]^2)
        if(.diag < ltol){
            L[,j] <- 0.0
        } else {
            ok[j] <- 1
            L[j,j] <- sqrt(.diag)
            if(j < d){
                for(i in (j+1):d){
                    L[i,j] <- (1 / L[j,j]) * (x[i,j] - sum( L[i,1:(j-1)]*L[j,1:(j-1)]) )
                }
            }
        }
    }
    
    ## stopifnot(sum(ok) == .rank)
    L <- L[,which(ok > 0)]
    stopifnot( ma((L %*% t(L)) - x) < rtol)
    if(check_target_rank > 0){
        stopifnot(check_target_rank == sum(ok))
    }
    L
}
