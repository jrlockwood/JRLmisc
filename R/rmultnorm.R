rmultnorm <- function(n, mu, vmat, stol = 1e-07, ptol = 1e-07){
    ## Random sample of size n from MVN(mu,vmat) distribution
    ## stol is a threshold used to test symmetry of vmat
    ## ptol is a threshold used to test positive definiteness of vmat
    ## Result is (n x length(mu)) matrix with iid rows
    p <- ncol(vmat)
    if(length(mu) != p) {
        stop("Length of mu and dimension of vmat are not compatible")
    }
    if(max(abs(vmat - t(vmat))) > stol)
        stop("vmat is not symmetric")
    if(any(eigen(vmat)$values < ptol))
        stop("vmat is not positive definite")
    vs <- svd(vmat)
    vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
    ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
    return(sweep(ans, 2, mu, "+"))
}

ldmultnorm <- function(x, mu, vmat, stol = 1e-07, ptol = 1e-07){
    ## Calculates log density at x under MVN(mu,vmat) distribution
    ## stol is a threshold used to test symmetry of vmat
    ## ptol is a threshold used to test positive definiteness of vmat
    p <- ncol(vmat)
    if(length(mu) != p) {
        stop("Length of mu and dimension of vmat are not compatible")
    }
    if(length(x) != p) {
        stop("Length of x and dimension of vmat are not compatible")
    }
    if(max(abs(vmat - t(vmat))) > stol)
        stop("vmat is not symmetric")
    evals <- eigen(vmat)$values
    if(any(evals < ptol))
        stop("vmat is not positive definite")
    return(-0.5*(p*log(2*pi)+sum(log(evals))+t(x-mu)%*%solve(vmat)%*%(x-mu)))
}
