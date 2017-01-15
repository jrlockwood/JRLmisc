rmultt <- function(n, df, mu, vmat, stol = 1e-07, ptol = 1e-07){
    ## Random sample of size n from MVT(df) distribution with location mu and scale vmat
    ## stol is a threshold used to test symmetry of vmat
    ## ptol is a threshold used to test positive definiteness of vmat
    ## Result is (n x length(mu)) matrix with iid rows
    p <- ncol(vmat)
    if(length(mu) != p)
        stop("Length of mu and dimension of vmat are not compatible")
    if(max(abs(vmat - t(vmat))) > stol)
        stop("vmat is not symmetric")
    if(any(eigen(vmat)$values < ptol))
        stop("vmat is not positive definite")
    if(df<=0)
        stop("df must be positive")
    vs <- svd(vmat)
    vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
    denom <- sqrt(rchisq(n, df)/df)
    spherical.mtx <- sweep(matrix(rnorm(n * p), nrow = n), 1, denom, "/")
    ans <- spherical.mtx %*% vsqrt
    return(sweep(ans, 2, mu, "+"))
}

ldmultt <- function(x, df, mu, vmat, stol = 1e-07, ptol = 1e-07){
    ## Calculates log density at x under MVT(df) distribution with location mu and scale vmat
    ## stol is a threshold used to test symmetry of vmat
    ## ptol is a threshold used to test positive definiteness of vmat
    p <- ncol(vmat)
    if(length(mu) != p)
        stop("Length of mu and dimension of vmat are not compatible")
    if(length(x) != p)
        stop("Length of mu and dimension of vmat are not compatible")
    if(max(abs(vmat - t(vmat))) > stol)
        stop("vmat is not symmetric")
    if(df<=0)
        stop("df must be positive")
    evals <- eigen(vmat)$values
    if(any(evals < ptol))
        stop("vmat is not positive definite")
    return(lgamma((df+p)/2)-lgamma(df/2)-((p/2)*log(df*pi))-0.5*sum(log(evals))-(((df+p)/2)*log(1+((1/df)*t(x-mu)%*%solve(vmat)%*%(x-mu)))))
}
