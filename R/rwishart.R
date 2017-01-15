rwishart <- function(df, b, stol = 1e-07, ptol = 1e-07){
    ## Random sample of size 1 from Wishart (df) distribution with scale matrix b
    ## stol is a threshold used to test symmetry of b
    ## ptol is a threshold used to test positive definiteness of b
    ## This is based on two results about Wishart random matrices:
    ##
    ## 1) If W is lower triangular with standard normals on the subdiagonals, and
    ##    the ith diagonal the square root of a chi-squared (df+1-i) variable
    ##    (all independent), then WW' ~ Wishart(df,I)
    ##
    ## 2) If W ~ Wishart(df,b) of dimension p and a is any (p x p) p.d. symmetric
    ##    matrix, then aWa' ~ Wishart(df,aba')
    ##
    ## Hence if W is constructed as in 1), and bsqrt is any symmetric square root
    ## of b, then (bsqrt%*%W)(bsqrt%*%W)' has the required distribution
    ##
    p <- ncol(b)
    if(max(abs(b - t(b))) > stol)
        stop("b is not symmetric")
    if(any(eigen(b)$values < ptol))
        stop("b is not positive definite")
    if(df < p)
        stop("Degrees of freedom must be at least as large as dimension")
    bs <- svd(b)
    bsqrt <- t(bs$v %*% (t(bs$u) * sqrt(bs$d)))
    W <- matrix(0,ncol=p,nrow=p)
    diag(W) <- sqrt(rchisq(n=p,(df+1-(1:p))))
    for(i in 2:p)
        W[i,1:(i-1)]<-rnorm(i-1)
    int<-bsqrt%*%W
    return(int%*%t(int))
}

ldwishart <- function(W, df, b, stol = 1e-07, ptol = 1e-07){
    ## Calculates log density at W under Wishart (df) distribution with scale matrix b
    ## stol is a threshold used to test symmetry of b and W
    ## ptol is a threshold used to test positive definiteness of b and W
    p.W <- ncol(W)
    p.b <- ncol(b)
    if(p.W != p.b)
        stop("Dimensions of W and b are not compatible")
    if(max(abs(W - t(W))) > stol)
        stop("W is not symmetric")
    if(max(abs(b - t(b))) > stol)
        stop("b is not symmetric")
    evals.b<-eigen(b)$values
    evals.W<-eigen(W)$values
    if(any(evals.W < ptol))
        stop("W is not positive definite")
    if(any(evals.b < ptol))
        stop("b is not positive definite")
    if(df<p.W)
        stop("Degrees of freedom too small")
    return(((df-p.W-1)/2)*sum(log(evals.W))-(df*p.W/2)*log(2)-(p.W*(p.W-1)/4)*log(pi)-(df/2)*sum(log(evals.b))-sum(lgamma((df+1-1:p.W)/2))-0.5*sum(diag(solve(b)%*%W)))
}
