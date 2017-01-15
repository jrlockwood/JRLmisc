rdirichlet <- function(n = 1, alpha){
    ## Random sample of size n from non-degenerate Dirichlet (alpha) distribution
    ## Result is (n x length(alpha)) matrix with iid rows
    if(any(alpha<=0))
        stop("Invalid parameter")
    gmat <- matrix(rgamma(n * length(alpha), alpha), nrow = n, byrow = T)
    rowsums <- apply(gmat, 1, sum)
    return(gmat/rowsums)
}

lddirichlet <- function(x, alpha){
    ## Calculates log density at x under non-degenerate Dirichlet (alpha) distribution
    if(any(alpha<=0))
        stop("Invalid parameter")
    if((length(x)!=length(alpha)))
        stop("x and alpha dimensions do not match")
    return(lgamma(sum(alpha))-sum(lgamma(alpha))+sum((alpha-1)*log(x)))
}
