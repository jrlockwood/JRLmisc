rmultinomial <- function(n = 1, size, p, ptol = 1e-07){
    ## Random sample of size n from non-degenerate multinomial (size, p) distribution
    ## ptol is a threshold used to ensure that sum of p is sufficiently close to 1
    ## Result is (n x length(p)) matrix with iid rows
    if((any(p<=0))|(abs(sum(p)-1.0)>ptol))
        stop("Invalid probability vector")
    ncells<-length(p)
    mat<-matrix(sample(1:ncells, size = n * size, prob = p, replace = T),nrow=n)
    return(t(apply(mat,1,tabulate,nbins=ncells)))
}

ldmultinomial <- function(x, p, ptol = 1e-07){
    ## Calculates log "density" at x under non-degenerate multinomial (sum(x), p) distribution
    ## ptol is a threshold used to ensure that sum of p is sufficiently close to 1
    if((any(p<=0))|(abs(sum(p)-1.0)>ptol))
        stop("Invalid probability vector")
    if((length(x)!=length(p)))
        stop("x and p dimensions do not match")
    return(lgamma(sum(x)+1)-sum(lgamma(x+1))+sum(x*log(p)))
}
