rdirmult <- function(n = 1, size, alpha){
  ## Random sample of size n from non-degenerate Dirichlet-Multinomial (size, alpha) distribution
  ## Requires rdirichlet.fun() and rmultinomial.fun()
  ## Result is (n x length(alpha)) matrix with iid rows
  pmat<-rdirichlet(n,alpha)
  return(t(apply(pmat,1,rmultinomial,n=1,size=size)))
}

lddirmult <- function(x, alpha){
    ## Calculates log "density" at x under non-degenerate Dirichlet-multinomial (sum(x), alpha) distribution
    if(any(alpha<=0))
        stop("Invalid parameter")
    if((length(x)!=length(alpha)))
        stop("x and alpha dimensions do not match")
    size<-sum(x)
    adot<-sum(alpha)
    return(lgamma(size+1)-sum(lgamma(x+1))+lgamma(adot)-lgamma(size+adot)+sum(lgamma(x+alpha)-lgamma(alpha)))
}
