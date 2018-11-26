univariate.vc <- function(thetahat, se, fixedmean = NULL, X = NULL, optim.reltol = NULL){
  #########################################################################  
  ## thetahat are estimates with standard errors se.
  ## estimates variance component tausq of theta under the model
  ##
  ## theta[i]    ~ N(mu, tausq)
  ## thetahat[i] ~ N(theta[i], se[i]^2)
  ##
  ## uses ML to estimate.  see meta.analyze.random for related function
  ##
  ## possible structures for "mu":
  ## !is.null(fixedmean) & !is.null(X)       error
  ## !is.null(fixedmean) &  is.null(X)       estimates tausq assuming fixed mean
  ##  is.null(fixedmean) &  is.null(X)       estimates intercept only
  ##  is.null(fixedmean) & !is.null(X)       estimates coeffs for all columns of X
  ##                                         NOTE: if intercept is desired it must be included in X
  #########################################################################
  stopifnot( (length(thetahat) > 1) && (length(thetahat) == length(se)) )
  stopifnot( all(!is.na(thetahat)) && all(!is.na(se)) )
  if(!is.null(X)){  stopifnot( is.matrix(X) && is.null(fixedmean) ) }
  if(!is.null(fixedmean)){ stopifnot( !is.matrix(fixedmean) && (length(fixedmean)==1) ) }
  v <- se^2

  ## negative log likelihood, depending on mean structure
  if( is.null(fixedmean) & is.null(X) ){
    negll <- function(par){
      0.5 * ( sum(log(v + exp(par[1]))) + sum((thetahat - par[2])^2/(v + exp(par[1]))) )
    }
    par0 <- c(log(var(thetahat)), mean(thetahat))
  } else if( !is.null(fixedmean) & is.null(X) ){
    negll <- function(par){
      0.5 * ( sum(log(v + exp(par))) + sum((thetahat - fixedmean)^2/(v + exp(par))) )
    }
    par0 <- log(var(thetahat))
  } else if( is.null(fixedmean) & !is.null(X) ){
    negll <- function(par){
      0.5 * ( sum(log(v + exp(par[1]))) + sum((thetahat - (X %*% par[-1]) )^2/(v + exp(par[1]))) )
    }
    par0 <- c(log(var(thetahat)), as.vector(coef(lm(thetahat ~ X - 1))))
  }
  
  ## do the optimization
  .control <- list(maxit=500, trace=5, REPORT=1)
  if(!is.null(optim.reltol)){
    .control$reltol <- optim.reltol
  }
  out <- optim(par = par0, fn = negll, method  = "BFGS", control=.control)
  stopifnot(out$convergence==0)
  tausq <- exp(out$par[1])
  
  if( is.null(fixedmean) & is.null(X) ){
    beta <- numeric(0)
    mu   <- out$par[2]
  } else if( !is.null(fixedmean) & is.null(X) ){
    beta <- numeric(0)
    mu   <- fixedmean
  } else if( is.null(fixedmean) & !is.null(X) ){
    beta <- out$par[-1]
    mu   <- as.vector(X %*% beta)
  }
  blup <- mu + (tausq/(tausq + v) * (thetahat - mu))
  
  return(list(tausq = tausq, beta = beta, mu = mu, blup = blup, ll = -out$value))
}
