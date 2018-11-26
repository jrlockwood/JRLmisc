multivariate.vc.general <- function(d, obs, Y, X, vE, make.vT, kV, kB, thetaV, thetaB = NULL, det.tol = 1e-10, do.hessian = TRUE, do.blup = FALSE, optim.reltol = NULL){  
  ## estimate parameters "thetaV" and "thetaB" by maximum likelihood given known error variance components.
  ## this subsumes univariate.vc and multivariate.* variance component functions, including ability to estimate FH model
  ##
  ## MODEL
  ## Y = XB + T + E
  ## Y (d x 1) vector
  ## E (d x 1) ~ N(0, vE), vE is observation-specific known covariance matrix
  ## T (d x 1) ~ N(0, vT), vT structured covariance in kV parameters "thetaV"
  ## X (d x kB) design matrix for mean structure with parameters B = "thetaB"
  ##
  ## INPUTS:
  ## d - full dimension of problem
  ## obs - list of length n with vectors of observed indices - i.e. which of 1:d is observed, in order
  ## Y - list of length n with observed values, each a vector of length length(obs[[i]])
  ## X - list of design matrix elements with dimension (length(obs[[i]]) * kB)
  ## vE  - list of known error covariance matrices, of dimensions length(obs[[i]])*length(obs[[i]])
  ## make.vT: function to construct fully expanded vT as a function of thetaV
  ## kV - length of "thetaV"
  ## kB - length of "thetaB". NOTE: if kB == 0, we are assuming fixed means provided in the values of X
  ## thetaV starting values (required)
  ## thetaB starting values (optional)                                  
  ## det.tol: smallest determinant we allow for matrices
  ## do.hessian: TRUE/FALSE
  ## do.blup: TRUE/FALSE
  ## optim.reltol: relative convergence tolerance for optim
  n <- length(Y)
  .dims <- sapply(Y, length)
  stopifnot( (length(obs) == n) && (length(X) == n) && (length(vE) == n) )
  stopifnot(max(.dims) <= d)
  stopifnot( all(sapply(obs, is.vector)) && all(sapply(X, is.matrix)) && all(sapply(vE, is.matrix)) && all(sapply(Y, is.vector)) )
  stopifnot( all(sapply(obs, length) == .dims) && all(sapply(X, nrow) == .dims) && all(sapply(vE, nrow) == .dims) && all(sapply(vE, ncol) == .dims) )
  stopifnot(all(!sapply(obs, function(x){ any(duplicated(x)) })))
  stopifnot( all(!is.na(unlist(Y))) && all(!is.na(unlist(obs))) && all(!is.na(unlist(X))) && all(!is.na(unlist(vE))) )
  stopifnot( is.function(make.vT) )
  if(kB == 0){ stopifnot(all(sapply(X, ncol) == 1)) }
  stopifnot( length(thetaV) == kV )

  ## if we are estimating means and they are not provided, set decent starting values with regression
  if( (kB > 0) && is.null(thetaB)){
    .y <- unlist(Y)
    .x <- do.call("rbind", X)
    thetaB <- as.vector(coef(lm(.y ~ .x - 1)))
    rm(.y, .x)
  }

  ##  negative log likelihood, depending on whether we are estimating mean or not.
  ##  we calculate variance matrices, then log likelihood at each observation, and return negative sum.
  ##  note we deal with obs here, restricting each vector to its observed components.
  ##  return "Inf" when parameters are out of range we normally want.
  ##  use determinant() function rather than log eigen values because it was faster
  
  if(kB > 0){  ## we are estimating mean parameters
    negll <- function(theta){
      thetaV <- theta[1:kV]
      thetaB <- theta[(kV+1):(kV+kB)]
      vT <- make.vT(thetaV)
      G  <- mapply(function(ve, o){ vT[o,o,drop=FALSE] + ve }, vE, obs, SIMPLIFY=FALSE)
      d2 <- lapply(G, determinant, logarithm=TRUE)
      d2m <- lapply(d2, function(x){ x$modulus })
      if( any(sapply(d2, function(x){ x$sign }) < 0) || any(unlist(d2m) < log(det.tol))){
        return(Inf)
      } else {
        M <- lapply(X, function(x){ x %*% thetaB })
        ## SLOWER: return( 0.5 * sum(mapply( function(y,g,m){ sum(log(eigen(g)$values)) + (t(y-m) %*% solve(g) %*% (y-m)) }, Y, G, M)))
        return( 0.5 * sum(mapply( function(y,g,m,ld){ ld + (t(y-m) %*% solve(g) %*% (y-m)) }, Y, G, M, d2m)))
      }
    }
    theta0 <- c(thetaV, thetaB)
  } else { ## fixed means in column of X
    negll <- function(theta){
      vT <- make.vT(theta)
      G  <- mapply(function(ve, o){ vT[o,o,drop=FALSE] + ve }, vE, obs, SIMPLIFY=FALSE)
      d2 <- lapply(G, determinant, logarithm=TRUE)
      d2m <- lapply(d2, function(x){ x$modulus })
      if( any(sapply(d2, function(x){ x$sign }) < 0) || any(unlist(d2m) < log(det.tol))){
        return(Inf)
      } else {
        ## SLOWER: return( 0.5 * sum(mapply( function(y,g,m){ sum(log(eigen(g)$values)) + (t(y-m) %*% solve(g) %*% (y-m)) }, Y, G, X)))
        ## FLAG: why do i call determinant again below?  this seems like a mistake, should just use "ld"
        return( 0.5 * sum(mapply( function(y,g,m,ld){ determinant(g, logarithm=TRUE)$modulus + (t(y-m) %*% solve(g) %*% (y-m)) }, Y, G, X, d2m)))
      }
    }
    theta0 <- thetaV
  }
  
  ## do the optimization
  .control <- list(maxit=1000, trace=5, REPORT=1)
  if(!is.null(optim.reltol)){
    .control$reltol <- optim.reltol
  }
  out <- optim(par = theta0, fn = negll, method  = "BFGS", control=.control, hessian=do.hessian)
  stopifnot(out$convergence==0)
  thetaV.hat <- out$par[1:kV]
  vT.hat <- make.vT(thetaV.hat)
  if(kB > 0){
    thetaB.hat <- out$par[(kV+1):(kV+kB)]
  } else {
    thetaB.hat <- NULL
  }
  if(do.hessian){
    theta.hessian <- out$hessian
  } else {
    theta.hessian <- NULL
  }

  ## calculate BLUPS of (XB + T)
  ## XB + vT(vT+vE)^{-1}(Y - XB)
  if(do.blup){
    if(kB > 0){
      M <- lapply(X, function(x){ x %*% thetaB.hat })
    } else {
      M <- X
    }
    blup <- mapply(function(o,y,m,ve){ m + (vT.hat[o,o,drop=FALSE] %*% solve(vT.hat[o,o,drop=FALSE] + ve)) %*% (y - m) }, obs, Y, M, vE, SIMPLIFY = FALSE)
  } else {
    blup <- NULL
  }

  ## return
  return(list(thetaV = thetaV.hat, thetaB = thetaB.hat, vT = vT.hat, theta.hessian = theta.hessian, ll = -out$value, blup = blup))
}
