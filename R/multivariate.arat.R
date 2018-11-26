multivariate.arat <- function(thetahat, se, fixedmean = NULL, optim.reltol = NULL){
  ## thetahat is a (n x p) matrix of estimates and se a (n x p) matrix
  ## of standard errors.  missing values are allowed (as long as they align)
  ## and are removed from the MLE calculation
  ##
  ## we assume the following model for the estimates
  ## thetahat[i,j] = mu[j] + a[j]*(theta[i] + eta[i,j]) + e[i,j]
  ##
  ## where
  ## 1) theta[i] are iid normal with mean zero and unknown variance "tausq"
  ## 2) eta[i,j] are     normal with mean zero and AR covariance matrix (nusq/(1-rhosq))*(rho^|i-j|)
  ## 3) e[i,j]   are iid normal with mean zero and   known variance se[i,j]^2
  ## 4) a[1] = 1 and a[j] are (p-1) positive scaling constants
  ##
  ## this implies that thetahat[i,] are iid N(mu, G = A(tausq*J + AR(nusq,rho))A' + diag(se^2[i,]))
  ## where A is symmetric diag(1, a[2], a[3], ..)
  ##
  ## we estimate(tausq, nusq, rho, a, mu) by ML.
  ## if !is.null(fixedmean), uses the passed value as the known mean and only estimates variance components
  if(all(dim(thetahat)!=dim(se))){
    stop("est and se different dimensions")
  }
  if(min(dim(thetahat)) < 2){
    stop("doesn't work with only one column or row")
  }
  if(any(is.na(thetahat) != is.na(se))){
    stop("est and se must have same missing values")
  }
  v <- se^2
  np <- ncol(thetahat)
  stopifnot(np >= 3) ## else model is not identified

  ## create bad MOM starting values for tausq and nusq, start rho at 0.
  ## note we use psi = log( (1 - rho) / (1 + rho) ) like the KF stats in med paper.
  ## start "a" at 1
  m <- var(thetahat,use="p") - diag(apply(v,2,mean,na.rm=T))
  tausq0 <- mean(m[lower.tri(m)])
  nusq0  <- mean(diag(m)) - tausq0
  rho0   <- 0.0  
  mu0    <- apply(thetahat, 2, mean, na.rm=T)
  a0     <- rep(1, np-1)
  p0 <- c(ifelse(tausq0 < 0, 0, log(tausq0)), ifelse(nusq0 < 0, 0, log(nusq0)), log( (1 - rho0) / (1 + rho0) ), log(a0) )

  ## turn thetahat and v into lists of vectors, and create observation indicators
  thetahat <- as.list(as.data.frame(t(thetahat)))
  v        <- as.list(as.data.frame(t(v)))
  obs      <- lapply(thetahat, function(x){ !is.na(x) })

  ##  negative log likelihood, depending on whether fixedmean supplied.
  ##  we calculate G matrices, then log likelihood at each observation, and return negative sum.
  ##  note we deal with obs here, restricting each vector to its observed components.
  ##  note we had convergence problems here because optimizer was trying crazy values of
  ##  the parameters so we return "Inf" when parameters are out of range we normally want.
  if(is.null(fixedmean)){
    negll <- function(p){
      tausq <- exp(p[1])
      nusq  <- exp(p[2])
      rho <- (1 - exp(p[3])) / (1 + exp(p[3])) ## inverse transformation of psi
      at  <- exp(p[4:(4+np-2)])
      mu <- p[(4+np-1):(4+np-1+np-1)]
      if( (rho > 0.99) || (rho < -0.99) || (tausq < 1e-8) || (tausq > 1e8) || (nusq < 1e-8) || (nusq > 1e8) || (min(at) < 1e-8) || (max(at) > 1e3) ){
        return(Inf)
      } else {
        A <- diag(c(1,at))
        G <- lapply(v, function(x){ A %*% (matrix(tausq, ncol=np, nrow=np) + (nusq/(1 - rho^2))*AR(rho, np)) %*% A + diag(x) })
        G <- mapply(function(g, x){ g[x,x,drop=FALSE] }, G, obs, SIMPLIFY=FALSE)
        return( 0.5 * sum(mapply( function(g,x,o){ sum(log(eigen(g)$values)) + (t(x[o] - mu[o]) %*% solve(g) %*% (x[o] - mu[o])) }, G, thetahat, obs)))
      }
    }
    p0 <- c(p0, mu0)
  } else {
    negll <- function(p){
      tausq <- exp(p[1])
      nusq  <- exp(p[2])
      rho <- (1 - exp(p[3])) / (1 + exp(p[3]))
      at  <- exp(p[4:(4+np-2)])
      if( (rho > 0.99) || (rho < -0.99) || (tausq < 1e-8) || (tausq > 1e8) || (nusq < 1e-8) || (nusq > 1e8) || (min(at) < 1e-8) || (max(at) > 1e3) ){      
        return(Inf)
      } else {
        A <- diag(c(1,at))
        G <- lapply(v, function(x){ A %*% (matrix(tausq, ncol=np, nrow=np) + (nusq/(1 - rho^2))*AR(rho, np)) %*% A + diag(x) })        
        G <- mapply(function(g, x){ g[x,x,drop=FALSE] }, G, obs, SIMPLIFY=FALSE)
        return( 0.5 * sum(mapply( function(g,x,o){ sum(log(eigen(g)$values)) + (t(x[o] - fixedmean[o]) %*% solve(g) %*% (x[o] - fixedmean[o])) }, G, thetahat, obs)))
      }
    }
  }
  
  ## do the optimization
  .control <- list(maxit=500, trace=5, REPORT=1)
  if(!is.null(optim.reltol)){
    .control$reltol <- optim.reltol
  }
  out <- optim(par = p0, fn = negll, method  = "BFGS", control=.control, hessian=TRUE)
  stopifnot(out$convergence==0)
  tausq <- exp(out$par[1])
  nusq  <- exp(out$par[2])
  rho <- (1 - exp(out$par[3])) / (1 + exp(out$par[3]))
  at  <- exp(out$par[4:(4+np-2)])
  if(is.null(fixedmean)){
    mu <- out$par[(4+np-1):(4+np-1+np-1)]
  } else {
    mu <- fixedmean
  }
  ## Gj <- lapply(v, function(x){ matrix(tausq, ncol=np, nrow=np) + diag(x) + nusq*diag(np)})
  ## Gj <- mapply(function(gj,o){gj[o,o,drop=FALSE]}, Gj, obs, SIMPLIFY=FALSE)

  ## getting confidence intervals for tausq, nusq, and rho
  V <- solve(out$hessian)[1:3,1:3]
  tausqci <- exp(out$par[1] + (sqrt(V[1,1]) * c(-1.96, +1.96)))
  nusqci  <- exp(out$par[2] + (sqrt(V[2,2]) * c(-1.96, +1.96)))
  tmp <- rev(out$par[3] + (sqrt(V[3,3]) * c(-1.96, +1.96)))
  rhoci <- (1 - exp(tmp)) / (1 + exp(tmp))
  ## ratio uses delta method.  if psi1 = log(tausq) and psi2 = log(nusq) then
  ## ratio is 1 / (1 + exp(psi2 - psi1))
  ## psi <- out$par[1:2]
  ## del <- matrix(exp(psi[2] - psi[1]) * ( (1 + exp(psi[2] - psi[1]))^{-2} ) * c(1, -1), ncol=1)
  ## sdrat <- sqrt( t(del) %*% V %*% del )
  ## ratci <- rat + (sdrat * c(-1.96, +1.96))
  
  return(list(tausq = tausq, tausqci = tausqci, nusq = nusq, nusqci = nusqci, rho = rho, rhoci = rhoci, at = at, mu = mu, ll = -out$value))
}
