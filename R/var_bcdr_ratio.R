var_bcdr_ratio <- function(.data, mu, delta, alpha, tau, mvars, pvars, soltol=1e-4){
  stopifnot(all(c(pvars,mvars) %in% names(.data)) && all(pvars == names(alpha)[-1]) && all(mvars == names(delta)[-1]))
 
  ## calculate psi, store as (n x K) matrix
  ## order of parameters is mu, delta, alpha, tau
  psi <- t(apply(as.matrix(.data), 1, function(u){
    wzm <- c(1, u[mvars]) ## add intercept
    xdelta <- sum(wzm * delta)
    wzp <- c(1, u[pvars]) ## add intercept
    xalpha <- sum(wzp * alpha)
    p    <- pcauchy(xalpha)
    pdot <- dcauchy(xalpha)
    psi1 <- (u["r"]*(u["y"] - xdelta) / p) + (tau * (xdelta - mu))
    psi2 <-  u["r"]*(u["y"] - xdelta) * wzm
    psi3 <- ( (u["r"] - p)*pdot/(p*(1-p))) * wzp
    psi4 <- u["r"]/p - tau
    return( c(psi1, psi2, psi3, psi4))
  }))

  ## check that (mu,delta,alpha,tau) solves estimating equations
  stopifnot(ma(apply(psi, 2, mean)) < soltol)

  ## get estimate of E[psi psi']
  B <- crossprod(psi) / nrow(.data)
  
  ## get estimate of mean derivative of psi w.r.t. params
  ## store as vector of column 1 corresponding to mu, column 2 is intercept of mean model, etc.
  ## in "d" blocks below, 1 = mu, 2 = delta, 3 = alpha, 4 = tau
  psidot   <- t(apply(as.matrix(.data), 1, function(u){
    wzm    <- c(1, u[mvars]) ## add intercept
    xdelta <- sum(wzm * delta)
    wzp    <- c(1, u[pvars]) ## add intercept
    xalpha <- sum(wzp * alpha)
    p      <- pcauchy(xalpha)
    pdot   <- dcauchy(xalpha)

    d11    <- -tau
    d21    <- rep(0, length(wzm))
    d31    <- rep(0, length(wzp))
    d41    <- 0
    
    d12    <- (tau - u["r"]/p) * wzm
    d22    <- -u["r"] * (wzm %*% t(wzm))
    d32    <- matrix(0, nrow=length(wzp), ncol=length(wzm))
    d42    <- rep(0, length(wzm))

    d13    <- (-u["r"]*(u["y"] - xdelta)*pdot / (p^2)) * wzp
    d23    <- matrix(0, nrow=length(wzm), ncol=length(wzp))
    d33    <- -( u["r"] * ( 1/p^2 + 2*pi*xalpha/p) + (1 - u["r"]) * ( 1 / (1-p)^2 - (2*pi*xalpha / (1-p)))) * (pdot^2) * ((wzp) %*% t(wzp))
    d43    <- (-u["r"]*pdot/ (p^2)) * wzp

    d14    <- xdelta - mu
    d24    <- rep(0, length(wzm))
    d34    <- rep(0, length(wzp))
    d44    <- -1

    as.vector( cbind(c(d11, d21, d31, d41), rbind(d12, d22, d32, d42), rbind(d13, d23, d33, d43), c(d14, d24, d34, d44)) )
  }))
  Ainv <- solve(matrix(apply(psidot, 2, mean), ncol=ncol(psi), byrow=F))

  ## return variance matrix
  res <- Ainv %*% B %*% t(Ainv) / nrow(.data)
  colnames(res) <- rownames(res) <- c("mu","m_intercept",paste("m",mvars,sep="_"),"p_intercept",paste("p",pvars,sep="_"),"tau")
  return(res)
}
