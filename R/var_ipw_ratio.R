var_ipw_ratio <- function(.data, mu, alpha, pvars, soltol=1e-4){
  stopifnot(all(pvars %in% names(.data)) && all(pvars == names(alpha)[-1]))
  
  ## calculate psi, store as (n x (p+1)) matrix
  psi <- t(apply(as.matrix(.data), 1, function(u){
    wz <- c(1, u[pvars]) ## add intercept
    xalpha <- sum(wz * alpha)
    p    <- pcauchy(xalpha)
    pdot <- dcauchy(xalpha)
    return( c(u["r"]*(u["y"] - mu) / p, ((u["r"] - p)*pdot / (p * (1-p)))*wz))
  }))

  ## check that (mu,alpha) solves estimating equations
  stopifnot(ma(apply(psi, 2, mean)) < soltol)

  ## get estimate of E[psi psi']
  B <- crossprod(psi) / nrow(.data)

  ## get estimate of mean derivative of psi w.r.t. params
  ## store as vector of column 1 corresponding to mu, column 2 is intercept of GLM, etc.
  psidot   <- t(apply(as.matrix(.data), 1, function(u){
    wz     <- c(1, u[pvars])
    xalpha <- sum(wz * alpha)
    p      <- pcauchy(xalpha)
    pdot   <- dcauchy(xalpha)
    d11    <- -u["r"]/p
    d21    <- rep(0, length(wz))
    d12    <- (-u["r"] * (u["y"] - mu) / p^2) * pdot * wz
    d22    <- -( u["r"] * ( 1/p^2 + 2*pi*xalpha/p) + (1 - u["r"]) * ( 1 / (1-p)^2 - (2*pi*xalpha / (1-p)))) * (pdot^2) * ((wz) %*% t(wz))
    as.vector(cbind(c(d11, d21), rbind(d12, d22)))
  }))
  Ainv <- solve(matrix(apply(psidot, 2, mean), ncol=ncol(psi), byrow=F))

  ## return variance matrix
  res <- Ainv %*% B %*% t(Ainv) / nrow(.data)
  colnames(res) <- rownames(res) <- c("mu","intercept",pvars)
  return(res)
}
