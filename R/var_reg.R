var_reg <- function(.data, mu, delta, mvars, soltol=1e-4){
  stopifnot(all(mvars %in% names(.data)) && all(mvars == names(delta)[-1]))
  
  ## calculate psi, store as (n x (p+1)) matrix
  psi <- t(apply(as.matrix(.data), 1, function(u){
    wz <- c(1, u[mvars]) ## add intercept
    xdelta <- sum(wz * delta)
    return( c(xdelta - mu, u["r"]*(u["y"] - xdelta)*wz))
  }))

  ## check that (mu,delta) solves estimating equations
  stopifnot(ma(apply(psi, 2, mean)) < soltol)

  ## get estimate of E[psi psi']
  B <- crossprod(psi) / nrow(.data)

  ## get estimate of mean derivative of psi w.r.t. params
  ## (first column is mu, second column is intercept of regression, etc)
  wz <- as.matrix(cbind(1, .data[,mvars]))
  Ainv <- solve( rbind(c(-1, apply(wz, 2, mean)), cbind(0, -t(.data$r * wz) %*% (.data$r * wz) / nrow(.data))))

  ## return variance matrix
  res <- Ainv %*% B %*% t(Ainv) / nrow(.data)
  colnames(res) <- rownames(res) <- c("mu","intercept",mvars)
  return(res)
}
