blp <- function(Y, lambda, mu, SigmaX, SigmaU, etol=1e-6){
    ## X has E[X] = mu, var[X] = SigmaX
    ## Y = X + U where E[U|X] = 0, var[U] = SigmaU
    ## compute BLP of lambda'X from Y, along with MSE
    ##
    ## algebra:
    ## Q = SigmaX(SigmaX + SigmaU)^{-1}
    ## BLP = lambda'[ (I-Q)mu + QY ]
    ## MSE = tr(lambda lambda' (Q SigmaU Q' + (I-Q)SigmaX(I-Q)'))

    ## checks
    stopifnot(is.numeric(Y) && is.numeric(lambda) && is.numeric(mu) && is.numeric(SigmaX) && is.numeric(SigmaU) && is.numeric(etol))
    stopifnot(all(!is.na(c(Y,lambda,mu,SigmaX,SigmaU,etol))))
    p <- length(Y)
    if(p==1){
        dim(SigmaX) <- dim(SigmaU) <- c(1,1)
    }    
    stopifnot(is.matrix(SigmaX) && is.matrix(SigmaU))
    stopifnot( (length(lambda) == p) && (length(mu) == p) && (nrow(SigmaX) == p) && (ncol(SigmaX) == p) && (nrow(SigmaU) == p) && (ncol(SigmaU) == p) )
    stopifnot(max(abs(SigmaX - t(SigmaX))) < 1e-8)
    stopifnot(max(abs(SigmaU - t(SigmaU))) < 1e-8)
    dim(Y) <- dim(lambda) <- dim(mu) <-  c(p,1)
    stopifnot(all(eigen(SigmaX + SigmaU)$values > etol))
    
    ## computations
    est.direct <- as.vector(t(lambda) %*% Y)
    mse.direct <- as.vector(t(lambda) %*% SigmaU %*% lambda)
    Q          <- SigmaX %*% solve(SigmaX + SigmaU)
    IminusQ    <- diag(p) - Q
    est.blp    <- as.vector(t(lambda) %*% ( (IminusQ %*% mu) + (Q %*% Y) ))
    part       <- (Q %*% SigmaU %*% t(Q)) + (IminusQ %*% SigmaX %*% t(IminusQ))
    mse.blp    <- sum(diag( (lambda %*% t(lambda)) %*% part))

    return(c(est.direct = est.direct, mse.direct = mse.direct, est.blp = est.blp, mse.blp = mse.blp, prmse.blp = (1 - mse.blp/mse.direct)))
}
    