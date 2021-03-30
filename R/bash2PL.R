#############################
## slow likelihood basher for 2PL
#############################
bash2PL <- function(X, mu, sigma, nquad = 81, ...){
    p <- ncol(X)
    tgrid = seq(from = mu - 4*sigma, to = mu + 4*sigma, length=nquad)
    stopifnot(all(na.omit(c(X)) %in% c(0L,1L)))
    ptheta <- dnorm(tgrid, mu, sigma)
    negll <- function(psi){
        a     <- exp(psi[1:p])
        b     <- psi[(p+1):(2*p)]
        logp  <- t(sapply(tgrid, function(.t){ log( (1 / (1 + exp(-a*(.t - b)))) )}))
        log_1_minus_p <- t(sapply(tgrid, function(.t){ log( 1 - (1 / (1 + exp(-a*(.t - b))))) }))
        -1.0 * sum(apply(X, 1, function(r){
            wh <- which(!is.na(r))
            .l <- exp( (logp[,wh,drop=F] %*% r[wh]) + (log_1_minus_p[,wh,drop=F] %*% (1 - r[wh])) )
            log(sum(.l * ptheta))
        }))
    }
    res <- optim(rep(0,2*p), negll, method="BFGS", control=list(maxit=2000), ...)
    return(cbind(a = exp(res$par[1:p]), b = res$par[(p+1):(2*p)]))
}

