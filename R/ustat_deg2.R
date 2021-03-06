ustat_deg2 <- function(x, psi){
    ## compute degree k=2 U-statistic based on (n x p) data matrix "x" and
    ## function "psi", along with analytical and asymptotic variance.
    ## NOTE: variance calculations are degenerate unless n >= 4.
    ## NOTE: x may have missing values, but in such case, psi() must be
    ## defined such that it returns a non-missing value. If missing values
    ## of psi result, the function stops
    
    if(!is.numeric(x)){
        stop("x must be numeric")
    }
    
    if(!is.matrix(x)){
        x <- matrix(x, ncol=1)
    }
    n <- nrow(x)
    
    if(n < 2){
        stop("there must be at least two observations")
    }
    
    if(!is.function(psi)){
        stop("psi must be a symmetric function accepting two arguments concordant with rows of x")
    }
    
    if(abs(psi(x[1,],x[2,]) - psi(x[2,],x[1,])) > 1e-8){
        stop("psi does not appear to be a symmetric function based on first two elements of x")
    }
    
    ij      <- t(combn(1L:n, 2))
    psivals <- apply(ij, 1, function(.ij){ psi(x[.ij[1],], x[.ij[2],]) })
    if(any(is.na(psivals))){
        stop("Evaluation of psi on pairs resulted in missing values")
    }
    
    ## compute U-statistic
    U       <- mean(psivals)
    
    ## compute terms needed to estimate exact and asymptotic variance.
    ##
    ## NOTE: the terms involving pairs with exactly one element in common are tricky to
    ## track.  we loop over choose(n,2) pairs and for each such target pair, we
    ## find the 2*(n-2) pairs that share exactly one element in common with it.
    ## thus there are a total of n(n-1)(n-2) products involved.
    
    ## "Epp"  denotes E[psi(X1,X2)*psi(X1,X3)] - product of pairs sharing one element.
    ## "Epsq" denotes E[psi(X1,X2)^2]
    Epp     <- 0.0
    ind     <- 0L
    for(istar in 1:(n-1)){
        for(jstar in (istar+1):n){
            ind <- ind + 1
            others <- setdiff(which( ((ij[,1] == istar) | (ij[,2] == istar)) | ((ij[,1] == jstar) | (ij[,2] == jstar)) ), ind)
            Epp <- Epp + sum(psivals[ind] * psivals[others])
        }
    }
    Epp  <- Epp / (n * (n-1) * (n-2))
    Epsq <- mean(psivals^2)
    
    ## get estimates of other pieces needed for variance calculation
    musq <- (1.0 / choose(n-2,2)) * ( (choose(n,2) * U^2) - Epsq - 2*(n-2)*Epp )
    ss1  <-  Epp - musq ## Cov(psi(X1,X2), psi(X1,X3)) - this drives asymptotic variance
    ss2  <- Epsq - musq ## Var(psi(X1,X2))
    
    return(list(psi         = psi,
                psivals     = psivals,
                U           = U,
                Epp         = Epp,
                Epsq        = Epsq,
                musq        = musq,
                ss1         = ss1,
                ss2         = ss2,
                varU.exact  = (1 / choose(n,2)) * ((2 * (n-2) * ss1) + ss2),
                varU.asymp  = 4*ss1 / n))
}
