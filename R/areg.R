areg <- function(d, y, xit, idvar){
    ## fits fixed effects model y_it = u_i + beta'xit + eit by
    ## absorption of fixed effects.  right now this is very limited
    ## and just returns the coefficients and the estimated fixed
    ## effects.
    ##
    ## d = dataframe (no missing values allowed)
    ## y = name of response variable
    ## xit = name of time-varying variables.  must all be numeric - factors
    ##       need to be expanded in advance
    ## idvar = name of id variable requiring fixed effects
    if( !all( c(y, xit, idvar) %in% names(d) ) ){
        stop("not all variables in d")
    }
    if( nrow(na.omit(d)) < nrow(d) ){
        stop("no missing values allowed")
    }
    if( !all(sapply(d[,xit], is.numeric)) ){
        stop("all predictors must be numeric")
    }
    ## since collection of fixed effects are coincident with the intercept
    ## make sure that model is identified
    if( qr(cbind(1,as.matrix(d[,xit])))$rank < (length(xit) + 1) ){
        stop("model is not of full rank")
    }
    ## calculate matrix of deviations (Y_it - Y_ibar) and (X_it - X_ibar)
    vars <- c(y,xit)
    gmcent <- apply( as.matrix(d[,vars]), 2, function(x){x - ave(x, d[, idvar])})
    if( any(apply(gmcent, 2, var) < 1e-7) ){
        stop("predictor or response not time-varying")
    }
    ## get betahat
    bhat <- coef(lm(gmcent[,y] ~ gmcent[,xit] - 1))
    names(bhat) <- xit
    ## recover estimates of fixed effects
    gm <- aggregate(d[,vars], list(d[,idvar]), mean)
    names(gm)[1] <- idvar
    gm$u <- as.vector(gm[,y] - (as.matrix(gm[,xit]) %*% bhat))
    u.named <- gm$u
    names(u.named) <- as.character(gm[,idvar])
    
    return(list(bhat = bhat, u = gm[,c(idvar,"u")], coef = c(u.named, bhat)))
}
