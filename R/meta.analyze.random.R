meta.analyze.random<-function(mi,sei,method="ML"){
    ## mi are individual estimates, sei are their standard errors
    ## originally written by R. Woodrow Setzer and modified by me.
    if(length(mi)!=length(sei)){
        stop("est and se different lengths")
    }
    if(any(is.na(mi)) | any(is.na(sei))){
        stop("missing values in est and/or se")
    }
    mll<-function (par, mi, vari){
        ## calculate -2 * log likelihood
        ## par[1] is the grand mean
        ## par[2] is the log of the between-group variance component
        sum(log(vari + exp(par[2]))) + sum((mi - par[1])^2/(vari + exp(par[2])))
    }
    mll.reml<-function (par, mi, vari){
        ## calculate -2 * log REML likelihood (see Lindstrom and Bates)
        sum(log(vari + exp(par[2]))) + sum((mi - par[1])^2/(vari +
                                                            exp(par[2]))) + log(sum(1/(vari + exp(par[2]))))
    }
    if(method!="ML" & method!="REML"){
        stop("Invalid value of method")
    }
    if (length(mi) == 1) {
        c(muhat = mi, muhat.se = sei, sdhat = NA)
    }
    else {
        vari <- sei^2
        start <- c(mean(mi), log(var(mi)))
        objfun <- if(method == "ML") mll else mll.reml
        out <- try(optim(par = start, fn = objfun, method
                         = "BFGS",control=list(maxit=500),mi = mi, vari = vari))
        
        if (inherits(out, "try-error")) {
            print(cbind(mi, sei))
            print(out)
            stop("Problem in optim")
        }
        if (out$convergence == 0)
            c(muhat = out$par[1], muhat.se = sqrt(1/sum(1/(vari +
                                                           exp(out$par[2])))), sdhat = sqrt(exp(out$par[2])))
        else {
            c(muhat = NA, muhat.se = NA, sdhat = NA)
        }
    }
}
