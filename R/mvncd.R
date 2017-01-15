mvncd <- function(mu, vmat, x, indicator, stol = 1e-07, ptol = 1e-07){
    ## MVN conditional distribution
    ##
    ## INPUTS:
    ## mu: mean vector
    ## vmat: covariance matrix
    ## x: vector (of same length of mu) which contains observed coordinates.
    ##    Only x[indicator>0] are used
    ## indicator: vector (of same length of mu) which contains positive values
    ##            for those coordinates that are being conditioned on and
    ##            non-positive values for the coordinates whose conditional
    ##            distribution is to be calculated.
    ## stol: a threshold used to test symmetry of vmat
    ## ptol: a threshold used to test positive definiteness of vmat
    ## 
    ## OUTPUTS: (a list)
    ##
    ## logmargdens: log marginal density of the observed coordinates
    ## nO: number of observed coordinates
    ## nU: number of unobserved coordinates
    ## cmu: conditional mean vector
    ## cvar: conditional covariance matrix
    ##
    p <- ncol(vmat)
    if(length(mu) != p) {
        stop("Length of mu and dimension of vmat are not compatible")
    }
    if(length(x) != p) {
        stop("Length of x and dimension of vmat are not compatible")
    }
    if(length(indicator) != p) {
        stop("Length of indicator and dimension of vmat are not compatible")
    }
    if(max(abs(vmat - t(vmat))) > stol)
        stop("vmat is not symmetric")
    if(any(eigen(vmat)$values < ptol))
        stop("vmat is not positive definite")
    if((all(indicator>0))|(all(indicator<=0)))
        stop("indicator implies degenerate conditioning")
    obs.ind<-indicator>0
    unobs.ind<-indicator<=0
    nO<-sum(obs.ind)
    nU<-sum(unobs.ind)
    muO<-mu[obs.ind]
    muU<-mu[unobs.ind]
    xmmuO<-x[obs.ind]-muO
    SigOO<-vmat[obs.ind,obs.ind]
    SigUU<-vmat[unobs.ind,unobs.ind]
    if(nU==1){
        SigUO<-matrix(vmat[unobs.ind,obs.ind],nrow=1)
    } else {
        SigUO<-vmat[unobs.ind,obs.ind]
    }
    SigOOinv<-solve(SigOO)    
    cmu<-muU+(SigUO%*%SigOOinv%*%xmmuO)
    reg<-SigUO%*%SigOOinv
    cvar<-SigUU-(SigUO%*%SigOOinv%*%t(SigUO))
    logmargdens<- -0.5*(nO*log(2*pi)+sum(log(eigen(SigOO)$values))+t(xmmuO)%*%SigOOinv%*%xmmuO)
    return(list(logmargdens=logmargdens,nO=nO,nU=nU,cmu=cmu,cvar=cvar,mvar=SigOO,reg=reg))
}
