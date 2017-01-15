meta.analyze.fixed<-function(est,se){
    ## given a vector of estimated effects and standard associated
    ## standard errors, calculates a precison-weighted overall
    ## average and its associated standard error.
    if(length(est)!=length(se)){
        stop("est and se different lengths")
    }
    if(any(is.na(est)) | any(is.na(se))){
        stop("missing values in est and/or se")
    }
    if(length(est)==1){
        return(c(est,se))
    } else{
        p<-se^(-2)
        w<-p/sum(p)
        ma.est<-sum(w*est)
        ma.se<-sqrt(1/sum(p))
        return(c(ma.est,ma.se))
    }
}
