cronbach <- function(x){
    ## calculates cronbach's alpha based on complete cases in
    ## (n x p) item data matrix x
    tmp <- as.matrix(na.omit(x))
    p <- dim(tmp)[2]
    (p/(p-1))*(1 - sum(apply(tmp,2,var))/var(apply(tmp,1,sum)))
}
