summarize <- function(x){
    ## Calculate summary statistics for a vector
    if( (!is.vector(x)) & (is.factor(x)) ){
        v<-c(length(x),nobs(x),NA,NA,NA,NA)
    }
    if( (is.factor(x))|(is.character(x)) ){
        v<-c(length(x),nobs(x),NA,NA,NA,NA)
    } else{
        v<-c(length(x),nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),mean(x,na.rm=TRUE),sd(x,na.rm=TRUE))
    }
    names(v)<-c("N","nobs","min","max","mean","sd")
    v
}
