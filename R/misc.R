l <- function(x){
    length(x)
}

lu <- function(x){
    length(unique(x))
}

lsu <- function(x){
    length(unique(sort(x)))
}

su <- function(x){
    sort(unique(x))
}

st <- function(x){
    sort(table(x, exclude=NULL))
}

ma <- function(x){
    max(abs(x),na.rm=TRUE)
}

nmis <- function(x){
    sum(is.na(x))
}

nobs <- function(x){
    sum(!is.na(x))
}

pmis <- function(x){
    sum(is.na(x))/length(x)
}

pobs <- function(x){
  sum(!is.na(x))/length(x)
}

csign <- function(x){
    res<-rep("",length(x))
    res[x<0]<-"-"
    res[x>0]<-"+"
    res
}

logit <- function(p){
    log(p/(1-p))
}

expit <- function(x){
    v<-exp(x)
    v/(1+v)
}

cs <- function(r,d){
    ## compound symmetry correlation matrix
    matrix(r,ncol=d,nrow=d)+(1-r)*diag(d)
}

ar <- function(phi,d){
    ## AR(1) correlation matrix
    matrix(apply(expand.grid(1:d,1:d),1,function(x){phi^{abs(x[2]-x[1])}}),nrow=d,ncol=d,byrow=F)
}

J <- function(dim){
    matrix(1,ncol=dim,nrow=dim)
}

more <- function(x) {
    ## taken from R-help
    tmp <- paste("/tmp/__R", floor(runif(1,0,1e6)), sep=".")
    sink(tmp)
    print(x)
    sink()
    file.show(tmp, delete.file=TRUE)
}

obsize<-function(){
    as.matrix(sort(sapply(objects(name=".GlobalEnv"), function(x){ object.size(get(x)) })))
}

logmeanexp <- function(x){
    ## Calculates log(mean(exp(x))) where elements of x might be very small
    m<-max(x)
    log(mean(exp(x-m)))+m
}

enclose<-function(x,rng){
    ## rng is the max and min values you want x to take.
    ## function sets all x: x<rng[1] to rng[1] and
    ## x: x>rng[2] to rng[2]
    v<-x
    v[v<rng[1]]<-rng[1]
    v[v>rng[2]]<-rng[2]
    v
}

scatterpoints<-function(x,y,f.lowess=1/3,pch=".",...){
    if(length(x)!=length(y)){
        stop("scatterpoints(): lengths of x and y differ")
    }
    ntot<-length(x)
    obs<-!(is.na(x)|is.na(y))
    n<-sum(obs)
    x.obs<-x[obs]
    y.obs<-y[obs]
    points(x.obs,y.obs,pch=pch)
    abline(lm(y.obs~x.obs),col=124,lwd=2)
    lines(lowess(y.obs~x.obs,f=f.lowess),col=33,lwd=2)
    ttle<-paste(n,"/",ntot,"; r=",round(cor(x.obs,y.obs),digits=2),sep="")
    title(main=ttle)
}

scatterplot<-function(x,y,xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),f.lowess=1/3,pch=19,cex=0.5,col="gray",...,do.lowess=TRUE){
    if(length(x)!=length(y)){
        stop("scatterplot(): lengths of x and y differ")
    }
    ntot<-length(x)
    obs<-!(is.na(x)|is.na(y))
    n<-sum(obs)
    x.obs<-x[obs]
    y.obs<-y[obs]
    plot(x.obs,y.obs,xlab=xlab,ylab=ylab,pch=pch,cex=cex,col=col,...)
    abline(lm(y.obs~x.obs),col=124,lwd=2)
    if(do.lowess){
        lines(lowess(y.obs~x.obs,f=f.lowess),col=33,lwd=2)
    }
    ttle<-paste(n,"/",ntot,"; r=",round(cor(x.obs,y.obs),digits=2),sep="")
    title(main=ttle)
}

cfill<-function(x,yl,yu,color="lightgray"){
    ## fill the area between two curves with "color" and
    ## assumes plot.new has been called with the right axes
    n<-length(x)
    if( (length(yl)!=n)|(length(yu)!=n) ){
        stop("cfill: Incorrect lengths")
    }
    if(!all(yl<yu)){
        stop("cfill: Lower curve not below upper curve")
    }
    for(i in 1:(n-1)){
        polygon(x=c(x[i],x[i],x[i+1],x[i+1],x[i]),y=c(yl[i],yu[i],yu[i+1],yl[i+1],yl[i]),col=color,border=NA)
    }
}

isgap <- function(x){
    ## tests for gaps in a vector, i.e. is there a missing value sandwiched between observed values
    stopifnot(is.vector(x))
    if( all(is.na(x)) || all(!is.na(x)) ){
        return(0)
    } else{
        r <- range(which(!is.na(x)))
        return( ifelse( any(is.na(x[r[1]:r[2]])), 1, 0) )
    }
}

read.csvnogarb <- function(x,...){
    ## x is path to csv file.  gets rid of garbage chars before reading
    ## ...: additional arguments to read.csv
    system(paste("cat ",x,"| tr -d \"#$&'()*:/~\"  > /tmp/__.csv"))
    read.csv("/tmp/__.csv",...)
}

rbzs <- function(x){
    ## calculates rank-based z-score of vector x (missing values allowed)
    qnorm(rank(x, na.last="keep") / (nobs(x) + 1))
}

na.show <- function(x){
    ## returns "mirror image" of na.omit(x) for a dataframe - i.e.
    ## a dataframe with only the records that were dropped by na.omit
    if(!is.data.frame(x)){
        stop("error in na.show: x must be a dataframe")
    }
    x[as.vector(attributes(na.omit(x))$na.action),]
}

put <- function(ofile, string){
    sink(ofile, append=TRUE)
    cat(string)
    sink(NULL)
}

descore <- function(x){
    ## x is character vector.  returns x with underscores removed (useful for parsing SAS names)
    gsub("_","",x)
}

remap <- function(x,orig,repl){
    ## takes the vector x and replaces each value in "orig" with the corresponding
    ## values in "repl".  missing values are not allowed
    if( any(is.na(x)) | any(is.na(orig)) | any(is.na(repl)) | (length(orig)!=length(repl)) ){
        stop("error in remap")
    }
    xnew <- x
    for(i in 1:length(orig)){
        xnew[which(x==orig[i])] <- repl[i]
    }
    xnew
}

chartoform <- function(respname,xname){
    as.formula(paste(respname," ~ ",paste(xname,collapse=" + ")))
}

freq <- function(x){
    ## taken from R-help
    xmat<-as.matrix(x)
    ifelse (ncol(xmat)==1,{
        Count<-table(x)
        Total<-sum(Count)
        Prcnt<-100*(Count/Total)
        x1<-cbind(Count,Prcnt)
        x2<-cbind(Total,sum(Prcnt))
        Frequency.Table<-as.data.frame(rbind(x1,x2))
        c<-nrow(Frequency.Table)
        rownames(Frequency.Table)[c]<-"Total"
        return(Frequency.Table)},return("To use this function across multiple columns use apply"))
}

ci <- function(v,formean=TRUE){
    v<-v[!is.na(v)]
    v.l<-length(v)
    if(v.l<2){
        return(c(NA,NA,NA))
    }
    v.m<-mean(v)
    v.sd<-ifelse(formean,sd(v)/sqrt(v.l),sd(v))
    return(c(v.m-2*v.sd,v.m,v.m+2*v.sd))
}
