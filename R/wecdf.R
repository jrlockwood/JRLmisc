wecdf <- function(x, w=NULL,lin.interp=TRUE){
    if(is.null(w))
        w = rep(1,length(x))/length(x)
    
    s.x = sort(x)
    s.w = w[order(x)]
    x.vals = unique(s.x)
    
    if(length(x.vals) != length(s.x)){
        x.mass = tapply(s.w, s.x, sum)
    }
    else
        x.mass = s.w
    x.cs = cumsum(x.mass)
    if(!lin.interp)
        ret = stepfun(x.vals, c(0,x.cs/tail(x.cs,1)))
    else{
        ret = list(fn =approxfun(x.vals, x.cs/tail(x.cs,1), yleft=0, yright=1))
        ret$fn.inverse <- approxfun(x.cs/tail(x.cs,1), x.vals, rule=2)
    }
    
    ret
}
