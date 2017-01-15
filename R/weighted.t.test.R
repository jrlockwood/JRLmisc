weighted.t.test <- function(x,w,group){
    ## two-sample t-test with weighted data - uses normal approximation
    ## x is full vector of values, w is vector of weights, group is group indicator (must be 0/1)
    ## result is group1 - group0  
    if( !all( sort(unique(group)) == c(0,1) )){
        stop("invalid group")
    }
    
    x0 <- x[which(group==0)]
    w0 <- w[which(group==0)]
    n0 <- length(x0)
    x1 <- x[which(group==1)]
    w1 <- w[which(group==1)]
    n1 <- length(x1)
    
    xbar0 <- weighted.mean(x0,w0)
    xbar1 <- weighted.mean(x1,w1)
    xbar.diff <- xbar1 - xbar0
    squared.deviations <- c( (x0 - xbar0)^2, (x1 - xbar1)^2 )
    sigsq.hat <- weighted.mean(squared.deviations,c(w0,w1))
    v.diff <- (sum(w0^2)/(sum(w0))^2 + sum(w1^2)/(sum(w1))^2)
    deff <- v.diff/( 1/n0 + 1/n1 )
    se.diff <- sqrt(v.diff * sigsq.hat)
    t.stat <- xbar.diff/se.diff
    p.val <- 2*(1 - pnorm(abs(t.stat)))
    return(c(n0 = n0, n1 = n1, xbar0 = xbar0, xbar1 = xbar1, xbar.diff = xbar.diff, sigsq.hat = sigsq.hat, deff = deff, se.diff = se.diff, t.stat = t.stat, p.val = p.val))
}
