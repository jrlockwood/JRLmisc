partial_integrate <- function(fx, x){
    ## returns vector of partial integrals of fx along x grid using trapezoid rule.
    ## uses convention that first partial integral = 0
    stopifnot(length(fx)==length(x))
    cumsum(c(0, diff(x) * (fx[-length(fx)] + fx[-1]) / 2))
}
## test:
## x <- seq(from=-1, to=2, length=10000)
## 
## fx <- x
## plot(x, partial_integrate(fx, x))
## lines(x, x^2/2 - x[1]^2/2, col="red")
## 
## fx <- x^2
## plot(x, partial_integrate(fx, x))
## lines(x, x^3/3 - x[1]^3/3, col="red")
## 
## fx <- dnorm(x, sd=2)
## plot(x, partial_integrate(fx, x))
## lines(x, pnorm(x,sd=2) - pnorm(-1,sd=2), col="red")
