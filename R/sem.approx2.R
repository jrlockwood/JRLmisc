sem.approx2 <- function(x, deg=5, plotfile="/tmp/sem.pdf", sem=TRUE){
  ## simpler version of sem.approx that uses a different objective function
  x <- as.data.frame(x)
  names(x) <- c("u","sem")
  x$prec <- 1.0 / (x$sem^2)
  if(sem){
    x$objective <- x$sem
  } else {
    x$objective <- x$prec
  }
  stopifnot(nrow(x) == lenu(x$u))
  x <- x[order(x$u),]

  pows <- matrix(0:deg, ncol=(deg + 1), nrow=nrow(x), byrow=TRUE)
  u    <- matrix(x$u,   ncol=(deg + 1), nrow=nrow(x), byrow=FALSE)
  X    <- u^pows
  .x   <- x

  if(sem){
    f <- function(beta){
      sum( ((X %*% beta - x$objective) / x$sem)^2 )
    }
  } else {
    f <- function(beta){
      sum( ((X %*% beta - x$objective) * sqrt(x$prec))^2 )
    }
  }
    
  ## start optimizer at lm fit
  par0 <- coef(lm(chartoform("objective", paste("I(u^",1:deg,")",sep="")), data=x))
  r <- optim(par0, f, method = "BFGS", control = list(maxit=500, trace=5, REPORT=1))
  stopifnot(r$convergence == 0)
  x$fitted <- X %*% r$par
  
  ## make plots of fits
  .n <- 1000
  useq <- seq(from = min(x$u) - 1.0*sd(x$u), to = max(x$u) + 1.0*sd(x$u), length=.n)
  u    <- matrix(useq,   ncol=(deg+1),   nrow=.n, byrow=FALSE)
  pows <- matrix(0:deg, ncol=(deg + 1), nrow=.n, byrow=TRUE)
  X    <- u^pows
  if(sem){
    fits.sem <- X %*% r$par
    fits.prec <- 1.0 / (fits.sem^2)
  } else {
    fits.prec <-  X %*% r$par
    fits.sem  <- 1.0 / sqrt(fits.prec)
  }
  pdf(plotfile)
  plot(x$u, x$sem,  xlim = range(useq), ylim = range(na.omit(c(fits.sem, x$sem))), ylab="SEM"); lines(useq, fits.sem, col="blue")
  plot(x$u, x$prec, xlim = range(useq), ylim = range(na.omit(c(fits.prec, x$prec))), ylab="PRECISION"); lines(useq, fits.prec, col="blue")
  dev.off()
  
  return(list(p = r$par, x = x))
}
