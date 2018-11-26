sem.approx <- function(x, deg=5, plotfile="/tmp/sem.pdf", sem=TRUE){
  ## x is matrix of unique (score, SEM) pairs.  this finds a polynomial of degree
  ## "deg" that minimizes discrepancies after accounting for discreteness of
  ## reported function.  if sem==TRUE, finds approximating function to SEM.
  ## if sem==FALSE, finds approximating function to precision = 1/(SEM^2)
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
  ## create groups of identical SEMs corresponding to consecutive scale scores
  x$semgroup <- 1
  for(j in 2:nrow(x)){
    if(abs(x$sem[j] - x$sem[j-1]) < 0.00001){
      x$semgroup[j] <- x$semgroup[j-1]
    } else {
      x$semgroup[j] <- x$semgroup[j-1] + 1
    }
  }
  x$w <- ave(rep(1, nrow(x)), x$semgroup, FUN = sum)
  ## create function that calculates mean, ACROSS GROUPS, of squared
  ## discrepancies between reported objective value (sem or prec) and the
  ## average of the approximating function within the group.  we weight by the
  ## number of observations in each SEM group.  we use relative discrepancies
  ## where we upweight low SEM/high precision cases so we need separate function
  pows <- matrix(0:deg, ncol=(deg + 1), nrow=nrow(x), byrow=TRUE)
  u    <- matrix(x$u,   ncol=(deg + 1), nrow=nrow(x), byrow=FALSE)
  X    <- u^pows
  .x   <- x

  if(sem){
    f <- function(beta){
      .f <- X %*% beta
      .x$diff <- (.f - x$objective) / x$sem
      .x$meandiff <- ave(.x$diff, .x$semgroup, FUN = mean)
      tmp <- unique(.x[,c("semgroup","sem","w","meandiff")])
      weighted.mean(tmp$meandiff^2, w = tmp$w)
    }
  } else {
    f <- function(beta){
      .f <- X %*% beta
      .x$diff <- sqrt(x$prec) * (.f - x$objective)
      .x$meandiff <- ave(.x$diff, .x$semgroup, FUN = mean)
      tmp <- unique(.x[,c("semgroup","sem","w","meandiff")])
      weighted.mean(tmp$meandiff^2, w = tmp$w)
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
