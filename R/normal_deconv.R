## FLAG REMEMBER THIS ALLOWS FOR A CSEM FUNCTION SO IS NOT COMPLETELY
## OBSOLETE AT THIS STAGE, unlike the homoskedastic case which is available
## in closed form
normal_deconv <- function(X, csem, umin = -5, umax = +5, ngrid=2000, discrete = FALSE, quietly=FALSE){
  ## estimation of parameters of normal distribution of error free variables "U"
  ## given error-prone variable "X" where X = U + e and e|U ~ N(U, csem^2(U))

  ## discrete: if TRUE, treat X as discrete and use pxgu() to evaluate probabilities -
  ## see below for how this is made efficient.  we tried various things and settled on approxfun
  ##
  ## NOTE: sometimes had trouble with it wanting tausq too big so we start with a smaller value
  stopifnot( all(!is.na(X)) && is.function(csem) && (umin < umax))

  ## NOTE: we collapse over repeat X values when possible, using "counts" to contribute to likelihood
  .mean <- mean(X)
  .var  <- var(X)
  stopifnot(.var < 100000)
  .X <- su(X)
  nX <- length(.X)
  .counts <- sapply(.X, function(x){ sum(X == x) })
  stopifnot(sum(.counts) == length(X))
  X <- .X

  ## NOTE: if discrete, it is extremely inefficient to evaluate the likelihood
  ## by calling pxgu() over and over.  so we build pre-made evaluations on a
  ## grid and then pull pieces out of it by using whatever are the nearest
  ## values of U that integrate() wants to use.
  ## NOTE: even this turned out to be too slow, so we use approxfun instead to build functions
  if(discrete){
    grid <- seq(from = umin, to = umax, length=ngrid)
    tmp <- lapply(grid, function(u){ pxgu(u, csem(u), .X)})
    stopifnot(all(sapply(tmp, function(x){ all(x[,1] == .X) })))
    fxu <- do.call("cbind", lapply(tmp, function(x){ x[,2] }))

    fxufun <- vector(nX, mode="list")
    for(i in 1:nX){
      fxufun[[i]] <- approxfun(grid, fxu[i,])
    }
    rm(fxu)
    gc()
  }
  
  ## define negative log likelihood depending on discrete==TRUE/FALSE  
  if(!discrete){
    negll <- function(p){
      mu <- p[1]
      tausq <- exp(p[2])
      if( (tausq < 0.000001) || (tausq > 100000) ){
        return(1e40)
      }
      like.i <- sapply(X, function(x){
        integrate(f=function(u){dnorm(x, mean = u, sd = csem(u)) * dnorm(u, mean = mu, sd = sqrt(tausq))},lower = umin, upper=umax, stop.on.error = FALSE)$value
      })
      -sum(.counts * log(like.i))
    }
  } else {
    negll <- function(p){
      mu <- p[1]
      tausq <- exp(p[2])
      if( (tausq < 0.000001) || (tausq > 100000) ){
        return(1e40)
      }
      like.i <- sapply(X, function(x){
        integrate(f=function(u){
          ## pxgu accepts only scalars but integrate() requires f evaluate on a vector so we need to build a loop in here
          ## NOTE THIS IS VERY INEFFICIENT!
          ## res <- rep(0, length(u))
          ## for(i in 1:length(u)){
          ##  tmp <- pxgu(u[i], csem(u[i]), X)
          ##  res[i] <- tmp[which(tmp[,1]==x),"p"]
          ## }
          ##
          ## BETTER WAYS
          ## res <- fxu[which(X==x), sapply(roundto(u,grid), function(z){ which(grid==z) })]
          ## res <- fxu[which(X==x), sapply(u, function(z){ which.min(abs(grid-z)) })]
          ## BEST WAY
          fxufun[[which(X==x)]](u) * dnorm(u, mean = mu, sd = sqrt(tausq))}, lower = umin, upper=umax, stop.on.error = FALSE)$value
      })
      -sum(.counts * log(like.i))
    }
  }

  .trace <- 5*(1 - as.numeric(quietly))
  r <- optim(par=c(.mean, log(0.8* .var)), fn = negll, method = "BFGS", control = list(maxit=500, trace=.trace, REPORT=1))
  stopifnot(r$convergence == 0)
  return(list(mu = r$par[1], tausq = exp(r$par[2]), ll = -r$value))
}
