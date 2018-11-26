## FLAG NEED TO CHANGE TO CARROLL NOTATION FOR REST
npmle_deconv_biv <- function(X, csem, umin = c(-5,-5), umax = c(5,5), ngrid = c(100,100), lambda = 0.0005, lltol = 1e-7, psmall = 0.000005, discrete=FALSE, quietly = FALSE, startmiddle=FALSE){
  ## bivariate nonparametric MLE of distribution of error-free variables "(U1,U2)" given
  ## error-prone variables "(X1,X2)" where Xi = Ui + ei and ei|Ui ~ N(Ui, csem_i^2(Ui))
  ##
  ## see npmle_deconv for more details
  ##
  ## NOTE: csem is passed a list of functions csem[[i]] is csem fxn for variable i
  ## NOTE: reduced default on psmall because now we are dealing with joint probabilities
  ## NOTE: reduced the .05 on new value to 0.01 because it was taking very long
  ## NOTE: startmiddle==TRUE mean algorithm will start at values near middle of grid
  stopifnot( all(!is.na(X)) && is.function(csem[[1]]) && is.function(csem[[2]]) && all(umin < umax) && (lambda > 0) && (lltol > 0) )

  ## NOTE: we collapse over repeat X values when possible, using "counts" to contribute to likelihood
  .X <- unique(X)
  .X <- .X[order(.X[,1], .X[,2]),]
  nX <- nrow(.X)
  .counts <- sapply(1:nX, function(i){ sum( (X[,1] == .X[i,1]) & (X[,2] == .X[i,2]) ) })
  stopifnot(sum(.counts) == nrow(X))
  X <- .X

  ## functions to map probabilities to reduced-dimension unconstrained scale, and inverse
  theta_to_p <- function(theta){
    if(length(theta) == 0){
      return(1)
    } else {
      p <- exp(theta) / (1 + sum(exp(theta)))
      return(c(p, 1-sum(p)))
    }
  }
  
  p_to_theta <- function(p){
    stopifnot( all(p > 0) && (abs((1 - sum(p))) < 1e-9) )
    if(length(p) == 1){
      return(numeric(0))
    } else {
      last <- length(p)
      return(log(p[1:(last-1)] * (1 + ((1-p[last])/p[last])) ))
    }
  }
  
  ## CHECKS:
  ## p <- c(0.1, 0.2, 0.4, 0.3); print(ma(p - theta_to_p(p_to_theta(p))))
  ## theta <- c(-3.2, 1.0, -1.0); print(ma(theta - p_to_theta(theta_to_p(theta))))
  ## print(ma(1 - theta_to_p(p_to_theta(1))))
  ## p_to_theta(theta_to_p(numeric(0)))

  ## build (nX x ngrid) matrix of conditional densities by multipling due to local independence, and
  ## then vec'ing the grid.  NOTE there is probably a better way to do this
  grid <- expand.grid(seq(from=umin[1], to=umax[1], length=ngrid[1]), seq(from=umin[2], to=umax[2], length=ngrid[2]))
  names(grid) <- c("u1","u2")
  ngrid <- nrow(grid)
  grid$index <- 1:ngrid
  grid$csem1 <- csem[[1]](grid$u1)
  grid$csem2 <- csem[[2]](grid$u2)  

  fxu <- matrix(0, nrow=nX, ncol=ngrid)
  if(!discrete){
    for(i in 1:nX){
      fxu[i,] <- dnorm(rep(X[i,1],ngrid), mean=grid$u1, sd = grid$csem1) * dnorm(rep(X[i,2],ngrid), mean=grid$u2, sd = grid$csem2)
    }
  } else { ## uses pxgu() and conditional independence of (X1,X2) given (U1,U2)
    ## FLAG FIXME this is inefficient
    .tab1 <- su(.X[,1])
    .tab2 <- su(.X[,2])
    p1 <- lapply(1:ngrid, function(i){ pxgu(grid$u1[i], grid$csem1[i], .tab1) })
    p2 <- lapply(1:ngrid, function(i){ pxgu(grid$u2[i], grid$csem2[i], .tab2) })
    p1 <- do.call("cbind", lapply(p1, function(x){ x[,"p"] })) ## sorted unique scale scores on rows, conditional probs at each grid point on columns
    p2 <- do.call("cbind", lapply(p2, function(x){ x[,"p"] })) ## sorted unique scale scores on rows, conditional probs at each grid point on columns
    for(i in 1:nX){
      fxu[i,] <- p1[which(.tab1 == X[i,1]),] * p2[which(.tab2 == X[i,2]),]
    }
    rm(p1,p2)
    gc()
  }
  
  ## negative log likelihood for a set of probabilities given ".inds" which are
  ## indices of "grid" and which are continually updated
  ## NOTE: updated to use counts
  negll <- function(theta){
    -sum(.counts * log(as.vector(fxu[,.inds] %*% theta_to_p(theta))))
  }

  ## find initial best grid point.  NOTE this will often want the u with the
  ## biggest CSEM because a single point is inconsistent with most X unless the
  ## CSEM is large. the exception is if we pick ugrid to be very large, then it
  ## will find interior point.  however even if it picks an extreme starting
  ## point, that point will be dropped later if it is too far in the tails.
  ll <- apply(fxu, 2, function(x){ sum(.counts * log(x)) })
  if(startmiddle){
    .inds <- which.min( (grid$u1 - mean(grid$u1))^2 + (grid$u2 - mean(grid$u2))^2 )
  } else {
    .inds <- which.max(ll)
  }
  .probs   <- 1
  .eligible <- setdiff(1:ngrid, .inds)
  ll.current <- ll[.inds]
  
  .history <- list()
  tmp <- grid[.inds,]
  tmp$p <- 1
  tmp$ll <- ll.current
  tmp$eu1 <- tmp$u1
  tmp$eu2 <- tmp$u2
  tmp$varu1 <- tmp$varu2 <- 0
  tmp$cor   <- 0
  .history[[1]] <- tmp

  ## now add grid points per the algorithm until there is no improvement
  done <- FALSE  
  while(!done){
    ## evaluate each eligible point with a weight "lambda"
    part.static   <- as.vector(fxu[,.inds,drop=FALSE] %*% ((1 - lambda)*.probs))
    part.total   <- matrix(part.static, ncol=length(.eligible), nrow=nX) + (fxu[,.eligible] * lambda)
    ll.candidate <- apply(part.total, 2, function(x){ sum(.counts * log(x)) })
    if(all(ll.candidate - ll.current < lltol)){
      done <- TRUE
    } else {
      w <- which.max(ll.candidate - ll.current)
      .inds <- c(.inds, .eligible[w])
      .eligible <- setdiff(.eligible, .eligible[w])
      
      ## set starting value: mass 0.01 at new value, normalized p on existing values
      o <- optim(p_to_theta(c(0.99*(.probs/sum(.probs)), 0.01)), negll, method  = "BFGS", control = list(maxit=1000, trace=5*(1 - as.numeric(quietly)), REPORT = 1))
      .probs <- theta_to_p(o$par)
      ll.current <- -o$value
      
      ## sometimes we might have picked up an early grid point but now it has virtually no probability, so dump it
      w <- which(.probs < psmall)
      if(length(w) > 0){
        .probs <- .probs[-w]
        .probs <- .probs/sum(.probs)
        .inds  <- .inds[-w]
      }
      
      ## reorder the grid if need be
      o <- order(.inds)
      .inds <- .inds[o]
      .probs <- .probs[o]

      ## summarize the state
      tmp <- grid[.inds,]
      tmp$p <- .probs
      tmp$ll <- ll.current
      tmp$eu1 <- sum(tmp$u1 * tmp$p)
      tmp$eu2 <- sum(tmp$u2 * tmp$p)
      tmp$varu1 <- sum(tmp$u1^2 * tmp$p) - (su(tmp$eu1))^2
      tmp$varu2 <- sum(tmp$u2^2 * tmp$p) - (su(tmp$eu2))^2
      tmp$cor   <- (sum(tmp$u1 * tmp$u2 * tmp$p) - (su(tmp$eu1) * su(tmp$eu2))) / sqrt( su(tmp$varu1) * su(tmp$varu2))
      .history[[length(.history)+1]] <- tmp
      
      if(!quietly){
        print(.history[[length(.history)]])
        cat("\n\n\n")      
      }
    }
  }

  ## calculate nonparametric estimate of conditional mean fxn
  ## WARNING- extremely noisy
  pu <- .history[[length(.history)]]
  tmp <- data.frame(u1 = su(pu$u1), eu2 = 0)
  for(i in 1:nrow(tmp)){
    z <- subset(pu, u1 == tmp$u1[i])
    tmp$eu2[i] <- weighted.mean(z$u2, w = z$p)
  }
  ## a different way using moving average to smooth
  ## tmp <- data.frame(u1 = seq(from=min(pu$u1), to=max(pu$u1), length=1000), eu2 = 0)
  ## for(i in 1:nrow(tmp)){
  ##  z <- subset(pu, (u1 > tmp$u1[i] - 0.6) & (u1 < tmp$u1[i] + 0.6))
  ##  tmp$eu2[i] <- weighted.mean(z$u2, w = z$p)
  ## }
  
  ## return
  return(list(.history = .history, pu = pu, eu2gu1 = tmp))
}
