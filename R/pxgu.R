pxgu <- function(u, csem, ss){
  ## evaluates conditional probability of each possible score in table "ss" under a model
  ## Z|U = u ~ N(u, csem), csem a scalar = conditional SEM at u
  ## X = f(Z) rounding Z to the nearest scale score in "ss".
  ##
  ## NOTE: had problems with the likelihood function induced by this, it because
  ## impossible to distinguish extreme U's because of the censoring.
  ## tried fixing by create probabilities by normalizing the densities at each "ss"
  ## but that led to even worse behavior in the NPMLE.  so we leave the original way
  stopifnot(is.numeric(u) && is.numeric(csem) && (csem > 0) && is.vector(ss) && all(!is.na(ss)) && (length(ss) >= 3))
  ss <- sort(ss)
  k <- length(ss)
  ## THIS IS PROBLEMATIC
  p <- pnorm(ss[-k] + diff(ss)/2, mean = u, sd = csem)
  tab <- cbind(grid=ss, p = c(p[1], diff(p), 1-p[(k-1)]))
  ##
  ## THIS FAILED
  ## p <- dnorm(ss, mean = u, sd = csem)
  ## tab <- cbind(grid=ss, p = p/sum(p))
  stopifnot(ma(sum(tab[,"p"]) - 1) < 1e-10)
  return(tab)
}

## testing discrete case - generate data from the discrete model and try to recover distribution of U
## NOTE: this is temperamental, dependent on grid size and tends to want to overestimate var(u)
## NOTE: as a result we ended up changing pxgu -- see comments there
## 
## n <- 25000
## u <- as.vector(scale(rnorm(n)))
## csem <- function(u){ .30 + 0.1*u^2 }
## ss   <- seq(from=-3, to = 3, length=50)
## X <- rep(0,n)
## for(i in 1:n){
##   tmp <- pxgu(u[i], csem(u[i]), ss)
##   X[i] <- sample(tmp[,"grid"], size=1, prob=tmp[,"p"])
## }
