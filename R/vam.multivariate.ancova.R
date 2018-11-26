vam.multivariate.ancova <- function(dlinks, dscores, yvar, xvars, xvars.min = NULL, thresh.minpatsize, thresh.mincnobs, thresh.K, sid.name, tid.name, plotfile, use.gls=FALSE, sparse=FALSE, cluster=NULL, irlstol=1e-8, rcntol = c(warn=1e-5, die=1e-10), ...){
  ## function to perform multivariate ancova VAM estimator for a cross-section
  ## of teachers using sum to zero teacher effects within reference collections
  ## and pattern-mixture design matrix for covariates
  ##
  ## INPUTS:
  ## dlinks = dataframe consisting of (at least) sid.name, tid.name, frac, reff
  ## dscores = dataframe consisting of (at least) sid.name, outcome "yvar", and prior scores/covariates "xvars"
  ## xvars.min: argument for create.pmixdesign.uni
  ## thresh.minpatsize = minimum allowable pattern count
  ## thresh.mincnobs   = minimum allowable number of observed xvars in a pattern
  ## thresh.K = cfreq must be >= (cnobs + 1)*thresh.K
  ## sid.name = field name in data corresponding to student id (must be character)
  ## tid.name = field name in data corresponding to teacher id (must be character)
  ## plotfile = path to diagnostic plot file
  ## use.gls: if TRUE, estimates model with GLS using separate residual variances by pattern
  ## sparse: if TRUE, uses sparse matrices.
  ## cluster: if non-NULL, indicates a single student-level variable at which to use random-effects.
  ##          the cluster variable must be in "dlinks" and must be unique at the student level
  ## irlstol: if use.gls && sparse, the tolerance used for IRLS algorithm
  ## rcntol: conditioning tolerance for [XGF]'[XGF] (minimum reciprocal condition number)
  ## ...: optional arguments passed to gls() or lme()
  ##
  ## OUTPUTS: a list with components
  ## est = estimates and standard errors
  ## coef = all model coeffs and standard errors
  ## modelstat = model statistics like R^2 etc
  ## upat = pattern table from pattern mixture function (if(!is.null(xvars)))
  ## callinfo = information from call
  ## preds = predicted scores with and without teacher effects.
  ##         note predictions "without" tx effx do still include reff means.
  ##
  ################################################################
  ## basic checks on inputs
  ################################################################
  stopifnot( sid.name %in% names(dlinks) )
  stopifnot( sid.name %in% names(dscores) )
  stopifnot( tid.name %in% names(dlinks) )
  stopifnot( is.character(dscores[,sid.name]) && all(sapply(dlinks[,c(sid.name, tid.name)], is.character)) )
  stopifnot( identical(su(dlinks[,sid.name]), su(dscores[,sid.name])) )
  stopifnot( nrow(dscores) == lenu(dscores[,sid.name]) )
  stopifnot( identical(1:nrow(dscores), order(dscores[,sid.name])) )
  stopifnot( all(c(yvar,xvars) %in% names(dscores)) )
  stopifnot( !any(is.na(dscores[,yvar])) )  
  stopifnot( !(is.null(xvars) && use.gls) )
  stopifnot( !any(strsplit(paste(c(xvars, unique(as.character(dlinks[,tid.name]))), collapse=""),"")[[1]] %in% c("-","+","%","/")) )
  use.cluster <- ifelse(is.null(cluster), FALSE, TRUE)
  if(use.cluster){
    stopifnot(!sparse) ## FLAG REMOVE THIS LATER
    stopifnot( cluster %in% names(dlinks) )
    dcluster <- unique(dlinks[,c(sid.name,cluster)])
    dcluster <- dcluster[order(dcluster[,sid.name]),]
    stopifnot( lenu(dlinks[,sid.name]) == nrow(dcluster) )
  }
  
  ################################################################
  ## create and check design matrices
  ################################################################
  ## teacher effect design
  dlinks$iblock <- 1
  rgf <- create.stzdesign(dlinks[,c(sid.name,tid.name,"iblock","frac","reff")], sid.name=sid.name, tid.name=tid.name, sparse=sparse)
  stopifnot( all(rgf$usid == paste(dscores[,sid.name], 1, sep="-")) )
  tmap <- rgf$tmap
  gfnames <- colnames(rgf$gf)
  stopifnot( nrow(rgf$gf) == nrow(dscores) )
  rcn.gf <- rgf$rcn
  
  ## prior score design.  NOTE: we will need to remove one column from
  ## cbind(X,GF) because they each contain something equilvant to an intercept:
  ## pattern means in X and reff means in GF.  we choose to remove the first
  ## pattern mean column.  this does not affect the teacher effect estimates or
  ## standard errors.
  ##
  ## 12/23/2010: addition of .px = indices of columns used to create "X only" predictions
  ## which include all X's, and teacher reff means but no individual teacher effects.
  ## if is.null(xvars), includes ONLY the reff means.
  if(!is.null(xvars)){
    rx <- create.pmixdesign.uni(d = dscores, xvars = xvars, minsize = thresh.minpatsize, xvars.min = xvars.min, sparse=sparse)
    rcn.x <- rx$rcn
    
    wh <- which(rx$info$cnobs < thresh.mincnobs)
    if(length(wh) > 0){
      cat(paste(wh,"\n"))
      stop("Pattern collapsing produced records not meeting thresh.mincnobs (indices above)")
    }
    
    tmp <- unique(rx$upat[,c("cpat","cpatlab","cnobs","cfreq")])
    tmp <- subset(tmp, cfreq < ( (cnobs + 1) * thresh.K))
    if(nrow(tmp) > 0){
      print(tmp)
      stop("Pattern collapsing produced collapsed patterns not meeting thresh.K (collapsed patterns above)")
    }
    
    upat <- rx$upat
    pat <- factor(rx$info$cpatlab)
    firstpat <- sort(unique(upat$cpatlab))[1]
    xnames <- setdiff(colnames(rx$X), firstpat)
    if(sparse){
      XGF <- cBind(rx$X[,xnames], rgf$gf)
    } else {
      XGF <- cbind(rx$X[,xnames], rgf$gf)
    }
    rm(rx)
  } else {
    rcn.x <- 1.0    
    upat <- NULL
    xnames <- character(0)
    XGF <- rgf$gf
  }
  rm(rgf)
  .n <- nrow(XGF)
  .p <- ncol(XGF)
  ## FLAG kludgy: define x-only positions by diffing against named teacher columns
  tmp <- subset(tmap, tidn > 0)
  .px <- which( !(colnames(XGF) %in% as.character(tmp$tid)) )
  
  #################################################################################
  ## fit model and create diagnostics.  fork on use.gls, sparse, and cluster
  ## NOTE: calculate rsq by hand because it is distorted in lm() by excluding intercept
  #################################################################################
  .Y <- dscores[,yvar]
  
  ## initial conditioning check
  cXGF <- crossprod(XGF)
  rcn.xgf <- 1.0 / condest(cXGF)$est
  stopifnot(rcn.xgf > rcntol["die"])
  if(rcn.xgf < rcntol["warn"]){
    warning("[XGF]'[XGF] MAY BE ILL-CONDITIONED")
  }
  
  ###################
  ## eight-way fork #
  ###################
  if(!sparse){
    m <- lm(.Y ~ XGF - 1, singular.ok=TRUE)
    stopifnot( !any(is.na(coef(m))) ) ## stop if singular

    if(!use.gls && !use.cluster){
      .coef  <- summary(m)$coef
      .fits  <- fitted(m)
      .fitsx <- as.vector(XGF[,.px,drop=FALSE] %*% coef(m)[.px])
      .resid <- resid(m)
      .sigma <- summary(m)$sigma
      .vbetahat <- (.sigma^2) * summary(m)$cov.unscaled
      rm(m)
    } else{
      tmp <- data.frame(.Y, pat, XGF)
      names(tmp) <- c(".Y","pat", paste("XGF",c(xnames,gfnames),sep=""))
      form <- as.formula(paste(".Y ~ ", paste(setdiff(names(tmp), c(".Y","pat")), collapse=" + "), "-1"))

      if(use.gls && !use.cluster){
        print(m  <- gls(model = form, data=tmp, weights=varIdent(form=~1|pat), ...))
      } else {
        tmp$cluster <- factor(dcluster[,cluster])
        if(!use.gls && use.cluster){
          print(m  <- lme(fixed = form, data=tmp, random = ~1|cluster, ...))
        } else {
          print(m  <- lme(fixed = form, data=tmp, random = ~1|cluster, weights=varIdent(form=~1|pat), ...))
        }
        m$varBeta <- m$varFix
      }
      
      .coef  <- summary(m)$tTable
      .fits  <- as.vector(fitted(m))
      .fitsx <- as.vector(XGF[,.px,drop=FALSE] %*% .coef[.px,1])
      .resid <- as.vector(resid(m))
      .sigma <- summary(m)$sigma
      .vbetahat <- m$varBeta
      rm(tmp, m)
    }
  } else if (!use.gls & sparse){  ## FLAG PICK UP HERE TO ADD CLUSTERING
    fac <- Cholesky(cXGF)
    bhat <- solve(fac, crossprod(XGF, .Y))

    .coef  <- as.matrix(data.frame(est = as.vector(bhat), se = -999.0))
    .fits <- as.vector(XGF %*% bhat)
    .fitsx <- as.vector(XGF[,.px] %*% bhat[.px,,drop=FALSE])
    .resid <- .Y - .fits
    .sigsq <- as.vector(crossprod(.resid)) / (.n - .p)
    .sigma <- sqrt(.sigsq)
    .vbetahat <- .sigsq * solve(fac)
    .coef[,2] <- sqrt(diag(.vbetahat))
    rownames(.vbetahat) <- colnames(.vbetahat) <- rownames(.coef) <- paste("XGF",c(xnames, gfnames),sep="")
    rm(fac)
  } else{
    fac <- Cholesky(cXGF)
    bcur <- solve(fac, crossprod(XGF, .Y))
    .resid <- .Y - as.vector(XGF %*% bcur)
    bhat <- bcur + 1

    cat("\n\nIRLS ITERATIONS (SPARSE GLS):\n")
    while(ma(bhat - bcur) > irlstol){
      bcur <- bhat
      Whalf <- Diagonal(.n, x = sqrt(ave(.resid^2, pat, FUN = function(x){ 1 / mean(x) })))
      WhalfXGF <- Whalf %*% XGF
      Whalfy <- Whalf %*% .Y
      fac <- Cholesky(crossprod(WhalfXGF))
      bhat <- solve(fac, crossprod(WhalfXGF, Whalfy))
      .resid <- .Y - as.vector(XGF %*% bhat)
      v <- tapply(.resid^2, pat, mean)
      cat(paste( paste(c(ma(bhat - bcur), as.vector(sqrt(v/v[1]))), collapse="   "), "\n") )
    }
    
    .coef  <- as.matrix(data.frame(est = as.vector(bhat), se = -999.0))
    .fits  <- as.vector(XGF %*% bhat)
    .fitsx <- as.vector(XGF[,.px] %*% bhat[.px,,drop=FALSE])    
    .resid <- .Y - .fits
    .sigma <- sqrt(as.vector(crossprod(.resid)) / (.n - .p))
    
    .r     <- Whalfy - as.vector(WhalfXGF %*% bhat)
    .vbetahat <- (as.vector(crossprod(.r)) / (.n - .p)) * solve(fac)
    .coef[,2] <- sqrt(diag(.vbetahat))
    rownames(.vbetahat) <- colnames(.vbetahat) <- rownames(.coef) <- paste("XGF",c(xnames, gfnames),sep="")
    rm(fac, Whalf, WhalfXGF, Whalfy)
  }
  rm(cXGF)
  tmp <- gc()

  .rsq   <- cor(.fits, .Y)^2
  .stats <- c(rsq = .rsq, arsq = 1 - ((1 - .rsq)*(.n - 1)/(.n - .p - 1)), sigma = .sigma)

  ## diagnostic plots
  postscript(plotfile)
  hist(.resid)
  qqnorm(.resid, main="normal QQ plot of residuals"); qqline(.resid)
  plot(.fits, .resid, xlab="fitted values", ylab="residuals"); abline(h=0, col="gray")
  plot(.fits, .Y, xlab="fitted values", ylab=yvar, main=round(.rsq,3)); abline(lm(.Y ~ .fits))

  ######################################################################
  ## get teacher estimates and se's
  ######################################################################
  tmapest        <- subset(tmap, tidn > 0)
  tmapest$cnames <- paste("XGF",tmapest$tid,sep="")

  est <- .coef
  est   <- as.data.frame(est[tmapest$cnames,1:2])
  names(est) <- c("est","se")
  est$tid <- substring(rownames(est),4)
  stopifnot( all(sort(est$tid) == sort(tmapest$tid)) )
  est    <- merge(tmap, est, by="tid", all.x=TRUE) ## includes reff hold-outs
  
  ## create estimates and se's for reff hold-out teachers.
  ## NOTE: if teacher is only teacher in their reff, estimate and SE ends up zero
  for(i in unique(est$reff)){
    est$est[ (est$reff == i) & (est$tidn==0) ] <- -sum(est$est[ (est$reff == i) & (est$tidn > 0) ])
    xn <- tmapest$cnames[tmapest$reff == i]
    one   <- matrix(1, nrow=1, ncol=length(xn))
    est$se[ (est$reff == i) & (est$tidn==0) ] <- as.vector(sqrt(one %*% .vbetahat[xn,xn] %*% t(one)))
  }
  stopifnot( ma(tapply(est$est, est$reff, sum)) < 0.00000001 )

  ## make some diagnostic plots for estimates (continuation of previous plotfile)
  hist(est$est, main="Estimates")
  caterpillar(x = est[,c("est","se")], main="Estimates with Error Bars")
  dev.off()

  ################################################################
  ## return results
  ################################################################
  return(list(est       = est,
              coef      = .coef,
              modelstat = .stats,
              upat      = upat,
              rcn       = c(x = rcn.x, gf = rcn.gf, xgf = rcn.xgf),
              callinfo  = list(yvar = yvar, use.gls = use.gls, sparse = sparse, cluster = cluster, irlstol = irlstol, rcntol = rcntol),
              preds     = data.frame(yhat = .fits, yhat.xonly = .fitsx)))
}
