create.pmixdesign.uni <- function(d, xvars, minsize, xvars.min = NULL, sparse=FALSE, drop.collinear.vars=TRUE, freqonly=FALSE){
  ## create design matrix appropriate for fitting regression with missing values
  ## using pattern mixture approach.
  ##
  ## d: dataframe or matrix
  ## xvars: names of the variables to be included in the patterns and
  ##   final design matrix.  NOTE: must be numeric
  ## minsize: smallest allowable pattern frequency; smaller ones are collapsed  
  ## xvars.min: names of variables that will be part of every pattern
  ##   including taking the place of the null pattern when collapsing fails.
  ##   must be subset of xvars and must be fully observed
  ## sparse: if TRUE, returns sparse version of matrix
  ## drop.collinear.vars: if TRUE, removes redundant columns of the design matrix
  ## freqonly: if TRUE, just returns table of pattern frequencies in data
  ##
  ## NOTE: the logic we take here, which might not be optimal in all settings,
  ## is that we collapse all small patterns back to existing big patterns.  we
  ## do not attempt to combine small patterns to make a new big pattern.  this
  ## logic is appropriate for test scores where there are dominant patterns
  ##
  ## result: list(X     = pattern expanded design matrix in order of d,
  ##              upat  = summary of all pattern info
  ##              info  = record-level information that might be useful
  ##              rcn   = reciprocal condition number of X'X
  d <- as.data.frame(d)
  
  stopifnot(all(xvars %in% names(d)))
  if(!is.null(xvars.min)){
    stopifnot(all(xvars.min %in% xvars))
    stopifnot(max(sapply(d[,xvars.min], nmis)) == 0)
  }

  infonames <- c("index","pat","nobs","cpat","cnobs","cpatlab")
  if(any(infonames %in% names(d))){
    stop("variables created in this function interfere with variables in d")
  }
  d$index <- 1:nrow(d)

  if(!all(sapply(d[,xvars], is.numeric))){
    stop("xvars must be numeric. convert factors to indicators")
  }
  d$pat   <- apply(as.matrix(d[,xvars]), 1, function(x){paste(as.numeric(!is.na(x)), collapse="")})
  
  ################################################
  ## gather unique patterns and their frequencies
  ################################################
  upat     <- as.data.frame(table(d$pat), stringsAsFactors=FALSE)
  names(upat) <- c("pat","freq")
  stopifnot( sum(upat$freq) == nrow(d) )
  upat$pat <- as.character(upat$pat)  
  upat <- upat[order(upat$freq),]
  rownames(upat) <- 1:nrow(upat)
  upat$nobs <- sapply(strsplit(upat$pat,""), function(x){sum(as.numeric(x))})
  if(freqonly){
    return(list(X = NULL, upat = upat, info = NULL, rcn = NULL))
  }
  if(max(upat$freq) < minsize){
    stop("no patterns with sufficient frequency; reduce minsize")
  }
  
  ################################################
  ## collapse patterns
  ## NOTE: changes made 12/23/2010 regarding xvars.min.  and note if there is a pattern
  ## containing only xvars.min that is below minsize, it gets called "collapsed" which is misnomer
  ################################################
  upat$collapsed <- 0
  upat$cpat <- upat$pat
  upat$cpat[upat$freq < minsize] <- ""
  
  big    <- subset(upat, freq >= minsize)
  small  <- subset(upat, freq <  minsize)
  nsmall <- nrow(small)
  
  if(nsmall > 0){
    small$collapsed <- 1
    for(i in 1:nrow(small)){
      p <- as.numeric(unlist(strsplit(small$pat[i],"")))
      ## to be eligible, the potential target cannot have a "1" anywhere the small
      ## has a "0".  put another way the number of "1"s in the potential target
      ## that overlap with the "1"s of the small must equal nobs for the target
      elig <- big
      elig$noverlap <- sapply(elig$pat, function(x){ sum( as.numeric(unlist(strsplit(x,""))) * p) })
      elig <- subset(elig, nobs==noverlap)
      ## if nothing is eligible, assign a pattern of unobserved or xvars.min.
      ## otherwise find the best pattern which maximizes nobs, and if ties, maximizes freq (for efficiency),
      ## and if further ties, arbitrarily takes the first
      if( (nrow(elig) == 0) && is.null(xvars.min) ){
        small$cpat[i] <- paste(rep(0, length(xvars)), collapse="")
      } else if( (nrow(elig) == 0) && !is.null(xvars.min) ){
        small$cpat[i] <- paste(as.integer(xvars %in% xvars.min), collapse="")
      } else {
        elig <- subset(elig, nobs == max(elig$nobs))
        small$cpat[i] <- subset(elig, freq == max(elig$freq))[1,"cpat"]
      }
    }
  }
  
  upat <- rbind(big, small)
  upat$cfreq <- ave(upat$freq, upat$cpat, FUN=sum)

  ## create cpat labels sorted by cfreq then by cpat and put onto upat
  tmp <- unique(upat[,c("cpat","cfreq"),drop=FALSE])
  if( (sum(tmp$cfreq) != nrow(d)) || (sum(upat$freq) != nrow(d)) ){ stop("freqs off")}
  tmp$cnobs <- sapply(tmp$cpat, function(x){sum( as.numeric(unlist(strsplit(x,""))))})
  tmp <- tmp[order(tmp$cfreq, tmp$cpat),]
  tmp$cpatlab <- rev(paste("p",lex(nrow(tmp)),sep=""))
  upat <- merge(upat, tmp[,c("cpat","cnobs","cpatlab")], by="cpat", all.x=TRUE)
  upat <- upat[order(-1*upat$cfreq, -1*upat$freq),c("pat","nobs","freq","collapsed","cpat","cpatlab","cnobs","cfreq")]
  rownames(upat) <- 1:nrow(upat)
  if( any(upat$cnobs > upat$nobs) ){stop("cnobs too big")}
  
  ## merge cpat and related info onto d
  d <- merge(d, upat[,c("pat","nobs","cpat","cnobs","cpatlab")], by="pat", all.x=TRUE)
  d <- d[order(d$cpatlab, d$index),]
  
  ###########################################
  ## create design matrix
  ###########################################
  u <- unique(upat[,c("cpat","cpatlab","cnobs")])
  X <- vector(nrow(u), mode="list")
  names(X) <- u$cpatlab
  d$.one <- 1.0
  
  for(i in 1:nrow(u)){
    xvarspat <- xvars[unlist(strsplit(u$cpat[i],""))=="1"] ## character(0) for all-missing pattern
    X[[i]] <- as.matrix(d[,c(".one",xvarspat),drop=F]) * as.numeric(d$cpatlab == u$cpatlab[i])  ## uses recycling rules
    X[[i]][is.na(X[[i]])] <- 0.0
    if(ncol(X[[i]])==1){ ## necessary because of how paste deals with character(0)
      colnames(X[[i]]) <- u$cpatlab[i]
    } else {
      colnames(X[[i]]) <- c(u$cpatlab[i], paste(u$cpatlab[i], xvarspat, sep="."))
    }
  }
  
  if(any(sapply(X, ncol) != (u$cnobs + 1))){stop("X failure")}
  Xnames <- lapply(X, colnames)
  X <- do.call("cbind",X)
  
  nobsx <- apply(X, 1, function(x){ sum(abs(x)!=0) })
  if( (max(nobsx) > (length(xvars) + 1)) || (min((d$nobs + 1) - nobsx) < 0) ){ stop("X failure") }
  if( ma(apply(X[,u$cpatlab, drop=FALSE], 1, sum) - 1) > 0.0000000000001 ) { stop("X failure") }

  ## put X back onto d.  this is inefficient but allows some checking.
  ## NOTE: if we remove this we will need to resort X correctly!!!
  if( any(names(X) %in% names(d)) ){ stop("X/d name interference") }
  keep <- colnames(X)
  d <- data.frame(d, X)
  d <- d[order(d$index),]
  rm(X)
  tmp <- gc()
  
  ## if there is more than one non-null pattern, then elements for the others
  ## should be appropriately zeroed
  allnames <- as.vector(unlist(Xnames))
  if(length(Xnames) > 1){
    for(i in 1:length(Xnames)){
      shouldbezero <- unique(c(as.matrix(subset(d, cpatlab==u$cpatlab[i], select=setdiff(allnames, Xnames[[i]])))))
      if( (length(shouldbezero) > 1) || (abs(shouldbezero) > 0.00000001) ){
        stop("X/d failure")
      }
    }
  }
  
  X <- as.matrix(d[,keep])

  ## collinearity check, calculate condition number and drop collinear vars if needed
  m <- lm.fit(X, rep(1,nrow(X)), singular.ok=TRUE)
  stopifnot(all(names(m$coef) == colnames(X)))
  v <- which(is.na(m$coef))
  if(length(v) == 0){        ## non-singular
    rcn <- 1.0 / condest(crossprod(X))$est
  } else {                   ## singular
    if(drop.collinear.vars){
      cat(paste("\ncreate.pmixdesign.uni: drop.collinear.vars dropping:",paste(names(v), collapse=" "), "\n\n"))
      X <- X[,-as.vector(v)]
      rcn <- 1.0 / condest(crossprod(X))$est
    } else {
      rcn <- 0.0
    }
  }
  
  if(sparse){
    X <- Matrix(X, sparse=TRUE)
  }
  rownames(X) <- NULL
  
  return(list(X = X, upat = upat, info = d[,infonames], rcn = rcn))
}
