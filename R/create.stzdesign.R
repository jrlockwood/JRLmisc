create.stzdesign <- function(d, sid.name, tid.name, normalize.sfrac=TRUE, sparse=FALSE, use.sparse.model.matrix=FALSE){
  ######################################################################
  ## function to create a sum to zero constrained effects design matrix,
  ## with fractions.
  ##
  ## INPUTS:
  ## d = data.frame with student id (sid.name), teacher id (tid.name),
  ## which must be characters,
  ## iblock, reff, and fraction.  "iblock" is instructional block
  ## indicating at which level we are accounting student instruction and
  ## could be things like subject, year, or subject*year.  "reff" is
  ## teacher reference collection and will often be the same as iblock but
  ## this is not required
  ##
  ## normalize.sfrac = logical indicating how to proceed if sum of
  ## fractions within (sid.name, iblock) exceeds one.  this occasionally
  ## happens due to overlap of beginning and ends of spells.  if TRUE will
  ## normalize these cases to sum to 1; if FALSE it will abort.
  ##
  ## if sparse, uses sparse matrices.  further, if use.sparse.model.matrix,
  ## uses more efficient process to build design matrix; latter only
  ## available for case of no fractional linkage and one observation
  ## per student per iblock
  ##
  ## OUTPUT: a list with components
  ##
  ## gf        = design matrix with one record per (sid.name, iblock) pair
  ##             occuring in the data, sorted by those factors, parameterized with
  ##             sum-to-zero constrained effects and accounting for frac linkages
  ## usid      = unique (sid.name,iblocks)
  ## tmap      = how teachers are mapped to numeric ids (should not be needed)
  ## maxsfrac  = maximum observed sfrac within (sid.name, iblock)
  ## rcn       = reciprocal condition number of gf'gf
  ######################################################################
  if(use.sparse.model.matrix && !sparse){
    stop("use.sparse.model.matrix requires sparse==TRUE")
  }
  
  ## check and homogenize names
  stopifnot(all(c(sid.name,tid.name,"iblock","frac","reff") %in% names(d)))
  stopifnot(all(sapply(d[,c(sid.name, tid.name)], is.character)))
  if( tid.name != "tid" ){
    d <- d[,setdiff(names(d), "tid")] ## kill any existing "tid"
    names(d)[names(d) == tid.name] <- "tid"
  }
  if( sid.name != "sid" ){
    d <- d[,setdiff(names(d), "sid")] ## kill any existing "sid"
    names(d)[names(d) == sid.name] <- "sid"
  }

  ## checks on validity of d
  if(min(sapply(d, pobs)) < 1){
    stop("missing values not allowed")
  }
  if( min(d$frac < 0) || max(d$frac > 1) ){
    stop("invalid fractions")
  }
  if( lenu(d$tid) != nrow(unique(d[,c("tid","reff")]))){
    stop("invalid reff")
  }

  ## checks on validity of d if use.sparse.model.matrix==TRUE
  ## need to have one observation per student per iblock and the fractions
  ## must all be 1.
  if(use.sparse.model.matrix){
    if(nrow(unique(d[,c("sid","iblock")])) != nrow(d)){
      stop("use.sparse.model.matrix requires one observation per student per iblock")
    }
    if( (min(d$frac) < 0.999999) || (max(d$frac) > 1.000001) ){
      stop("use.sparse.model.matrix requires all fractions == 1")
    }
  }
  
  ## create "sfrac" which is sum of fractions within (sid.name, iblock)
  ## and make sure sfrac is valid
  d$sidi  <- paste(d$sid, d$iblock, sep="-")
  nsidi   <- lenu(d$sidi)
  d$sfrac <- ave(d$frac, d$sidi, FUN=sum)
  maxsfrac <- max(d$sfrac)
  if( (maxsfrac > 1) & !normalize.sfrac){
    stop("maxsfrac exceeds one and normalization declined")
  }
  if( (maxsfrac > 1) &  normalize.sfrac){
    wh <- d$sfrac > 1
    d$frac[wh]  <- d$frac[wh] / d$sfrac[wh]
    d$sfrac[wh] <- 1.0
  }

  ## replace reff with 1:K mapping and move "reff" to "reff.orig"
  tmp <- unique(d[,"reff",drop=FALSE])
  tmp <- tmp[order(tmp$reff),,drop=FALSE]
  tmp$reffi <- 1:nrow(tmp)
  d <- merge(d, tmp)
  names(d)[names(d)=="reff"]  <- "reff.orig"
  names(d)[names(d)=="reffi"] <- "reff"
  
  ## create mappings of teachers to positional ID "tidn" within reff.
  ## tidn == 0 is assigned to first teacher in each reff and they become holdout
  tmp <- unique(d[,c("tid","reff")])
  ntch <- nrow(tmp)
  tmp <- split(tmp, tmp$reff)
  nreff <- l(tmp)
  reffcounts <- sapply(tmp, nrow)
  if(ntch != sum(reffcounts)){ stop() }
  
  tmp <- lapply(tmp, function(x){
    x <- x[order(x$tid),]
    x$tidn <- as.integer(1:nrow(x) - 1)
    x$reffcount <- nrow(x)
    x
  })
  tmap <- do.call("rbind",tmp)
  rownames(tmap) <- 1:nrow(tmap)

  d <- merge(d, tmap[,c("tid","tidn","reffcount")], by="tid", all.x=TRUE)
  d <- d[order(d$sid, d$iblock, d$tid), c("sid","sidi","iblock","tid","tidn","frac","sfrac","reff","reffcount")]
  d$index <- as.integer(1:nrow(d))

  ## create "se" which has starting and ending index for each reff
  if(nreff == 1){
    se <- matrix(c(1, ntch), nrow=1)
  } else {
    .c <- cumsum(reffcounts)
    se <- cbind( 1 + c(0, .c[-nreff]), .c)
  }
  rownames(se) <- 1:nrow(se)
  colnames(se) <- c("start","end")
  
  ## assign [G,F] parts to each (sid.name, iblock)
  ##
  ## NOTE: this reflects decisions that dan and i made on 8/19/2009 as
  ## follows: let "s" = sum(frac) within (sid.name, iblock).  s==1 means we
  ## know how to allocate 100% of this student's instruction and s<1 means
  ## that some part of this student's instruction is missing data and we
  ## need to decide how to allocate.
  ##
  ## case 1: s==1.  we allocate to both G and F according to frac. sum(G)
  ## across reffs is thus 1.
  ##
  ## case 2: s < 1 and single reff: G gets 1 and F gets the allocation
  ## according to frac.  this is consistent with the assumption that the
  ## (1-s) unknown portion of the student's instruction came from the
  ## average teacher in this reff
  ##
  ## case 3: s < 1 and multiple reff: G for the linked reffs get their
  ## fraction divided by s, and F gets the allocation according to
  ## frac.  this is consistent with the assumption that the (1-s) of the
  ## student's instruction came from average teachers in the linked reffs
  ## in proportion to the amount of instruction we know they got in each
  ## reff
  ##
  ## the upshot of this is that teachers always get frac and the reff
  ## indicator always gets frac/s so that sum across reffs for each
  ## kid is 1 regardless of s
  gfrownames <- su(d$sidi)
  gfcolnames <- as.character(tmap$tid)
  gfcolnames[tmap$tidn==0] <- paste("g",1:nreff,sep="")
  
  if(!use.sparse.model.matrix){
    d <- split(d[,c("tidn","reff","reffcount","frac","sfrac")], d$sidi)
    ## gf will be list containing either rows of final matrix in dense form (if
    ## !sparse), or the pieces necessary to construct a sparse representation of
    ## it (if sparse).
    gf <- lapply(d, function(x){
      thisgf <- rep(0.0, ntch)
      s <- unique(x$sfrac)
      for(i in 1:nrow(x)){
        thisgf[ se[x$reff[i],1] ] <- thisgf[ se[x$reff[i],1] ] + (x$frac[i] / s)
        if(x$reffcount[i] > 1){ ## non-trivial reff; more to do
          if( (x$tidn[i] > 0) ){
            thisgf[ se[x$reff[i],1] + x$tidn[i] ] <- thisgf[ se[x$reff[i],1] + x$tidn[i] ] + x$frac[i]
          } else {
            thisgf[ (se[x$reff[i],1] + 1):(se[x$reff[i],2]) ] <- thisgf[ (se[x$reff[i],1] + 1):(se[x$reff[i],2]) ] - x$frac[i]
          }
        }
      }
      if(sparse){ ## reduce thisgf to key information
        j <- as.integer(which(abs(thisgf) > .Machine$double.eps))
        thisgf <- list(nzcols = j, nz = thisgf[j])
      }
      thisgf
    })

    if(!sparse){
      gf <- do.call("rbind", gf)
    } else {
      ## note we build the transpose of the matrix and then transpose to match
      ## Matrix's column-oriented sparse format, since we calculated gf by rows.
      tmp   <- Matrix(0, nrow=ntch, ncol=l(d))
      tmp@x <- as.vector(unlist(lapply(gf, function(x){x$nz})))
      tmp@p <- as.integer(c(0, cumsum(as.vector(sapply(gf, function(x){ l(x$nz) })))))
      ## check: ma(diff(tmp@p) - sapply(gf, function(x){ l(x$nz)}))
      tmp@i <- as.vector(unlist(lapply(gf, function(x){x$nzcols - as.integer(1)})))
      ## note in tmp@a - need to make them zero-based indices, hence subtraction of 1
      gf <- t(tmp)
      rm(tmp)
      tmp <- gc() ## to avoid printing
    }
  } else{
    ## use.sparse.model.matrix==TRUE; use entirely different logic to get gf
    f <- function(x, reffname){
      if(lenu(x$tid) > 1){
        x$tid <- as.character(x$tid)
        sortid <- su(x$tid)
        firstid <- sortid[1]
        x$tid[x$tid==firstid] <- paste("ZZZZ",firstid,sep="")
        x$tid <- factor(x$tid)
        mr <- sparse.model.matrix(~tid, data=x, contrasts.arg = list(tid = "contr.sum"))
        colnames(mr) <- c(paste("g",reffname,sep=""), sortid[-1])
      } else {
        mr <- Matrix(1.0, ncol=1, nrow=nrow(x))
        colnames(mr) <- c(paste("g",reffname,sep=""))
      }
      return(mr)
    }

    tmp <- split(d[,c("index","reff","tid")], d$reff)
    o   <- order(unlist(lapply(tmp, function(x){x$index})))
    tmp <- mapply(f, tmp, names(tmp), SIMPLIFY=FALSE)
    gf   <- bdiag(tmp)[o,]
    if(!all(as.vector(unlist(lapply(tmp, colnames))) == gfcolnames)){ stop("error in names when use.sparse.model.matrix") }
    rm(tmp)
    gc()
  }

  if( (nrow(gf) != nsidi) | (ncol(gf) != ntch) ){ stop("gf invalid") }
  rownames(gf) <- gfrownames
  colnames(gf) <- gfcolnames

  ## make sure reff columns all sum to one by row
  v <- range(apply(gf[,se[,1],drop=FALSE], 1, sum))
  stopifnot( (min(v) > 0.99999) && (max(v) < 1.00001) )

  ## reciprocal condition number
  rcn <- try( 1.0 / condest(crossprod(gf))$est, silent=TRUE)
  if(class(rcn) == "try-error"){
    rcn <- 0
    warning("GF IS NOT OF FULL RANK")
  }
  
  return(list(gf = gf, usid = gfrownames, tmap=tmap, maxsfrac = maxsfrac, rcn = rcn))
}
