test.stzdesign <- function(d, sid.name, tid.name, verbose = FALSE){
  ##############################################################################
  ## takes "d" which is like "dlinks" in vam.multivariate.ancova and determines
  ## if all effects are identified.
  ## if verbose==TRUE, prints additional diagnostics
  ##
  ## value: a list
  ## ok:     TRUE if all effects identified, FALSE otherwise
  ## nbad:   number of undefined effects         (0 if ok == TRUE)
  ## badsid: vector of sid.names in singular set (NULL if ok==TRUE)
  ## badtid: vector of tid.names in singular set (NULL if ok==TRUE)
  ##############################################################################

  ## check and homogenize names
  stopifnot( all( c(sid.name,tid.name,"iblock","frac","reff") %in% names(d) ) )
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
  stopifnot( min(sapply(d, pobs)) == 1 )
  stopifnot( min(d$frac >= 0) && max(d$frac <= 1) )
  stopifnot( lenu(d$tid) == nrow(unique(d[,c("tid","reff")])) )

  ## check GF and branch on whether it is full rank or not
  gf <- create.stzdesign(d, sid.name = "sid", tid.name = "tid")$gf
  b <-  coef(lm(rep(1,nrow(gf)) ~ gf - 1))
  whichbmis <- sort(which(is.na(b)))

  if(all(!is.na(b))){
    res <- list(ok=TRUE, nbad = 0, droppedtid = NULL, badsid = NULL, badtid = NULL, tidclusters = NULL)
  } else {
    n <- substr(names(b)[whichbmis], 3, 100000)
    nbad <- length(n)

    badcols <- tidclusters <- vector(nbad, mode="list")
    names(badcols) <- names(tidclusters) <- n
    for(i in 1:nbad){
      w <- whichbmis[i]
      cols.touse <- sort(setdiff(1:w, whichbmis[1:i]))
      m <- lm(gf[,w] ~ gf[,cols.touse] - 1)
      stopifnot(summary(m)$r.sq == 1.0)
      cf <- coef(m)
      names(cf) <- colnames(gf[,cols.touse])
      badcols[[i]]     <- as.matrix(cf[which(abs(cf) > 0.00000001)])
      tidclusters[[i]] <- c(n[i], rownames(badcols[[i]]))
    }

    if(verbose){
      print(badcols)
    }

    ## determine badtid and badsid
    badtid <- su(unlist(tidclusters))
    tmp <- subset(d, tid %in% badtid)
    badsid <- sort(unique(tmp$sid))
    res <- list(ok=FALSE, nbad = nbad, droppedtid = n, badsid = badsid, badtid = badtid, tidclusters = tidclusters)
    
    ## verify that after deleting all records for these students and teachers, gf is full rank
    tmp <- subset(d, !(tid %in% badtid) & !(sid %in% badsid))
    if(lenu(tmp$tid) > 1){
      gf <- create.stzdesign(tmp, sid.name = "sid", tid.name = "tid")$gf    
      stopifnot( ncol(gf) == qr(gf)$rank )
    }
  }
  return(res)
}
