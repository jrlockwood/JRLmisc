get.strata <- function(d, tid.name, sid.name){
  ## d is dataset with teacher/student pairings.  function finds all
  ## strata and returns list where each component is a list giving the
  ## tid's and sid's of that stratum.  function calls get.stratum()
  ## defined within to get a single stratum.  tid.name is the name of
  ## the teacher id variable and sid.name is the name of the student
  ## id variable.

  get.stratum <- function(dsub){
    tid.set <- dsub$tid[1] ## take first teacher in dsub as seed
    sid.set <- unique(dsub$sid[ dsub$tid %in% tid.set ])
    l.old <- c(length(tid.set), length(sid.set))
    
    tid.set <- unique(dsub$tid[ dsub$sid %in% sid.set ])
    sid.set <- unique(dsub$sid[ dsub$tid %in% tid.set ])
    
    while( max( c(length(tid.set), length(sid.set)) - l.old ) > 0 ){
      l.old <- c(length(tid.set), length(sid.set))
      tid.set <- unique(dsub$tid[ dsub$sid %in% sid.set ])
      sid.set <- unique(dsub$sid[ dsub$tid %in% tid.set ])
    }
    return(list(tid = sort(tid.set), sid = sort(sid.set)))
  }

  ## clean up data
  stopifnot( all( c(tid.name, sid.name) %in% names(d) ) )
  if( tid.name != "tid" ){
    d <- d[,setdiff(names(d), "tid")] ## kill any existing "tid"
    names(d)[names(d) == tid.name] <- "tid"
  }
  if( sid.name != "sid" ){
    d <- d[,setdiff(names(d), "sid")] ## kill any existing "sid"
    names(d)[names(d) == sid.name] <- "sid"
  }
  d <- na.omit(d[,c("tid","sid")])

  ## calculate strata
  s <- vector(0, mode="list")
  wh <- 1
  s[[wh]] <- get.stratum(d)
  dsub <- subset(d, !(tid %in% s[[wh]]$tid)) ## equivalent to pulling out non-matching students

  while(nrow(dsub) > 0){
    wh <- wh + 1
    s[[wh]] <- get.stratum(dsub)
    dsub <- subset(dsub, !(tid %in% s[[wh]]$tid))
  }

  ## some checks on s, and sorting by the number of teachers
  ## (largest to smallest)
  ## sort s by the number of teachers (largest to smallest)
  nt <- sapply(s, function(x){ length(x$tid) })
  ns <- sapply(s, function(x){ length(x$sid) })
  if( (sum(nt) != length(unique(d$tid))) | (sum(ns) != length(unique(d$sid))) ){
    stop("strata error")
  }
  s <- s[order(-nt)]

  ## fix the names of s to match tid.name and sid.name and return:
  return(lapply(s, function(x){ names(x) <- c(tid.name, sid.name); x }))
}
