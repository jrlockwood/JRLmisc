schoolnames.match <- function(x1,x2){
  ## x1 and x2 are assumed to be regularized versions of school names
  ## from two different sources, obtained by schoolnames.regularize.
  ## this function starts by dropping words like "elementary" and then
  ## counts intersections of the remaining words and the intersections
  ## of the soundex transformation of those words
  ##
  ## return is a data frame:
  ## n1 = # non-trivial words in x1
  ## n2 = # non-trivial words in x2
  ## nexact = # of exact matches
  ## nsoundex = # of soundex matches
  if(length(x1) != length(x2)){
    stop("x1 and x2 different length")
  }
  x1 <- strsplit(as.character(x1),"_")
  x2 <- strsplit(as.character(x2),"_")
  f <- function(x){x[ !(x %in% c("elementary","middle","high","learning","academy","community","education")) ]}
  x1 <- lapply(x1,f)
  x2 <- lapply(x2,f)
  nexact <- mapply(function(a,b){length(intersect(a,b))}, x1, x2)
  nsoundex <- mapply(function(a,b){length(intersect(soundex(a),soundex(b)))}, x1, x2)
  return(data.frame(n1 = sapply(x1,length), n2 = sapply(x2,length), nexact = nexact, nsoundex = nsoundex))
}
