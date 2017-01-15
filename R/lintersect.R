lintersect <- function(s){
    ## s is a list of vectors.  this returns the set of elements that are
    ## common to all members of the list
    allv <- unique(unlist(s))
    allv[sapply(allv, function(v){ all(sapply(s, function(s1){ v %in% s1 })) })]
}
