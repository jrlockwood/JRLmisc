stringpad <- function(x, pad="0", trailing=FALSE){
    ## take a vector of character strings and makes each element the same nchar by padding with leading "pad".
    ## alternatively, puts "trailing" pad if trailing==TRUE
    if( nchar(pad) > 1 ){ stop("illegal pad") }
    if( !(is.vector(x) && is.character(x)) ){ stop("x must be a character vector") }
    .mis <- which(is.na(x))
    .obs <- which(!is.na(x))
    x[.mis] <- ""
    tot     <- max(nchar(x))
    if(trailing){
        res <- as.vector( sapply(x, function(z){ paste(z,paste(rep(pad,times=tot-nchar(z)),collapse=""),sep="") }) )
    } else {
        res <- as.vector( sapply(x, function(z){ paste(paste(rep(pad,times=tot-nchar(z)),collapse=""),z,sep="") }) )
    }
    is.na(res[.mis]) <- TRUE
    return(res)
}
