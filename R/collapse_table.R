collapse_table <- function(tab, K){
    stopifnot(all(!is.na(tab)) && is.vector(tab) && is.numeric(tab) && all(tab >= 0) && !is.null(names(tab)) && (length(tab) > K) && (K >= 2))
    newtab   <- tab
    L        <- length(newtab)
    
    while(L > K){
        .cands <- lapply(1:(L-1), function(i){
            if(i==1){
                .tmp <- c(sum(newtab[1:2]), newtab[3:L])
                names(.tmp) <- c(paste(names(newtab)[1:2],collapse=""), names(newtab)[3:L])
            } else if(i==(L-1)){
                .tmp <- c(newtab[1:(L-2)], sum(newtab[(L-1):L]))
                names(.tmp) <- c(names(newtab)[1:(L-2)], paste(names(newtab)[(L-1):L],collapse=""))
            } else {
                .tmp <- c(newtab[1:(i-1)], sum(newtab[i:(i+1)]), newtab[(i+2):L])
                names(.tmp) <- c(names(newtab)[1:(i-1)], paste(names(newtab)[i:(i+1)],collapse=""), names(newtab)[(i+2):L])
            }
            return(.tmp)
        })
        newtab <- .cands[[which.min(sapply(.cands, sd))]]
        L      <- L-1
    }
    return(newtab)
}

## ## try it
## tab <- as.vector(table(cut(rnorm(10000), breaks=20)))
## names(tab) <- paste0("c",1:length(tab))
## print(collapse_table(tab, K=10))
## print(collapse_table(tab, K=3))
## print(collapse_table(tab, K=2))
## 
## tab <- as.vector(table(cut(runif(10000), breaks=20)))
## names(tab) <- paste0("c",1:length(tab))
## print(collapse_table(tab, K=10))



