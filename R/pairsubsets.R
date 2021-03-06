pairsubsets <- function(R, K){
    ## R is number of units
    ## goal is to choose subset of size K < choose(R,2) pairs without replacement to achieve as close to marginal representation as possible
    npair <- choose(R,2)
    stopifnot( (K >= 1) && (K <= npair) && (R >= 2) )

    ## fork on K
    if(K==1){
        tab <- sort(sample(1:R, size=2, replace=F))
        dim(tab) <- c(1,2)
        colnames(tab) <- c("i","j")
    } else if(K == npair){
        tab <- t(combn(R,2))
        colnames(tab) <- c("i","j")
    } else {
        nper <- (2*K) %/% R
        nrem <- (2*K) %% R
        counts <- c(rep(nper, R - nrem), rep(nper + 1, nrem))
        stopifnot(sum(counts) == 2*K)
        counts <- sample(counts) ## available counts for raters 1...R
        
        done <- FALSE
        while(!done){
            curcounts <- counts
            assg      <- vector(R, mode="list")
            
            for(i in 1:(R-1)){
                needed    <- curcounts[i]
                if(needed > 0){
                    pad       <- c(rep(0,i), curcounts[(i+1):R])
                    if(needed <= sum(pad > 0)){
                        w <- which(pad > 0)
                        ## to avoid sample() bizarre behavior with a positive integer argument
                        if(length(w)==1){
                            take <- w
                        } else {
                            take <- sort(sample(w, size=needed, replace=FALSE))
                        }
                        assg[[i]] <- take
                        curcounts[take] <- curcounts[take] - 1
                    }
                }
            }
            
            wh <- which(sapply(assg, length) > 0)
            if(length(wh) > 1){
                tab <- do.call("rbind", mapply(function(i,j){ cbind(i,j)}, as.list((1:R)[wh]), assg[wh], SIMPLIFY=FALSE))
            } else {
                tab <- cbind(i=wh,j=assg[[wh]])
            }
            done <- all(sapply(1:R, function(r){ sum(tab == r)}) == counts)        
        }
    }
    return(tab)
}
