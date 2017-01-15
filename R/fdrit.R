fdrit <- function(pvec, fdr){
    ## Given the chosen false discovery rate "fdr" (e.g. 0.10), you
    ## rank all p-values from smallest (1) to largest (N).  then you reject
    ## all nulls up to the largest one "i" such that p <= (i/N)*fdr
    n <- length(pvec)
    tmp <- data.frame(p = sort(pvec), cut = fdr*(1:n)/n)
    tmp$reject <- ifelse(tmp$p <= tmp$cut, 1, 0)
    return(list(allp = tmp, ntests = n, nreject = sum(tmp$reject), fdr = fdr, pcut = max(tmp$p[tmp$reject==1])))
}
