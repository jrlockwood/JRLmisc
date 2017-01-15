lex <- function(N){
    ## produce vector of N lexicograpically ordered strings
    if(N<=0){
        stop("invalid N")
    }
    ndig <- nchar(N)
    substr(formatC((1:N)/10^ndig,digits=ndig,format="f"),3,10000000)
}
