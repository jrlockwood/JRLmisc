## given a numeric vector "x", compute pairwise differences of elements
## (across lower triangle of indices)
pairdiffs <- function(x){
    stopifnot(is.numeric(x))
    .ji <- t(combn(length(x),2))
    return(data.frame(j = .ji[,1], i = .ji[,2], dif = apply(.ji, 1, function(z){ x[z[2]] - x[z[1]] })))
}
