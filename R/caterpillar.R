caterpillar <- function(x,plot=TRUE,include.pdiff=TRUE,hline=0,...){
  ## creates "caterpillar" plot
  ## if x is (n x 2) uses mean (first column) and +/- 1.96 sd (second column)
  ## if x is (n x 3) uses middle column as pt estimate and first and third as limits
  x <- as.matrix(x)
  n <- dim(x)[1]
  k <- dim(x)[2]
  if(!is.numeric(x) | (!(dim(x)[2] %in% c(2,3))) ){
    stop("Invalid x")
  }
  if(k==3){
    z <- x
  } else {
    z <- cbind(x[,1] - 1.96*x[,2], x[,1], x[,1] + 1.96*x[,2])
  }
  z <- z[order(z[,2]),]
  if(plot){
    ylim <- range(z)
    plot(1:n,z[,2],type="n",ylim=ylim,...)
    for(j in 1:n){
      segments(x0=j,x1=j,y0=z[j,1],y1=z[j,3],lwd=1.0,col="darkgray")
    }
    points(1:n,z[,2])
    abline(h=hline)
  }
  x <- as.vector(apply(z,1,function(tmp){ as.numeric( (tmp[3] < 0) | (tmp[1] > 0))} ))
  pdiff <- round(sum(x)/n,digits=2)
  if(plot & include.pdiff){
    text(x=floor(0.9*n),y=ylim[1]+0.1*(ylim[2]-ylim[1]),labels=pdiff)
  }
  return(pdiff)
}
