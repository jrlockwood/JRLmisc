roundto<-function(x,tab){
    ## rounds/projects the elements of x to the nearest elements
    ## found in tab.
    if(any(is.na(tab))){
        stop("Error in roundto: tab has missing values")
    }
    n<-length(x)
    res<-rep(NA,n)
    for(i in 1:n){
        if(!is.na(x[i])){
            res[i]<-tab[which.min(abs(x[i]-tab))]
        }
    }
    res
}
