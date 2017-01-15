haversine <- function(locs,units="Statute Miles"){
    ## Author: J.R. Lockwood, RAND
    ## Last Modified 2/27/03
    ##
    ## Calculates matrix of pairwise distances between long/lat
    ## pairs specified in locs using Haversine Great Circle distance.
    ## locs is a (nloc X 2) matrix of the form
    ##
    ## (node 1 longitude) (node 1 latitude)
    ## (node 2 longitude) (node 2 latitude)
    ## ...
    ## Result is a symmetric matrix whose (i,j) element is the distance in
    ## miles between node i and node j
    ## Radii taken from http://jan.ucc.nau.edu/~cvm/latlongdist.html
    nloc<-dim(locs)[1]
    if(nloc<2){
        stop("Need at least two locations")
    }
    if(units=="Statute Miles"){
        R<-3963.1 ## in the past I used 3956
    } else if(units=="Nautical Miles"){
        R<-3443.9
    } else if(units=="Kilometers"){
        R<-6378
    } else{
        stop("Units must be \"Statute Miles\", \"Nautical Miles\", or \"Kilometers\"")
    }
    dists<-matrix(0,ncol=nloc,nrow=nloc)
    for(i in 1:nloc){
        for(j in 1:i){
            dlon<-(locs[j,1]-locs[i,1])*pi/180
            dlat<-(locs[j,2]-locs[i,2])*pi/180
            a<-sin(0.5*dlat)*sin(0.5*dlat)+cos(locs[i,2]*pi/180)*cos(locs[j,2]*pi/180)*sin(0.5*dlon)*sin(0.5*dlon)
            if(sqrt(a)>=1) c<-pi
            else c<-2*asin(sqrt(a))
            dists[i,j]<-R*c
        }
    }
    return(dists+t(dists))
}

