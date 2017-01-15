dfcmp<-function(df1,df2,detailed=FALSE){
    ## report of field by field comparison of dataframes
    if(!is.data.frame(df1)){
        stop("df1 must be a dataframe")
    }

    if(!is.data.frame(df2)){
        stop("df1 must be a dataframe")
    }

    equiv<-function(x,y){
        ## a less strict version of identical()
        if(length(x)!=length(y)){
            stop("error in equiv(): Lengths of x and y differ")
        }
        x.miss<-is.na(x)
        y.miss<-is.na(y)
        x.nmis<-sum(x.miss)
        y.nmis<-sum(y.miss)
        if(x.nmis!=y.nmis){
            return(list(x.nmis=x.nmis,y.nmis=y.nmis,agree.mis=FALSE,agree.obs=FALSE))
        } else {
            agree.mis<-all(which(x.miss)==which(y.miss))
            agree.obs<-all(x[!x.miss]==y[!y.miss])
            return(list(x.nmis=x.nmis,y.nmis=y.nmis,agree.mis=agree.mis,agree.obs=agree.obs))
        }
    }
    
    names.df1<-names(df1)
    names.df2<-names(df2)
    df1.n<-deparse(substitute(df1))
    df2.n<-deparse(substitute(df2))
    print(paste("Dimension of",df1.n,":",deparse(dim(df1))))
    print(paste("Dimension of",df2.n,":",deparse(dim(df2))))
    f.df1.notdf2<-names.df1[!(names.df1 %in% names.df2)]
    if(length(f.df1.notdf2)==0){
        f.df1.notdf2<-"NONE"
    }
    print(paste("Fields in",df1.n,"not in",df2.n,":"))
    print(f.df1.notdf2)
    f.df2.notdf1<-names.df2[!(names.df2 %in% names.df1)]
    if(length(f.df2.notdf1)==0){
        f.df2.notdf1<-"NONE"
    }
    print(paste("Fields in",df2.n,"not in",df1.n,":"))
    print(f.df2.notdf1)
    v<-c(names.df1,names.df2)
    v<-v[duplicated(v)]
    disagree<-character(0)
    for(i in 1:length(v)){
        ## Need to deal with the problem that if one database
        ## was subsetted from a larger one, then factor levels
        ## may not match up.  Thus convert factors to characters.
        if( (is.factor(df1[,v[i]])) & (is.factor(df2[,v[i]])) & (length(levels(df1[,v[i]]))!=length(levels(df2[,v[i]]))) ){
            tmp1<-as.character(df1[,v[i]])
            tmp2<-as.character(df2[,v[i]])
            res<-equiv(tmp1,tmp2)
        } else{
            res<-equiv(df1[,v[i]],df2[,v[i]])
        }
        if(!((res$agree.mis)&(res$agree.obs))) disagree<-c(disagree,v[i])
        if(detailed) print(c(v[i],res$agree.obs,res$agree.mis))
    }
    print("Common fields where they disagree:")
    print(disagree)
}
