dataX<-function(N, pBM,pstrata,surrogate) {

  #if (length(strata) != length(pstrata)){
    #return("incompatible dimensions of 'strata' and 'prevalence'")}

  if(is.null(pstrata) & is.null(surrogate)){
    x1<-rbinom(N, 1, prob=pBM)
    X<-cbind(x1)}

  if (!is.null(pstrata) & is.null(surrogate)) {
    x2  <-  rbinom(N, 1, prob=pstrata[1])
    x1<-rbinom(N, 1, prob=pBM[1])
    X<-cbind(x1,x2)
  }
  if(!is.null(pstrata) &  !is.null(surrogate)){
    x2  <-  rbinom(N, 1, prob=pstrata[1])
    x1<-rbinom(N, 1, prob=pBM[1])
    sp<-rbinom(N, 1, prob=1-surrogate[2])
    se<-rbinom(N, 1, prob=surrogate[1])
    x3<-ifelse(x1==1,se,sp)
    X<-cbind(x1,x2,x3)
  }

  if(is.null(pstrata) & !is.null(surrogate)){
    x1<-rbinom(N, 1, prob=pBM[1])
    sp<-rbinom(N, 1, prob=1-surrogate[2])
    se<-rbinom(N, 1, prob=surrogate[1])
    x3<-ifelse(x1==1,se,sp)
    X<-cbind(x1,x3)}

  return(X=X)
}
