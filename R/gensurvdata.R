
genSAMPLEexpWeiC<-function(N, lambda, k, betaBM, betastrata, X, cens, tau)
{
  u <- runif(N, min=0, max=1)

  if("x3" %in% colnames(X)) {
    beta<-c(betaBM, betastrata,0)
  }else{
    beta<-c(betaBM, betastrata)}

  #if (k=="1")
    #genecens onset time:
    #t<- -1/(lambda*exp(X %*% beta))*log(1-u)
  #t<- -1/(lambda*exp(x1*beta+x2*beta2+x3*beta3+x4*beta4))*log(1-u)
  #if (k!="1")
    t <- (-(log(u))/(lambda*exp(X %*% beta)))^(1/k)

  if(cens==0){
   cens.time<-rep(tau, N)
  }
  if(cens!=0){
  cens.time <- rexp(length(t), rate = cens)
  }
  cens.time <- ifelse(cens.time > tau, tau, cens.time)
  censor <- ifelse (t > cens.time, 0 , 1)

  time <- ifelse(t < cens.time, t, cens.time)
  id<-1:N
  return(data.frame(id, X, time, censor))
}
