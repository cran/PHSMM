# parameter transformation
n2wHSMM<-function(N,p_list,mu,sigma=NULL,omega=NULL,delta=NULL,
                  y_dist=c("norm","gamma","pois","bern"),stationary=TRUE,p_ref=2){
  tp_list<-unlist(lapply(p_list,function(p_h) return(log(p_h[-p_ref]/p_h[p_ref]))))
  tmu<-c(mu[1],diff(mu)) # to avoid label switching
  tsigma<-NULL
  tomega<-NULL
  if(N>2){ # values only needed for N>2
    tomega<-matrix(t(omega)[!diag(N)],N,N-1,byrow=TRUE)
    tomega<-as.vector(log(tomega/tomega[,1])[,-1])
  }
  if(stationary){
    tdelta=NULL
  }else{
    tdelta<-log(delta[-1]/delta[1])
  }
  if(y_dist=='gamma'){
    tmu<-log(tmu) # log transformation - means must be positive
    tsigma<-log(sigma) # log transformation - standard deviations must be positive
    parvect<-c(tp_list,tomega,tdelta,tmu,tsigma)
  }else if(y_dist=='norm'){
    tsigma<-log(sigma) # log transformation - standard deviations must be positive
    parvect<-c(tp_list,tomega,tdelta,tmu,tsigma)
  }else if(y_dist=='pois'){
    tmu<-log(tmu) # log transformation - means must be positive
    parvect<-c(tp_list,tomega,tdelta,tmu)
  }else{
    tmu<-qlogis(mu) # logit transformation - probabilities must lie in (0,1)
    tmu<-c(tmu[1],diff(tmu))
    parvect<-c(tp_list,tomega,tdelta,tmu)
  }
  return(parvect)
}
