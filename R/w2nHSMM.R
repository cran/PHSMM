# parameter back-transformation
w2nHSMM<-function(N,parvect,R_vec,y_dist=c("norm","gamma","pois","bern"),
                  stationary=TRUE,p_ref=2){
  cR_vec<-cumsum(c(0,R_vec))
  sR_vec<-sum(R_vec)
  p_list<-lapply(1:N,function(i){
    p_h<-numeric(R_vec[i]+1)
    p_h[-p_ref]<-parvect[cR_vec[i]+1:R_vec[i]]
    p_h<-c(exp(p_h))/(sum(exp(p_h)))
    return(p_h)
  })
  d_r<-lapply(p_list,function(p_h) return(p_h[-length(p_h)]))
  foo<-sR_vec
  if(N>2){ # only needed if N>2
    omega<-matrix(0,N,N)
    omega[!diag(N)]<-as.vector(t(matrix(c(rep(1,N),exp(parvect[foo+1:(N*(N-2))])),N,N-1)))
    omega<-t(omega)/apply(omega,2,sum)
    foo<-foo+(N*(N-2))
  }else{
    omega<-matrix(c(0,1,1,0),2,2)
  }
  Gamma<-tpmHMM(N,omega,d_r,R_vec) # tpm of the HMM representing the HSMM
  if(stationary){
    delta<-solve(t(diag(sR_vec)-Gamma+1),rep(1,sR_vec))
  }else{
    delta<-c(1,exp(parvect[foo+1:(N-1)]))
    delta<-delta/sum(delta)
    foo<-foo+N-1
  }
  mu2<-NULL
  sigma2<-NULL
  if(y_dist=="norm"){
    mu<-cumsum(parvect[foo+1:N])
    sigma<-exp(parvect[foo+N+1:N])
  }else if(y_dist=="gamma"){
    mu2<-cumsum(exp(parvect[foo+1:N]))
    sigma2<-exp(parvect[foo+N+1:N])
    mu<-mu2^2/sigma2^2 # shape parameter, needed for dgamma
    sigma<-mu2/sigma2^2 #rate parameter, needed for dgamma
  }else if(y_dist=='pois'){
    mu<-cumsum(exp(parvect[foo+1:N]))
    sigma<-NULL
  }else{
    mu<-plogis(cumsum(parvect[foo+1:N]))
    sigma<-NULL
  }
  return(list(p_list=p_list,d_r=d_r,omega=omega,Gamma=Gamma,delta=delta,mu=mu,sigma=sigma,mu2=mu2,sigma2=sigma2))
}
