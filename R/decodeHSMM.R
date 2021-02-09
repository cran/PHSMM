# HSMM state decoding
decodeHSMM<-function(y,mod){
  N<-length(mod$p_list)
  R_vec<-sapply(mod$p_list,length)-1
  T_y<-length(y)
  allprobs<-matrix(1,T_y,sum(R_vec))
  cR_vec<-cumsum(c(0,R_vec))
  ind_na<-which(!is.na(y)) # filter missing values
  y_probs<-rep(1,T_y)
  f_y<-get(paste("d",mod$y_dist,sep='')) # state-dependent distributions
  if(mod$y_dist=="norm"){
    for(j in 1:N){
      y_probs[ind_na]<-f_y(y[ind_na],mod$mu[j],mod$sigma[j])
      allprobs[,cR_vec[j]+1:R_vec[j]] <- y_probs
    }
  }else if(mod$y_dist=="gamma"){
    shape<-mod$mu^2/mod$sigma^2
    rate<-mod$mu/mod$sigma^2
    for(j in 1:N){
      y_probs[ind_na]<-f_y(y[ind_na],shape[j],rate[j])
      allprobs[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }else{
    for(j in 1:N){
      y_probs[ind_na]<-f_y(y[ind_na],mod$mu[j])
      allprobs[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }
  if(mod$stationary){
    delta<-mod$delta
  }else{
    delta<-rep(0,sum(R_vec))
    delta[cR_vec[1:N]+1]<-mod$delta
  }
  # traverse forward in time 
  y_t<-matrix(0,T_y,sum(R_vec))
  foo<-delta*allprobs[1,]
  y_t[1,]<-foo/sum(foo)
  for(t in 2:T_y){
    foo<-apply(y_t[t-1,]*mod$Gamma,2,max)*allprobs[t,]
    y_t[t,]<-foo/sum(foo)
  }
  # traverse backward in time 
  s_t<-numeric(T_y)
  s_t[T_y]<-which.max(y_t[T_y,])
  for(t in (T_y-1):1){
    s_t[t]<-which.max(mod$Gamma[,s_t[t+1]]*y_t[t,])
  }
  sHSMM<-rep(1,T_y)
  for(j in 2:N){
    sHSMM[s_t>cR_vec[j] & s_t<(cR_vec[j+1]+1)]<-j
  }
  return(sHSMM)
}  
