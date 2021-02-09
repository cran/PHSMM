# HSMM pseudo-residuals
pseudoResHSMM<- function(y,mod){
  if(!mod$y_dist%in%c("norm","gamma")){
    stop("this function requires normal or gamma distributions")
  }
  N<-length(mod$p_list)
  R_vec<-sapply(mod$p_list,length)-1
  T_y<-length(y)
  cR_vec<-cumsum(c(0,R_vec)) 
  ind_na<-which(!is.na(y))
  if(mod$stationary){
    delta<-mod$delta
  }else{
    delta<-rep(0,sum(R_vec))
    delta[cR_vec[1:N]+1]<-mod$delta
  }
  #1) calculate forward variables
  lalpha <- matrix (NA,sum(R_vec),T_y)
  allprobs <- matrix(1,T_y,sum(R_vec))
  f_y<-get(paste('d',mod$y_dist,sep=''))
  y_probs<-rep(1,T_y)
  if(mod$y_dist=="gamma"){
    mod$shape=mod$mu^2/mod$sigma^2
    mod$rate=mod$mu/mod$sigma^2
    for(j in 1:N){
      y_probs[ind_na]<-f_y(y[ind_na],mod$shape[j],mod$rate[j])
      allprobs[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }else{
    for (j in 1:N){
      y_probs[ind_na]<-f_y(y[ind_na],mod$mu[j],mod$sigma[j])
      allprobs[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }
  # forward algorithm, scaling to avoid numerical underflow
  foo<-delta*allprobs[1,]
  sumfoo<-sum(foo)
  lscale<-log(sumfoo)
  foo<-foo/sumfoo
  lalpha[,1]<-lscale+log(foo)
  for (t in 2:T_y){
    foo<-foo%*%mod$Gamma*allprobs[t,] 
    sumfoo<-sum(foo)
    lscale<-lscale+log(sumfoo)
    foo<-foo/sumfoo
    lalpha[,t]<-log(foo)+lscale
  }
  # 2) calculate one-step ahead forecast pseudo residuals
  # scaling to avoid numerical underflow
  Res<-rep(NA,T_y)
  F_mat<-matrix(NA,T_y,sum(R_vec))
  y_probs<-rep(NA,T_y)
  F_y<-get(paste('p',mod$y_dist,sep=''))
  if(mod$y_dist=="gamma"){
    for (j in 1:N){
      y_probs[ind_na]<-F_y(y[ind_na],mod$shape[j],mod$rate[j])
      F_mat[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }else{
    for (j in 1:N){
      y_probs[ind_na]<-F_y(y[ind_na],mod$mu[j],mod$sigma[j])
      F_mat[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }
  Res[1]<-qnorm(t(delta)%*%F_mat[1,]) 
  for(i in 2:T_y){
    c<-max(lalpha[,i-1])
    a<-exp(lalpha[,i-1]-c) 
    Res[i]<-qnorm(t(a)%*%(mod$Gamma/sum(a))%*%F_mat[i,])
  }
  return(Res)
}
