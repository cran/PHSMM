# negative penalised log-likelihood function
npllHSMM<-function(parvect,N,y,R_vec,lambda,order_diff,
                   y_dist=c("norm","gamma","pois","bern"),stationary=TRUE,T_y,p_ref=2){
  lpn<-w2nHSMM(N,parvect,R_vec,y_dist,stationary,p_ref)
  allprobs<-matrix(1,T_y,sum(R_vec)) # condensed P-matrices needed for the log-likelihood evaluation
  cR_vec<-cumsum(c(0,R_vec))
  ind_na<-which(!is.na(y)) # filter missing values
  y_probs<-rep(1,T_y)
  f_y<-get(paste('d',y_dist,sep='')) # state-dependent distribution
  if(y_dist%in%c("norm","gamma")){
    for(j in 1:N){
      y_probs[ind_na]<-f_y(y[ind_na],lpn$mu[j],lpn$sigma[j])
      allprobs[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }else{
    for(j in 1:N){
      y_probs[ind_na]<-f_y(y[ind_na],lpn$mu[j])
      allprobs[,cR_vec[j]+1:R_vec[j]]<-y_probs
    }
  }
  if(stationary){
    delta<-lpn$delta
  }else{
    delta<-rep(0,sum(R_vec))
    delta[cR_vec[1:N]+1]<-lpn$delta
  }
  lscale=nll_Rcpp(allprobs,lpn$Gamma,delta,T_y) # forward algorithm in C++
  penalty<-sum(lambda*sapply(lpn$d_r, function(d_r,order_diff) return(sum(diff(d_r,differences=order_diff)^2)),order_diff=order_diff))
  plscale<-lscale+penalty # penalised log-likelihood value
  return(plscale)
}
