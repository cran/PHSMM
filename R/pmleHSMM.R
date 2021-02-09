# HSMM penalised maximum likelihood estimation
pmleHSMM<-function(y,N,p_list,mu,sigma=NULL,omega=NULL,delta=NULL,lambda,
                   order_diff,y_dist=c("norm","gamma","pois","bern"),
                   stationary=TRUE,p_ref=2,print.level = 0,iterlim = 10000,
                   stepmax = NULL,hessian=FALSE,gradtol=10^(-6)){
  # check parameters
  if(!y_dist%in%c("norm","gamma","pois","bern")){
    stop('state-dependent distribution is not correctly specified')
  }
  if(!length(p_list)==N|!length(mu)==N){
    stop("number of parameters is wrong")
  }
  if(y_dist%in%c("gamma","pois")){
    if(sum(!mu>0)>0){
      stop('mu must be strictly positive')
    }
  }
  if(sum(!order(mu)==(1:N))>0){
    stop('mu must be sorted in an ascending order')
  }
  if(y_dist%in%c("norm","gamma")){
    if(!length(sigma)==N){
      stop("number of parameters is wrong")
    }
    if(sum(!sigma>0)>0){
      stop('sigma must be strictly positive')
    }
  }
  if(y_dist=="bern"){
    if(sum(mu<0|mu>1)>0){
      stop('elements in mu must lie between zero and one')
    }
  }
  if(stationary==FALSE & is.null(delta)){
    stop('delta is missing; set stationary=TRUE or define an initial delta')
  }
  if(N>2){
    if(is.null(omega)){
      stop('omega is missing')
    }
    if(!nrow(omega)==N|!ncol(omega)==N){
      stop('dimensions of omega are wrong; omega must be an NxN matrix')
    }
  }
  if(!length(lambda)==N){
    stop('length of lambda is wrong')
  }
  R_vec<-sapply(p_list,length)-1 # lengths of the unstructured starts
  T_y<-length(y) # length of observed time series
  parvect<-n2wHSMM(N,p_list,mu,sigma,omega,delta,y_dist,stationary,p_ref)
  stepmax<-ifelse(!is.null(stepmax),stepmax,max(1000*sqrt(sum((parvect/rep(1, length(parvect)))^2)), 1000))
  mod<-withCallingHandlers({nlm(npllHSMM,parvect,N,y,R_vec,lambda,order_diff,y_dist,stationary,T_y,p_ref,
                                print.level = print.level,iterlim = iterlim,stepmax = stepmax,gradtol=gradtol,hessian=hessian)},
                           warning=function(w) {
                             if(any(grepl("NA/Inf replaced by maximum positive value",w))){
                               invokeRestart("muffleWarning")
                             }
                           }) # optimisation using nlm, filters nlm warning "NA/Inf replaced by maximum positive value"
  pn<-w2nHSMM(N,mod$estimate,R_vec,y_dist,stationary,p_ref)
  if(y_dist%in%c("pois","bern")){
    if(hessian){
      return(list(mu=pn$mu,p_list=pn$p_list,delta=pn$delta,omega=pn$omega,Gamma=pn$Gamma,npll=mod$minimum,gradient=mod$grad,iterations=mod$iterations,code_conv=mod$code,w_estimates=mod$estimate,stationary=stationary,y_dist=y_dist))
    }else{
      return(list(mu=pn$mu,p_list=pn$p_list,delta=pn$delta,omega=pn$omega,Gamma=pn$Gamma,npll=mod$minimum,gradient=mod$grad,iterations=mod$iterations,code_conv=mod$code,w_estimates=mod$estimate,stationary=stationary,y_dist=y_dist))
    }
  }else{
    if(y_dist=="gamma"){
      pn$mu<-pn$mu2 # mean value
      pn$sigma<-pn$sigma2 # standard deviation
    }
    if(hessian){
      return(list(mu=pn$mu,sigma=pn$sigma,p_list=pn$p_list,delta=pn$delta,omega=pn$omega,Gamma=pn$Gamma,npll=mod$minimum,gradient=mod$grad,iterations=mod$iterations,code_conv=mod$code,w_estimates=mod$estimate,stationary=stationary,y_dist=y_dist))
    }else{
      return(list(mu=pn$mu,sigma=pn$sigma,p_list=pn$p_list,delta=pn$delta,omega=pn$omega,Gamma=pn$Gamma,npll=mod$minimum,gradient=mod$grad,iterations=mod$iterations,code_conv=mod$code,w_estimates=mod$estimate,stationary=stationary,y_dist=y_dist))
    }
  }
}