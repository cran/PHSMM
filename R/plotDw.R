# plot the fitted HSMM dwell-time distributions
plotDw<-function(mod,R_max,state='all',mfrow=NULL){
  nplots<-ifelse(state=='all',length(mod$mu),1)
  R_vec<-sapply(mod$p_list,length)-1
  oldmfrow<-par()$mfrow
  on.exit(par(mfrow = oldmfrow))
  if(is.null(mfrow)){
    par(mfrow=c(nplots,1))
  }else{
    par(mfrow=mfrow)
  }
  if(state=='all'){
    ptail<-list() # dwell time probabilities for r=1,..,R_max
    for(j in 1:nplots){
      pl<-mod$p_list[[j]]
      if(R_max<=R_vec[j]){
        ptail[[j]]<-pl[1:R_max]
      }else{
      ptail[[j]]<-pl[1:R_vec[j]]
      ptail[[j]]<-c(ptail[[j]],pl[R_vec[j]]*
                      (pl[R_vec[j]+1]/sum(pl[R_vec[j]+0:1]))^((R_vec[j]+1):R_max-R_vec[j]))
      }
    }
    ymax<-round(max(sapply(ptail,max))+0.01,2)
    ylabs<-substitute(d[idx](r), list(idx= 1))
    plot(1:R_max,ptail[[1]],type='h',main='state 1',xlab='dwell time r',
         ylab=ylabs,ylim=c(0,ymax),xlim=c(1,R_max))
    for(j in 2:nplots){
      ylabs<-substitute(d[idx](r), list(idx= j))
      plot(1:R_max,ptail[[j]],type='h',main=paste('state ',j,sep=''),xlab='dwell time r',
           ylab=ylabs,ylim=c(0,ymax),xlim=c(1,R_max))
    }
  }else{
    j<-state
    pl<-mod$p_list[[j]]
    if(R_max<=R_vec[j]){
      ptail[[j]]<-pl[1:R_max]
    }else{
    ptail<-pl[1:R_vec[j]]
    ptail<-c(ptail,pl[R_vec[j]]*(pl[R_vec[j]+1]/sum(pl[R_vec[j]+0:1]))^((R_vec[j]+1):R_max-R_vec[j]))
    }
    ymax<-round(max(ptail)+0.01,2)
    ylabs<-substitute(d[idx](r), list(idx= j))
    plot(1:R_max,ptail,type='h',main=paste('state ',j,sep=''),xlab='dwell time r',
         ylab=ylabs,ylim=c(0,ymax),xlim=c(1,R_max))
    }
  par(mfrow = oldmfrow)
}
