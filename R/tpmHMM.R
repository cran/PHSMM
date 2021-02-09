# HMM transition probability matrix corresponding to HSMM
tpmHMM<-function(N,omega,d_r,R_vec,eps=1e-10){
  G<-matrix(0,0,sum(R_vec)) # initialise tpm Gamma with empty matrix
  for (i in 1:N){
    Ri<-R_vec[i]
    F<-cumsum(c(0,d_r[[i]][-Ri])) # cumulative distribution
    ci<-ifelse(abs(1-F)>eps,d_r[[i]]/(1-F),1) # elements needed for Gamma_ij
    cim<-ifelse(1-ci>0,1-ci,0)   # elements needed for Gamma_ii
    Gi<-matrix(0,Ri,0) # block of Gamma_i.
    for(j in 1:N){
      if(i==j){
        if(Ri==1){ # 1. Gamma_ii
          Gi<-cbind(Gi,c(rep(0,Ri-1),cim))
        }else{
          Gi<-cbind(Gi,rbind(cbind(rep(0,Ri-1),diag(cim[-Ri],Ri-1,Ri-1)),
                           c(rep(0,Ri-1),cim[[Ri]])))
        }
      }else{
        if(Ri==1){ # 2. Gamma_ij
          Gi<-cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,R_vec[j]-1)),1))
        }else{
          Gi<-cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,R_vec[i],R_vec[j]-1)))
        }
      }
    }
    G<-rbind(G,Gi) #combine all block elements G_i.
  }
  return(G)
}