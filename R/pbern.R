# cumulative distribution function of the Bernoulli distribution
pbern<-function(y,prob){
  pbinom(y,size=1,prob=prob,log.p=FALSE)
}
