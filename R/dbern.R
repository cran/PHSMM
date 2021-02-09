# probability mass function of the Bernoulli distribution
dbern<-function(y,prob){
  dbinom(y,size=1,prob=prob,log=FALSE)
}

