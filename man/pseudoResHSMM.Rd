\name{pseudoResHSMM}
\alias{pseudoResHSMM}
\title{HSMM pseudo-residuals}
\description{
Pseudo-residuals based on the one-step ahead forecast distributions under the HSMM which was estimated using \code{\link{pmleHSMM}}. This function can only be used for HSMMs with state-dependent normal or gamma distributions. 
}
\usage{
pseudoResHSMM(y, mod)
}
\arguments{
  \item{y}{vector containing the observations.}
  \item{mod}{model object as returned by \code{pmleHSMM}.}
}
\value{
 Returns a vector containing the forecast pseudo-residuals.
}
\details{
A good model fit is indicated by standard normally distributed pseudo-residuals. 
}
\references{
For more details about pseudo-residuals in the context of HMMs, see:

Zucchini, W., MacDonald, I.L. and Langrock, R. (2016): \emph{Hidden Markov models for time series:
An introduction using R}. 2nd edition. Chapman & Hall/CRC, Boca Raton.
}
\examples{
\donttest{
# running this example might take a few minutes
#
# 1.) 2-state gamma-HSMM for hourly muskox step length
# with an unstructured start of length of 10
#
# initial values
p_list0<-list()
p_list0[[1]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
p_list0[[2]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
mu0<-c(5,150)
sigma0<-c(3,180)
#
# fit 2-state gamma-HSMM with lambda=c(100,100)
# and difference order 3
# estimation might take a few minutes
PHSMM<-pmleHSMM(y=muskox$step,N=2,p_list=p_list0,mu=mu0,
                sigma=sigma0,lambda=c(100,100),order_diff=3,
                y_dist='gamma')
#
# pseudo-residuals
pseudoRes<-pseudoResHSMM(y=muskox$step,PHSMM)
hist(pseudoRes,probability=TRUE)
z<-seq(-3,3,0.01)
lines(z,dnorm(z),col='blue')
}


\donttest{
# running this example might take a few minutes
#
# 2.) 3-state gamma-HSMM for hourly muskox step length
# with an unstructured start of length of 10
#
# initial values
p_list0<-list()
p_list0[[1]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
p_list0[[2]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
p_list0[[3]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
omega0<-matrix(0.5,3,3)
diag(omega0)<-0
mu0<-c(5,100,350)
sigma0<-c(3,90,300)
#
# fit 3-state gamma-HSMM with lambda=c(1000,1000,1000)
# and difference order 3
# estimation might take some minutes
PHSMM<-pmleHSMM(y=muskox$step,N=3,p_list=p_list0,mu=mu0,
                sigma=sigma0,omega=omega0,
                lambda=c(1000,1000,1000),
                order_diff=3,y_dist='gamma')
#
# pseudo-residuals
pseudoRes<-pseudoResHSMM(y=muskox$step,PHSMM)
hist(pseudoRes,probability=TRUE)
z<-seq(-3,3,0.01)
lines(z,dnorm(z),col='blue')
}
}
