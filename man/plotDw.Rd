\name{plotDw}
\alias{plotDw}
\title{Plot of the estimated HSMM dwell-time distributions.}
\description{
Plots the HSMM dwell-time distributions estimated using \code{\link{pmleHSMM}}.
}
\usage{
plotDw(mod, R_max, state='all', mfrow=NULL)
}
\arguments{
  \item{mod}{model object as returned by \code{pmleHSMM}.}
  \item{R_max}{integer, maximum dwell time for which the dwell-time probabilities are plotted.}
  \item{state}{value determining the states for which the distributions are plotted. Either "all" (default) for plotting the dwell-time distributions of all states, or positive integer in 1,..,N.}
  \item{mfrow}{If \code{NULL} (default) and \code{state="all"}, the probability mass functions are plotted one below the other. Otherwise, a vector of length 2 which determines the number of rows (first element) and the number of columns (second argument) of the matrix of plots.}
}
\value{
 Plot of the estimated HSMM dwell-time distributions.
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
# plot the estimated dwell-time distributions
# for dwell-times up to 12
plotDw(mod=PHSMM,R_max=12)
plotDw(mod=PHSMM,R_max=12,state=1)
plotDw(mod=PHSMM,R_max=12,mfrow=c(1,2))
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
# plot the estimated dwell-time distributions
# for dwell-times up to 15
plotDw(mod=PHSMM,R_max=15)
plotDw(mod=PHSMM,R_max=15,state=1)
plotDw(mod=PHSMM,R_max=15,mfrow=c(1,3))
}
}
