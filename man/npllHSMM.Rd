\name{npllHSMM}
\alias{npllHSMM}
\title{negative penalised HSMM log-likelihood function}
\description{
  Evaluates the negative penalised log-likelihood function of the HSMM (internal function, called by the function \code{\link{pmleHSMM}}).
}
\usage{
npllHSMM(parvect ,N, y, R_vec, lambda, order_diff,
        y_dist=c("norm","gamma","pois","bern"),
        stationary=TRUE, T_y,  p_ref=2)
}
\arguments{
  \item{parvect}{vector of unconstraint working parameter as returned by the function \code{\link{n2wHSMM}}.}
  \item{N}{number of states of the HSMM, integer greater than 1.}
  \item{y}{vector containing the observed time series.}
  \item{R_vec}{vector of length \code{N} containing the lengths of the unstructured starts of the dwell-time distributions.}
  \item{lambda}{vector of length \code{N} containing the smoothing parameter values to weight the penalty term.}
  \item{order_diff}{order of the differences used for the penalty term, positive integer which does not exceed the length of the unstructured starts.}
  \item{y_dist}{character determining the class of state-dependent distributions. Supported values are \code{"norm"} (normal distribution), \code{"gamma"} (gamma distribution), \code{"pois"} (Poisson distribution) and \code{"bern"} (Bernoulli distribution).}
  \item{stationary}{Logical, if \code{TRUE} (default), stationarity is assumed, if \code{FALSE}, the underlying state-sequence is assumed to enter a new state at time t=1.}
  \item{T_y}{length of the observed time series.}
  \item{p_ref}{positive integer determining the reference dwell-time probability used for the multinomial logit parameter transformation. Default value is 2. Only needs to be changed if the dwell-time probability for dwell time r=2 is estimated very close to zero in order to avoid numerical problems.}
}
\value{
  Returns the value of the negative penalised HSMM log-likelihood function for the given parameters and time series.
}
\details{
The penalised log-likelihood function relies on the exact HMM representation of the HSMM and is evaluated using the forward algorithm which is implemented in \code{C++} to speed up the calculation.
}
\references{
Pohle, J., Adam, T. and Beumer, L.T. (2021): Flexible estimation of the state dwell-time distribution in hidden semi-Markov models. arXiv:\url{https://arxiv.org/abs/2101.09197}.

Zucchini, W., MacDonald, I.L. and Langrock, R. (2016): \emph{Hidden Markov models for time series:
An introduction using R}. 2nd edition. Chapman & Hall/CRC, Boca Raton.
}
\examples{
# 3-state gamma HSMM and hourly muskox step length
# natural parameters
p_list0<-list()
p_list0[[1]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
p_list0[[2]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
p_list0[[3]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
omega0<-matrix(0.5,3,3)
diag(omega0)<-0
mu0<-c(5,100,350)
sigma0<-c(3,90,300)
R_vec<-sapply(p_list0,length)-1 # lengths of the unstructured starts
# working parameter vector
parvect<-n2wHSMM(N=3,p_list=p_list0,mu=mu0,sigma=sigma0,
  omega=omega0,y_dist='gamma')
# evaluate the negative penalised log-likelihood function
npllHSMM(parvect,N=3,muskox$step,R_vec=R_vec,lambda=c(1000,1000,1000),
  order_diff=2,y_dist="gamma",T_y=nrow(muskox))
}

