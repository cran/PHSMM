\name{w2nHSMM}
\alias{w2nHSMM}
\title{Parameter transformation from working to natural parameters}
\description{
  Transforms unconstraint HSMM working parameters back into (constraint) natural parameters. Not intended to be run by the user (internal function, called by the functions \code{\link{pmleHSMM}} and \code{\link{npllHSMM}}).
}
\usage{
w2nHSMM(N, parvect, R_vec, y_dist=c("norm","gamma","pois","bern"),
      stationary=TRUE, p_ref=2)
}
\arguments{
  \item{N}{number of states of the HSMM, integer greater than 1.}
  \item{parvect}{vector of unconstraint working parameter as obtained by the function \code{\link{n2wHSMM}}.}
  \item{R_vec}{vector of length \code{N} containing the lengths of the unstructured starts of the dwell-time distributions.}
  \item{y_dist}{character determining the class of state-dependent distributions used to model the observations. Supported values are \code{"norm"} (normal distribution), \code{"gamma"} (gamma distribution), \code{"pois"} (Poisson distribution) and \code{"bern"} (Bernoulli distribution).}
  \item{stationary}{Logical, if \code{TRUE} (default), stationarity is assumed, if \code{FALSE}, the underlying state-sequence is assumed to enter a new state at time t=1.}
  \item{p_ref}{positive integer determining the reference dwell-time probability used for the multinomial logit parameter transformation. Default value is 2. Only needs to be changed if the dwell-time probability for dwell time r=2 is estimated very close to zero in order to avoid numerical problems.}
}
\value{
  A list containing the natural parameters
  \item{p_list}{list containing the dwell-time distribution vectors for each state. Each of the \code{N} vectors contains the state dwell-time probabilities for the unstructured start and, as last element, the probability mass captured in the geometric tail. Thus, each vector sums to one.}
  \item{mu}{vector of length \code{N}. For \code{y_dist="norm"} and \code{y_dist="pois"}, it contains the state-dependent mean values, for \code{y_dist="gamma"}, it contains the state-dependent shape parameters, and for \code{y_dist="bern"}, it contains the state-dependent probabilities.}
  \item{sigma}{vector of length \code{N} containing the state-dependent standard deviations if \code{y_dist="norm"} and the state-dependent rates if \code{y_dist="gamma"}. \code{NULL} otherwise.}
  \item{omega}{conditional transition probability matrix of the HSMM.}
  \item{delta}{equilibrium distribution if \code{stationary=TRUE}, initial distribution of length N if \code{stationary=FALSE}.}
  \item{d_r}{list containing the dwell-time probabilities of the unstructured starts.}
  \item{Gamma}{transition probability matrix of the HMM which represents the HSMM.}
}
\details{
  The function reverses the transformation of the function \code{\link{n2wHSMM}} and back-transforms the unconstraint parameters into the constraint natural parameters. Note that if \code{y_dist="gamma"}, \code{mu} and \code{sigma} do not include the mean values and standard deviations, but the shape and rate parameters as required by the density functions \code{\link{dgamma}} and \code{\link{pgamma}}. The mean and standard deviations are then assigned to \code{mu2} and \code{sigma2}.
}
\examples{
# natural parameters for 2-state HSMM with state-dependent normal distributions
p_list0<-list() # list of dwell-time distribution vectors,
                # vector elements must sum to one
p_list0[[1]]<-c(dgeom(0:9,0.2),1-pgeom(9,0.2))
p_list0[[2]]<-c(dgeom(0:9,0.1),1-pgeom(9,0.1))
mu0<-c(-10,10) # mean values
sigma0<-c(3,5) # standard deviations
# parameter transformation:
parvect<-n2wHSMM(N=2,p_list=p_list0,mu=mu0,sigma=sigma0,y_dist='norm',stationary=TRUE)
# back-transformation:
npar<-w2nHSMM(N=2,parvect=parvect,R_vec=sapply(p_list0,length)-1,y_dist='norm')
}
