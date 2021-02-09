\name{n2wHSMM}
\alias{n2wHSMM}
\title{Parameter transformation from natural to working parameters}
\description{
  Parameter transformation from the natural (constraint) HSMM parameters into unconstraint working parameters which are used in the numerical maximum likelihood estimation. Not intended to be run by the user (internal function, called by the function \code{\link{pmleHSMM}}).
}
\usage{
n2wHSMM(N, p_list, mu, sigma=NULL, omega=NULL, delta=NULL, 
      y_dist=c("norm","gamma","pois","bern"), stationary=TRUE, p_ref=2)
}
\arguments{
  \item{N}{number of states of the HSMM, integer greater than 1.}
  \item{p_list}{list containing the parameters of the states' dwell-time distributions. The list consists of \code{N} probability vectors, i.e. one vector for each state. Each probability vector contains the state dwell-time probabilities for the unstructured start and, as last element, the probability mass captured in the geometric tail. Thus, each vector must sum to one and automatically determines the length of the unstructured start of the according state dwell-time distribution.}
  \item{mu}{vector of length \code{N} containing the state-dependent mean values if gamma, normal or Poisson distributions are chosen and the state-dependent probabilities if the Bernoulli distribution is chosen. The values must be sorted in an ascending order.}
  \item{sigma}{vector of length \code{N} containing the state-dependent standard deviations if gamma or normal distributions are chosen. \code{NULL} (default) otherwise.}
  \item{omega}{conditional transition probability matrix of the underlying semi-Markov chain. Only needed if the number of states is greater than 2, \code{NULL} (default) otherwise. In the former case, \code{omega} is a matrix with \code{N} rows and \code{N} columns, its diagonal elements must be zero and its rows must sum to one.}
  \item{delta}{vector of length N containing the initial distribution. Only needed if \code{stationary=FALSE}, otherwise delta is ignored and can be set to \code{NULL} (default).}
  \item{y_dist}{character determining the class of state-dependent distributions. Supported values are \code{"norm"} (normal distribution), \code{"gamma"} (gamma distribution), \code{"pois"} (Poisson distribution) and \code{"bern"} (Bernoulli distribution).}
  \item{stationary}{Logical, if \code{TRUE} (default), stationarity is assumed, if \code{FALSE}, the underlying state-sequence is assumed to enter a new state at time t=1 and it is necessary to define the initial distribution \code{delta}.}
  \item{p_ref}{positive integer determining the reference dwell-time probability used for the multinomial logit parameter transformation. Default value is 2. Only needs to be changed if the dwell-time probability for dwell time r=2 is estimated very close to zero in order to avoid numerical problems.}
}
\value{
  A vector of unconstraint working parameters characterising the HSMM.
}
\details{
  The transformation from natural to working parameters is needed to carry out an unconstraint optimisation. The function includes log-transformations for positive parameters and (multinomial) logit-transformations for probabilities, probability vectors and matrices.
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
n2wHSMM(N=2,p_list=p_list0,mu=mu0,sigma=sigma0,y_dist='norm',stationary=TRUE)
}

