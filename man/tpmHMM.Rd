\name{tpmHMM}
\alias{tpmHMM}
\title{Tranistion probability matrix of the HMM representation}
\description{
  Construction of the transition probability matrix corresponding to the HMM which exactly represents the HSMM. Not intended to be run by the user (internal function, called by the functions \code{\link{w2nHSMM}} and \code{\link{pmleHSMM}}).
}
\usage{
tpmHMM(N, omega, d_r, R_vec, eps=1e-10)
}
\arguments{
  \item{N}{number of states of the HSMM, integer greater than 1.}
  \item{omega}{conditional transition probability matrix of the HSMM with N rows and N columns. The diagonal elements must be zero and the rows must sum to one.}
  \item{d_r}{list of vectors containing the dwell-time probabilities of the unstructured starts.}
  \item{R_vec}{vector of length N containing the lengths of the unstructured starts of the state dwell-time distributions.}
  \item{eps}{to avoid negative probabilities due to numerical underflow. Default is \code{1e-10}.}
}
\value{
  Returns the transition probability matrix of the HMM which exactly represents the HSMM.
}
\references{
For details on the code and construction of the matrix, see:

Langrock, R. and Zucchini W. (2011): Hidden Markov models with arbitrary state dwell-time distributions. \emph{Computational Statistics and Data Analysis}, \strong{55}, p. 715--724.

Zucchini, W., MacDonald, I.L. and Langrock, R. (2016): \emph{Hidden Markov models for time series:
An introduction using R}. 2nd edition. Chapman & Hall/CRC, Boca Raton.
}
\examples{
# list of dwell-time probability vectors 
# (vector elements should not sum to one)
d_r<-list()
d_r[[1]]<-c(dgeom(0:9,0.2))
d_r[[2]]<-c(dgeom(0:9,0.1))
# tranistion probability matrix:
Gamma<-tpmHMM(N=2,omega=matrix(c(0,1,1,0),2,2),
  d_r=d_r,R_vec=sapply(d_r,length))
}

