\name{nll_Rcpp}
\alias{nll_Rcpp}
\title{Negative log-likelihood evaluation}
\description{Evaluation of the negative HMM log-likelihood function based on the forward algorithm, written in \code{C++} (internal function, called by the function \code{\link{npllHSMM}}).
}
\usage{
nll_Rcpp(allprobs, gamma, delta, T_y)
}
\arguments{
  \item{allprobs}{matrix containing the state-dependent distribution values for each observation (row) and state (column), respectively.}
  \item{gamma}{transition probability matrix.}
  \item{delta}{initial distribution.}
  \item{T_y}{number of observations.}
}
\value{
  Returns the negative log-likelihood value.
}
