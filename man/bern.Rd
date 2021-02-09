\name{Bernoulli distribution}
\alias{dbern}
\alias{pbern}
\title{Bernoulli distribution}
\description{
  Probability mass function and cumulative distribution function of the Bernoulli distribution.
}
\usage{
dbern(y, prob)
pbern(y, prob)
}
\arguments{
  \item{y}{a vector of zeros and ones.}
  \item{prob}{probability.}
}
\value{
 \code{dbern} returns the probability mass function, \code{pbern} returns the cumulative distribution function.
}
\details{
  The code relies on the functions \code{\link{dbinom}} and \code{\link{pbinom}} with \code{size=1} and \code{log=FALSE}.
}
\references{
  See the documentation for \code{\link{dbinom}} and \code{\link{pbinom}} for more details.
}
\examples{
dbern(c(0,1),0.4)
pbern(c(0,1),0.4)
}
