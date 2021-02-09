\docType{data}
\name{muskox}
\alias{muskox}
\title{Muskox movement data}
\description{
  Example data based on the movement of a muskox tracked in northeast Greenland.
}
\usage{
muskox
}
\format{
A data frame with 6825 GPS-based observations and 5 columns:
\itemize{
  \item \code{date}: date
  \item \code{tday}: time of day
  \item \code{x}: UTM easting coordinate
  \item \code{y}: UTM northing coordinate
  \item \code{step}: hourly step length calculated based on the coordinates of time t and t+1}
}
\references{
The data set is based on a subset of the data described and analysed in:

Beumer, L.T., Pohle, J., Schmidt, N.M, Chimienti, M., Desforges, J.-P., Hansen, L.H., Langrock, R., Pedersen, S.H., Stelvig, M. and van Beest, F.M. (2020). An application of upscaled optimal foraging theory using hidden Markov modelling: year-round behavioural variation in a large arctic herbivore". \emph{Movement Ecology}, \strong{8}, https://doi.org/10.1186/s40462-020-00213-x.
}
\source{
Schmidt, N.M. (2020). Data for \emph{"An application of upscaled optimal foraging theory using hidden Markov modelling: year-round behavioural variation in a large arctic herbivore"} [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3768080.
}