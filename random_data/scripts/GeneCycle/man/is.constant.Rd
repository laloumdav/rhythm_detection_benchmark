\name{is.constant}
\alias{is.constant}
\title{Simple Check for Constant Time Series}
\usage{
is.constant(x)
}
\arguments{
  \item{x}{vector or matrix with time series data (one time series per column)}
}
\description{
  \code{is.constant} is a utility function that
  checks whether a time series is constant.
}

\value{
  A vector with a boolean statement (\code{TRUE} or \code{FALSE}) for each time series.
}

\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}

\examples{
# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)

# any constant genes?
sum(is.constant(caulobacter))

# but here:
series.1 <- rep(1, 10)
series.2 <- seq(1, 10)
is.constant( cbind(series.1, series.2) )

}
\keyword{ts}
