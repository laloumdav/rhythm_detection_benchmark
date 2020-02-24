\name{dominant.freqs}
\alias{dominant.freqs}
\title{Dominant Frequencies in Multiple (Genetic) Time Series}
\usage{
dominant.freqs(x, m=1, \dots)
}
\arguments{
  \item{x}{multivariate (genetic) time series (each \bold{column} of this
          matrix corresponds to a separate variable/time series), or a vector with a single
	   time series}
  \item{m}{number of dominant frequences}
  \item{\dots}{arguments passed to \code{\link{periodogram}}}
}
\description{
  \code{dominant.freqs} returns the m dominant frequencies (highest peaks)
  in each of the periodogram computed for the individual time series.
}

\value{
  A matrix (or vector, if only 1 time series is considered) with the
  dominant frequencies. In a matrix, each column
  corresponds to one time series.
}

\author{
  Konstantinos Fokianos and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{periodogram}}, \code{\link{spectrum}}.}

\examples{
# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)

# how many samples and how many genes?
dim(caulobacter)

# first three dominant frequencies for each gene
dominant.freqs(caulobacter, 3)

# first four dominant frequencies for gene no. 1000
dominant.freqs(caulobacter[,1000], 4)
}
\keyword{ts}
