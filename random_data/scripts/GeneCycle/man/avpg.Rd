\name{avgp}
\alias{avgp}
\alias{avpg}
\title{Average Periodogram for Multiple (Genetic) Time Series}
\description{
  \code{avgp} calculates and plots the average periodogram as described in
  Wichert, Fokianos and Strimmer (2004).
}
\usage{
avgp(x, title = deparse(substitute(x)), plot = TRUE, angular = FALSE, \dots)
}
\arguments{
  \item{x}{multiple (genetic) time series data. Each \bold{column} of this
          matrix corresponds to a separate variable/time series}
  \item{title}{name of the data set (default is the name of the data object)}
  \item{plot}{plot the average periodogram?}
  \item{angular}{convert frequencies to angular frequencies?}
  \item{\dots}{arguments passed to \code{\link{plot}}
                and to \code{\link{periodogram}}}
}
\details{
  The average periodogram is simply the frequency-wise average of the spectral density (as estimated
  by the Fourier transform) over all times series. 
  To calculate the average periodogram the function \code{\link{periodogram}}
  is used. See Wichert, Fokianos and Strimmer (2004)
  for more details.
}
\value{
  A list object with the following components:
  \item{freq}{A vector with the discrete Fourier frequencies (see \code{\link{periodogram}}).
              If the option angular=TRUE  then the output are angular frequencies (2*pi*f).}
  \item{avg.spec}{A vector with the average power spectral density at each frequency.}
  \item{title}{Name of the data set underlying the average periodogram.}

  The result is returned invisibly if \code{plot} is true.
}
\author{
  Konstantinos Fokianos and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{
  Wichert, S., Fokianos, K., and Strimmer, K. (2004).  Identifying periodically expressed transcripts
   in microarray time series data. \emph{Bioinformatics} \bold{20}:5-20.
}
\seealso{\code{\link{periodogram}}, \code{\link{spectrum}}.}

\examples{
# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)

# how many samples and how many genes?
dim(caulobacter)

# average periodogram
avgp.caulobacter <- avgp(caulobacter, "Caulobacter")
avgp.caulobacter

# just compute and don't plot 
avgp(caulobacter, "Caulobacter", plot=FALSE)
}
\keyword{ts}
