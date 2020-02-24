\name{periodogram}
\alias{periodogram}

\title{Periodogram Power Spectral Density}
\usage{
periodogram(x, method = "builtin")
}
\arguments{
  \item{x}{vector or matrix containing the time series data
         (one time series per column)}
  \item{method}{a string that specifies which method should be used to
     compute the spectral density: "builtin" employs the function
     \code{\link{spectrum}} with the options 
     taper=0, plot=FALSE, fast=FALSE, detrend=FALSE, and demean=TRUE; 
     "clone" employs directly the Fourier transform function \code{\link{fft}}
     (with sames results as "builtin"); and "smooth" uses the
     function \code{\link{spectrum}} with options as above plus span=3.
     
     }
}
\description{
  \code{periodogram} is a wrapper function for \code{\link{spectrum}}
  with some special options set. It
   returns the power spectral density, i.e. the
  squared modulus of the Fourier coefficient divided by the length
  of the series, for multiple time series as well as the corresponding 
  Fourier frequencies. The frequencies range between 
   0 and the Nyquist critical frequency fc = \code{\link{frequency}}(x)/2. 
  
  \code{periodogram}  is used by the functions
  \code{\link{avgp}} and \code{\link{fisher.g.test}}.
  For general periodogram functions
  please refer to \code{\link{spectrum}}.
}

\value{

  A list object with the following components:
  \item{spec}{A vector or matrix with the estimated power spectral densities
             (one column per time series).}
  \item{freq}{A vector with frequencies f ranging from 0 to fc 
              (if the sampling rate \code{\link{frequency}}(x)) equals 1 then fc = 0.5).
              Angular frequencies may be obtained by multiplication with 2*pi
	      (i.e. omega = 2*pi*f).}
}


\author{
  Konstantinos Fokianos and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{spectrum}}, \code{\link{avgp}}, \code{\link{fisher.g.test}}.}

\examples{
# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)

# how many genes and how many samples?
dim(caulobacter)

# periodograms of the first 10 genes
periodogram(caulobacter[,1:10])
}
\keyword{ts}
