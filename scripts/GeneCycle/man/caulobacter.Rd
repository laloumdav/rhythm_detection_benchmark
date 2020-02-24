\name{caulobacter}
\alias{caulobacter}

\title{Microarray Time Series Data for 1444 Caulobacter  Crescentus 
       Genes}

\description{
  This data set describes the temporal expression of 1444 genes
  (open reading frames) in the cell cycle of the bacterium
  \emph{Caulobacter crescentus}.
}

\usage{
data(caulobacter)
}

\format{
  \code{caulobacter} is a \code{\link[longitudinal]{longitudinal}} object
  containing the data from the Laub et al. (2000) experiment.
  Essentially, this is a matrix with with 1444 columns (=genes)
  and 11 rows (=time points) 
}

\source{
  This data is described in Laub et al. (2000).
}


\references{
  Laub, M.T., McAdams, H.H., Feldblyum, Fraser, C.M., and Shapiro, L.
  (2000) Global analysis of the genetic network controlling
  a bacterial cell cycle. \emph{Science}, \bold{290}, 2144--1248.
}

\examples{
# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)
is.longitudinal(caulobacter)

# how many samples and how many genes?
dim(caulobacter)
summary(caulobacter)
get.time.repeats(caulobacter)

# plot first nine time series
plot(caulobacter, 1:9)
}
\keyword{datasets}

