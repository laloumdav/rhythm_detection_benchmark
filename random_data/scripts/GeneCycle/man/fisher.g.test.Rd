\name{fisher.g.test}
\alias{fisher.g.test}
\title{Fisher's Exact g Test for Multiple (Genetic) Time Series}
\description{
  \code{fisher.g.test} calculates the p-value(s) according to Fisher's
  exact g test for one or more time series. This test is useful to detect hidden
  periodicities of unknown frequency in a data set.  For an application to
  microarray data see Wichert, Fokianos, and Strimmer (2004). 
}
\usage{
fisher.g.test(x, \dots)
}
\arguments{
  \item{x}{vector or matrix with  time series data (one time series per column).}
  \item{\dots}{arguments passed to \code{\link{periodogram}}}
}
\details{
  Fisher (1929) devised an exact procedure to test the null hypothesis of Gaussian
  white noise against the alternative of an added deterministic periodic component
  of unspecified frequency.  The basic idea behind the test is to reject the 
  null hypothesis if the periodogram contains a value significantly larger
  than the average value (cf. Brockwell and Davis, 1991).  
  This test is useful in the context of microarray genetic time series
  analysis as a gene selection method - see Wichert, Fokianos and Strimmer (2004)
  for more details.   Note that in the special case of a constant time series
  the  p-value returned by \code{fisher.g.test} is 
  exactly 1 (i.e. the null hypothesis is not rejected).
}
\value{
  A vector of p-values (one for each time series).  Multiple testing
  may then be done using the the false discover rate approach
  (function \code{\link[fdrtool]{fdrtool}}).
}
\author{
  Konstantinos Fokianos and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{
  Fisher, R.A. (1929). Tests of significance in harmonic analysis. \emph{Proc. Roy. Soc. A}, \bold{125}, 54--59.

  Brockwell, P.J., and Davis, R.A.  (1991). Time Series: Theory and Methods (2nd ed).  Springer Verlag.
  (the g-test is discussed in section 10.2).

  Wichert, S., Fokianos, K., and Strimmer, K. (2004).  Identifying periodically expressed transcripts
   in microarray time series data. \emph{Bioinformatics} \bold{20}:5-20.
}
\seealso{\code{\link[fdrtool]{fdrtool}}.}

\examples{

# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)

# how many samples and and how many genes?
dim(caulobacter)

# p-values from Fisher's g test
pval.caulobacter <- fisher.g.test(caulobacter)
pval.caulobacter

# compute Fdr and fdr values
fdr.out <- fdrtool(pval.caulobacter, statistic="pvalue")

# how many significant?
sum(fdr.out$qval < 0.05) # tail area-based Fdr
sum(fdr.out$lfdr < 0.2)  # density-based local fdr

}
\keyword{htest}
