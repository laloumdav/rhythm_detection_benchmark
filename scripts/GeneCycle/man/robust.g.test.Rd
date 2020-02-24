\name{robust.g.test}
\alias{robust.g.test}
\alias{robust.spectrum}


\title{Robust g Test for Multiple (Genetic) Time Series}
\description{
  \code{robust.g.test} calculates the p-value(s) for a robust
   nonparametric version of Fisher's g-test (1929). Details
  of this approach are described in Ahdesmaki et al. (2005), along with
  an extensive discussion of its application to gene expression data.
 From GeneCycle 1.1.0 on the robust regression based method published 
 in Ahdesmaki et al. (2007) is also implemented (using Tukey's biweight
 based M-estimation/regression.)
  
  \code{robust.spectrum} computes a robust rank-based estimate
  of the periodogram/correlogram - see  Ahdesmaki et al. (2005)
  for details. Alternatively it can also be used (since GeneCycle 1.1.0)
  for evaluating the robust regression based spectral estimates,
  suitable for processing non-uniformly sampled data (unknown 
  periodicity time: return spectral estimates, known periodicity 
  time: return p-values).
}
\usage{
robust.g.test(y, index, perm = FALSE, x, noOfPermutations = 300, 
algorithm=c("rank", "regression"), t)
robust.spectrum(x, algorithm = c("rank", "regression"), t, 
periodicity.time = FALSE, noOfPermutations = 300)
}
\arguments{
  \item{y}{the matrix consisting of the spectral estimates
            as column vectors}
  \item{index}{an index to the spectral estimates (RANK BASED 
  APPROACH ONLY; for specifying a periodicity time 
 in the regression approach, see the parameter 
 periodicity.time) that is to be used in the 
 testing for periodicity. If \code{index} is 
 missing for the rank based approach, the maximum 
 component of the spectral estimate is used in 
 testing (regardless of the frequency of this 
 maximum)}
  \item{periodicity.time}{time (same units as in vector \code{t}) of 
    period where periodicity will be detected (ROBUST 
    REGRESSION BASED APPROACH ONLY) that is to be used 
    in the search for periodicity. If periodicity.time 
    is not given for the regression based approach, 
    the whole spectrum is evaluated (more time consuming)
    and the maximum periodogram ordinate will be investigated}
  \item{perm}{if \code{perm} is \code{FALSE}, 
              a simulated distribution for the g-statistic is used 
    (applies to the rank based approach only). 
    If per \code{perm} is \code{TRUE}, permutation 
    tests are used to find the distribution of the 
    g-statistic for each time series separately. 
    With the regression based approach (Ahdesmaki 
    et al. 2007) permutation tests will always be used}
  \item{x}{a matrix consisting of the time series as column
           vectors. In \code{robust.g.test} only needed if
           permutation tests are used}
  \item{noOfPermutations}{number of permutations that 
 are used for each time series (default = 300)}
  \item{algorithm}{\code{rank} corresponds to the rank based approach 
(Ahdesmaki et al. 2005) and 
 \code{regression} for the regression based approach (Ahdesmaki et al. 
2007), which is more suitable for time
 series with non-uniform sampling (default = rank)}
  \item{t}{sampling time vector (only for the regression based 
approach)}
}


\details{
  Application of \code{robust.g.test} can be very computer intensive, 
especially
  the production of the distribution of the test statistics may take a 
lot
  of time. Therefore, this distribution (dependening on the length of
  the time series) is stored in an external file to avoid recomputation
  (see example below). When applying permutation tests no external file 
is
  used but the computation time will always be high.
  
  For the general idea behind the Fisher's g test also see  
  \code{\link{fisher.g.test}} which implements an analytic approach for 
g-testing.
  This is faster but not robust and also assumes Gaussian noise.

  Note that when using the regression based approach there will regularly 
  be warnings about the non-convergence of the regression (iteration limit 
  default at 20 cycles in \link{rlm}).
  
}
\value{
 \code{robust.g.test} returns a list of p-values.
 \code{robust.spectrum} returns a matrix where the column vectors 
correspond
 to the spectra corresponding to each time series. As an exception, if 
the robust regression 
 based approach (Ahdesmaki et al. 2007) is used with a known periodicity 
time, the function 
 robust.spectrum returns p-values (computation will take a lot of time 
depending on how many
 permutations are used per time series and time series length).
}
\author{
  Miika Ahdesmaki (\email{miika.ahdesmaki@gmail.com}).
}
\references{
  Fisher, R.A. (1929). Tests of significance in harmonic analysis. 
\emph{Proc. Roy. Soc. A}, \bold{125}, 54--59.


  Ahdesmaki, M., Lahdesmaki, H., Pearson, R., Huttunen, H., and 
Yli-Harja O.
  (2005).  \emph{BMC Bioinformatics} \bold{6}:117. \url{http://www.biomedcentral.com/1471-2105/6/117}


  Ahdesmaki, M., Lahdesmaki, H., Gracey, A., Shmulevich, I., and 
Yli-Harja O.
  (2007).  \emph{BMC Bioinformatics} \bold{8}:233. \url{http://www.biomedcentral.com/1471-2105/8/233}

}
\seealso{\code{\link[fdrtool]{fdrtool}}, \code{\link{fisher.g.test}}.}

\examples{
\dontrun{

# load GeneCycle library
library("GeneCycle")

# load data set
data(caulobacter)

# how many samples and and how many genes?
dim(caulobacter)


# robust, rank-based spectral estimator applied to first 5 genes
spe5 = robust.spectrum(caulobacter[,1:5])

# g statistics can be computed from the spectrum (internal use mostly 
# but can be checked here)
## g.statistic(spe5)

# robust p-values, use Monte Carlo simulation (not permutation tests) 
# to estimate the null hypothesis distribution
pval = robust.g.test(spe5)  # generates a file with the name "g_pop_length_11.txt"
pval = robust.g.test(spe5)  # second call: much faster..

pval

# robust p-values, now look at index 4 (index can be anything from 1 
# (DC-level) to N (length of the time series and highest frequency))
pval = robust.g.test(spe5, 4)  # generates a file
pval = robust.g.test(spe5, 4)  # second call: much faster..


pval

# delete the external files 
unlink("g_pop_length_11.txt")
unlink("g_pop_length_11indexed.txt")

#
# Next let us see how the robust regression based approach can be 
# applied (Ahdesmaki et al. 2007)
# First: Unknown frequencies
t=c(0,15,30,45,60,75,90,105,120,135,150)
y = robust.spectrum(x=caulobacter[,1:5],algorithm="regression", t=t)
pvals = robust.g.test(y = y, perm=TRUE, x=caulobacter[,1:5], 
noOfPermutations = 50, algorithm = "regression", t=t)
 
pvals

#
# The following example illustrates how to use the regression based 
# method if we have prior knowledge about the frequency/period time
# of periodicity
t = 0:9 # time indices
t = t + runif(10)-0.5 # make time indices non-uniform
A = 0.5 * matrix(rnorm(50),10,5)       # create random time series (no outliers)
A[,5]=A[,5]+matrix(sin(0.5*pi*t),10,1) # superimpose a sinusoidal
periodicity.time=4                     # where to look for periodicity
# note that now the function robust.spectrum returns the p-values (in 
# all other cases it will return spectral estimates):
pvals=robust.spectrum(x=A,algorithm="regression", 
t=t,periodicity.time=periodicity.time, noOfPermutations=50)
pvals  # 5th p-value is smallish, as expected

}
}
\keyword{htest}

