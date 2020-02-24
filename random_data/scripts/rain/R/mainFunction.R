# main Functions to access Umbrella and hewitt test R hythmic A
# nalysis I ncoorperating N on-parametric methods Author: Thaben


#' Detection of rhythmic behavior in time series.
#' 
#' rain detects rhythms in time-series using non parametric methods.
#' It uses an extension of the rank test for Umbrella Alternatives 
#' (Mack & Wolfe, 1981), based on on the Jonckheere-Terpstra test, which tests
#' whether sets of groups have a trend or not.  The Umbrella method extends 
#' this to independent rising and falling sets.
#' @param x numeric array, containing the data. One row per time point, 
#' one column per sample. If more than one replicate is done, see 
#' nr.series for formatting.
#' @param deltat numeric: sampling interval. 
#' @param period numeric: Period to search for. 
#' The given period is mapped to the best matching number of measurements. 
#' A set of periods can be defined by period and period.delta.
#' @param period.delta numeric: width of period interval. 
#' A interval of different period-length to evaluate is defined by 
#' period $-$ period.delta and period $+$ period.delta. 
#' In this interval all possible numbers of time points according to the 
#' deltat are tested.
#' @param peak.border vector c(min,max): defines the different form of the 
#' peak. min and max have to be >0 and <1. The concrete interpretation depends
#' on the chosen method (see Details). (default = c(0.3,0.7))
#' @param nr.series numeric: Number of replicates of the whole time series. 
#' If using nr.series all series have to have the same length. 
#' These multiple time series contained in x must be organized timepoint by 
#' timepoint using the following format 
#' [r1t1, r2t1, r1t2, r2t2, ..., r1tn, r2tn] 
#' where ritj is the i'th repeat of the j'th time-point. 
#' @param measure.sequence numeric array: Numbers of replicates for each time
#' point. 
#' By using 'measure.sequence', irregular time courses may be evaluated.
#' A value of 0 is possible and handeled correctly. The array determines how 
#' many values are present for each time piont. 
#' The values are ordered in the same format as specified above. 
#' measure.sequence overwrites nr.series if both are set.
#' @param method string ('independent', 'longitudinal'): identify the method
#' to use (see Details).
#' @param na.rm boolean: calculate individual statistics for time series 
#' containign NAs. The time series of a sample containing NAs is treated
#' as if the time points with NA are not measured. Using this option increases 
#' calculation time.
#' @param adjp.method string (see \code{\link[multtest]{mt.rawp2adjp}}): 
#' select the method wich is used for the multiple testing correction for the 
#' different phases and periods and shapes tested
#' @param verbose status output
#' @return An array containing p-values and the description of the best
#' matching model. 
#' Each row apply to a sample in x.
#'    \item{pVal }{The p-Values}
#'    \item{phase }{The phase of the peak (see Details)}
#'    \item{peak.shape }{The shape of the curve depending on the method used
#'    (see Details)}
#'    \item{period }{The period length, same unit as deltat}
#' @details 
#' The method tests whether a the time course consists of alternating rising 
#' and falling slopes, repeated with a distinct period. The partitions of the 
#' rising part with respect to the whole period are given by 
#' peak.border = c(min, max). The value peak.shape specifies this partition in 
#' the best matching model. The phase is defined as the time point with the 
#' peak. There are two versions of umbrella:
#' \describe{
#' \item{independent}{Multiple periods are interpreted as repeats of one 
#' period.}
#' \item{logitudinal}{The whole time series remains unaffected. Partial slopes
#' in the beginning and end of the time series are evaluated as shorter 
#' slopes. This method implicitly rejects underlying trends.This should be 
#' used only with longitudinal samples, hat may contain strong trends}
#' }
#' @author Paul F. Thaben
#' @references 
#' Mack, G. A., & Wolfe, D. A. (1981). K-Sample Rank Tests for Umbrella
#' Alternatives. 
#' \emph{Journal of the American Statistical Association},
#' \bold{76(373)}, 175--181.
#' 
#' @examples 
#' # create a dataset with different noise levels
#' noise.levels <- c(1, 0.5, 0.2, 0.1, 0.05, 0.02)
#' period <- 15
#' testset <- apply(matrix(noise.levels, nrow = 1), 2, function(noise){
#'    timecourse = 1 + 0.4 * cos((1:30) / period * 2 * pi) + 
#'    rnorm(30, 0, noise)
#' })
#' 
#' 
#' results <- rain(testset, period=15, deltat=1, method='independent')
#' 
#' plot(-log(results$pVal) ~ noise.levels)
#' 
#' \dontrun{
#' # testing a biological dataset
#' data(menetRNASeqMouseLiver) 
#' menet.ossc <- rain(t( menetRNASeqMouseLiver ), deltat = 4, period = 24, 
#'    nr.series = 2, peak.border = c(0.3, 0.7), verbose=TRUE)
#' require('lattice')
#' 
#' best <- order(results$pVal)[1:10]
#' 
#' xyplot(as.matrix(menetRNASeqMouseLiver
#'    [best, (0:5 * 2 + rep(c(1, 2), each = 6))]) 
#'  ~rep(0:11 * 4 + 2, each = 10) |rownames(menetRNASeqMouseLiver)[best], 
#'  scales = list(y = list(relation = 'free')),
#'  layout = c(2, 5), type = 'b', pch = 16, xlab = 'time', 
#'  ylab = 'expression value', cex.lab = 1)
#' 
#' }
#' @keywords Statistics|nonparametric Statistics|ts
#' @export
rain <- function(x, deltat, period, period.delta = 0, peak.border = c(0.3, 
    0.7), nr.series = 1, measure.sequence = NULL, method = "independent", 
    na.rm = FALSE, adjp.method = "ABH", verbose = getOption("verbose"), 
    raw.pVal = FALSE) {
    
    # catch one dimensional vectors as interpreting them as one time
    # course
    if (is.null(dim(x))) 
        x <- matrix(x, ncol = 1)
    
    
    # create measure sequence if nessecary
    if (is.null(measure.sequence)) {
        measure.sequence = rep(nr.series, floor(nrow(x) / nr.series))
    }
    
    # check for correct data annotation
    if (!is.null(measure.sequence) && sum(measure.sequence) != nrow(x)) {
        stop("invalid annotation of time points")
    }
    
    # different behaviour for na.rm
    if (na.rm) {
        # look for na
        valmasks <- apply(x, 2, is.na)
        
        # make a list of the measurement sequence, for each sample column
        # (NA's are removed)
        seqLists <- apply(valmasks, 2, reduceByNA, mseq = measure.sequence)
        
        # get a List of unique measurement sequences occuring in Dataset
        codes <- apply(seqLists, 2, paste, collapse = "")
        uniqcode <- unique(codes)
        
        # empty resultframe as spaceholder
        result <- data.frame(pVal = rep(NA, ncol(x)), phase = rep(NA, 
            ncol(x)), peak = rep(NA, ncol(x)), period = rep(NA, ncol(x)))
        
        per <- floor((period - period.delta) / deltat):ceiling((period + 
            period.delta) / deltat)
        
        # some output
        if (verbose) {
            message("\r\ndeploying ", length(uniqcode), " statistics", 
                appendLF = TRUE)
            pb <- txtProgressBar(min = 0, max = ncol(x), style = 3, 
                width = 30)
        }
        counter = 0
        
        # go through all different variants of NA-Distributions
        for (code in uniqcode) {
            
            pos <- which(codes == code)
            counter = counter + length(pos)
            if (verbose) 
                setTxtProgressBar(pb, counter)
            
            subSet <- x[, pos]
            
            # if subset contains only 1 series, the matrix has to be forced
            if (is.null(dim(subSet))) 
                subSet <- matrix(subSet, ncol = 1)
            
            #if no measurements are left repeat 'empty' results
            if (all(is.na(subSet))) {
                len = ncol(subSet)
                result[pos, ] <- data.frame(pVal = rep(1, len), 
                    phase = rep(0, len), peak = rep(0, len), 
                    period = rep(period/deltat, len))
                next
            }
            
            # shrink the test set by removing the NA values
            subSet <- apply(subSet, 2, function(x) {
                x[!is.na(x)]
            })
            if (is.null(dim(subSet))) 
                subSet <- matrix(subSet, nrow = 1)
            
            # use the measure.sequence without the NA-measurements
            sub.measure.sequence <- seqLists[, pos[1]]
            
            # run rain
            if (method %in% c("umb1", "longitudinal")) {
                part.result <- umbrellaCirc(subSet, nr.series = 1, 
                    per, peak.border, type = 1, adjp.method = adjp.method, 
                    raw.pVal = raw.pVal,
                    sub.measure.sequence, verbose = FALSE)
            }
            if (method == "umb2") {
                part.result <- umbrellaCirc(subSet, nr.series = 1, 
                    per, peak.border, type = 2, adjp.method = adjp.method, 
                    raw.pVal = raw.pVal,
                    sub.measure.sequence, verbose = FALSE)
            }
            if (method %in% c("umb3", "independent")) {
                part.result <- umbrellaCirc(subSet, nr.series = 1, 
                    per, peak.border, type = 3, adjp.method = adjp.method, 
                    raw.pVal = raw.pVal,
                    sub.measure.sequence, verbose = FALSE)
            }
            
            # save the results in the result table
            result[pos, ] <- part.result
        }
    } else {
        # run rain in the choosen format
        if (method %in% c("umb1", "longitudinal")) {
            per <- floor((period - period.delta) / deltat):ceiling((period + 
                period.delta) / deltat)
            result <- umbrellaCirc(x, nr.series = 1, per, peak.border, 
                type = 1, adjp.method = adjp.method, measure.sequence, 
                raw.pVal = raw.pVal,
                verbose = verbose)
        }
        if (method == "umb2") {
            per <- floor((period - period.delta) / deltat):ceiling((period + 
                period.delta) / deltat)
            result <- umbrellaCirc(x, nr.series = 1, per, peak.border, 
                type = 2, adjp.method = adjp.method, measure.sequence, 
                raw.pVal = raw.pVal,
                verbose = verbose)
        }
        if (method %in% c("umb3", "independent")) {
            per <- floor((period - period.delta) / deltat):ceiling((period + 
                period.delta) / deltat)
            result <- umbrellaCirc(x, nr.series = 1, per, peak.border, 
                type = 3, adjp.method = adjp.method, measure.sequence, 
                raw.pVal = raw.pVal,
                verbose = verbose)
        }
    }
    
    colnames(result)[3] <- "peak.shape"
    result["phase"] <- result["phase"] * deltat
    result["period"] <- result["period"] * deltat
    result["peak.shape"] <- result["peak.shape"] * deltat
    if (!is.null(colnames(x))) 
        rownames(result) <- make.unique(colnames(x))
    return(result)
}


reduceByNA <- function(valmask, mseq) {
    valmask <- ifelse(valmask, 0, 1)
    apply(rbind(cumsum(mseq) - mseq + 1, cumsum(mseq)), 2, function(coord) {
        if (coord[1] <= coord[2] & coord[2] > 0) {
            return(sum(valmask[coord[1]:coord[2]]))
        } else {
            return(0)
        }
    })
} 
