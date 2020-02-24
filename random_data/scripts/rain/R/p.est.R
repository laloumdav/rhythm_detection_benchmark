# Different Methods to estimate the Propability-distribution using
# R.GMP for 'Multiple Precision Arithmetic' Author: thaben
library("gmp")

#' Folding arrays
#' 
#' computes multiplication of ploynomes from elements like (1-i^x)^n 
#' whith x e N and n e (+1,-1)
#' using gmp for exact numerics
#' @param arr array of polynomial coefficients
#' @param f number describing both indices. x <- abs(f) and n <- sign(f)
#' @param len length of arr
#' @return folding of the new element on the remaining array
#' 
#' @author thaben
#' @noRd
foldarr <- function(arr, f, len) {
    
    # decompose the information saved in f
    st <- abs(f)
    sig <- sign(f)
    
    # computation the to effect the series are skipped (should NEVER happen)
    if (st >= len) 
        return(arr)
    
    # seperate the cases sig = +1 and sig = -1
    if (sig == 1) {
        # construct an array where the initial array is shifted by st
        arrs <- matrix(c(as.bigz(rep(0, st)), arr[, 1:(len - st)]), 
            ncol = len)
        # and substract it
        arr2 <- arr - arrs
        return(arr2)
    } else {
        # this case is done by multiplication of a matrix (see Harding 
        # for details)
        seq <- c(1, rep(0, st - 1))
        mmat <- do.call("rbind", lapply(1:len, function(j) {
            c(rep(0, j - 1), rep(seq, len = len - j + 1))
        }))
        return(arr %*% mmat)
    }
}


# new version with l as a series of extrema

#' Exact Distributions 
#' 
#' Using the harding Procedure to generate exact Distributions for Umbrella 
#' and Rain Tests
#' @param mList List of group sizes 
#' @param l list of inflection points
#' @param cycl boolean: first point is a copy of the last one and has to be
#' treated specially. Is only evaluated if only one single inflection > 1 
#' is present. Parameter is meaningless without any inflection
#' @return list of pValues and the probability density
#' 
#' @author thaben
#' @noRd 
harding <- function(mList, l = NULL, cycl = FALSE) {
    
    len <- length(mList)
    N <- sum(mList)
    nvars <- 0
    maxval <- 0
    if (is.null(l)) {
        genfunc <- 1:N
        for (m in mList) {
            genfunc <- c(genfunc, -seq_len(m))
        }
        maxval <- (sum(mList)^2 - sum(mList^2))/2
    } else {
        # special treatment of statistic when a cyclic case is tested. In a
        # first step the first element is excluded and all variables are
        # adapted
        
        if (cycl && length(l) == 1 && l[1] != 1) {
            l <- l - 1
            special <- mList[1]
            mList <- mList[-1]
            len <- length(mList)
            N <- sum(mList)
        } else {
            cycl <- FALSE
            special <- 0
        }
        
        # very special case sometimes leeding to problems
        if(all(l == 1)){
            special <- 0
        }
        
        # including first and last positions
        lList <- c(1, l[l > 1 & l < len], len)
        
        # calculating the size of whole slopes
        Ns <- sapply(seq_len(length(lList) - 1), function(i) {
            sum(mList[lList[i]:lList[i + 1]])
        })
        Ns <- Ns[Ns > 0]
        
        # generating function for whole slopes
        genfunc <- do.call("c", lapply(Ns, seq_len))
        
        # generating function for all elements perhaps alowing non
        # complete time series
        for (m in mList) {
            genfunc <- c(genfunc, -seq_len(m))
        }
        
        # generating function for inflection points (have to be counted
        # twice)
        genfunc <- c(genfunc, do.call("c", lapply(l[l > 1 & l < len], 
            function(i) { 
                -seq_len(mList[i])
            })))
        
        # calculate maximum Value of test
        maxval <- sum(sapply(1:(length(lList) - 1), function(i) {
            (sum(mList[lList[i]:lList[i + 1]])^2 - 
                sum(mList[lList[i]:lList[i + 1]]^2)) / 2
        }))
        
        # second step for the cyclic case treatment genetrating function
        # gets the man-whitney test for the special element whith the
        # first slope exept the next inflection maxval is extendet by the
        # maximum possible counts for this test
        if (special != 0 && cycl) {
            addpos <- 1:(lList[2] - 1)
            addpos <- addpos[addpos > 0]
            remain <- sum(mList[addpos])
            genfunc <- c(genfunc, c(-seq_len(special)), c(seq_len(special)
                + remain))
            maxval <- maxval + special * remain
        }
        
    }
    
    # return for the 'intestable'
    if (maxval == 0) 
        return(list(pval = c(1), den = c(0)))
    
    ## as Distribution is symetric only the first half has to be
    ## calculated
    arr <- as.bigz(matrix(c(1, rep(0, ceiling(maxval / 2) - 1)), nrow = 1))
    len <- length(arr)
    
    # use the generating function
    for (gen in genfunc) {
        arr <- foldarr(arr, gen, len)
    }
    
    # combine to a full distribution, calculate density
    if (maxval %% 2 == 1) {
        arr <- cbind(arr, arr[, (ncol(arr) - 1):1])
    }
    if (maxval %% 2 == 0) {
        arr <- cbind(arr, arr[, ncol(arr):1])
    }
    den <- arr/sum(arr)
    
    # return
    return(list(pval = as.double((rev(cumsum(den)))), den = as.double(den)))
} 
