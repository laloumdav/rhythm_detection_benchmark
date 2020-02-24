# Main Function of the rain algorithm Author: thaben

#' String based code for distinct model settings
#' 
#' @param sequence sequence of measurement groups
#' @param extremas index of extremas
#' @param cycl boolean: statistic treats cyclic arrangement of samples
#' differently
#' @return code string
#' @noRd 
#' 
#' @author thaben
generateCode <- function(sequence, extremas, cycl = FALSE) {
    
    # different behavoiur for cyclic tests
    if (cycl & length(extremas) == 1 & extremas[1] != 1) {
        l <- extremas
        
        # seperate the inflections
        words <- sapply(list(sequence[1], sequence[l]), function(i) {
            paste(sort(i), ".", sep = "", collapse = "")
        })
        # part between start and inflection
        if (l > 2) {
            words <- c(words, paste(sort(sequence[2:(l - 1)]), collapse = "."))
        }
        # part between inflection and end
        if (l < length(sequence) - 1) {
            words <- c(words, paste(sort(sequence[(l + 1):(length(sequence) - 
                1)]), collapse = "."))
        }
        # collapse and return
        return(paste(sort(words), collapse = " "))
    }
    # write valid list of extremas
    ex <- c(1, extremas[extremas > 1 & extremas < length(sequence)], 
        length(sequence))
    # generate contentStrings for each slope
    words <- sapply(2:length(ex), function(i) {
        paste(sort(sequence[ex[i - 1]:ex[i]]), collapse = ".")
    })
    # sort slopes and combine them
    return(paste(sort(words), collapse = " "))
}

#' Search for rhytrhmicities by using umbrella alternatives
#' 
#' use rain() to access this function
#' @inheritParams rain
#' @param tSer numeric array; one row per time point, one colum per object of 
#' evaluation 
#' @param periods numeric array list of periods to lookup. 
#' @param peaks borders of peak shape
#' @param type numeric index of detection method
#' @return data frame containing best mtching phase, period, 
#' paek-shape and pValue
#' 
#' @noRd
#' @author thaben
#' @import gmp
#' @import multtest
umbrellaCirc <- function(tSer, nr.series = 1, periods = c(nrow(tSer)), 
    peaks = c(0.35, 0.65), type = 1, adjp.method, measure.sequence, 
    raw.pVal,
    verbose = TRUE) {
    
    # preparing lists
    if (is.null(measure.sequence)) {
        tp <- nrow(tSer) / nr.series
        measure.sequence <- rep(nr.series, tp)
    } else {
        tp <- length(measure.sequence)
    }
    peakpos <- c()
    relpeak <- c()
    periodlist <- c()
    test.cases <- 0
    hardings <- list()
    decoder <- data.frame(code = "null", harding = 0)
    compmat <- list()
    distris <- c()
    
    # there are several ways to calculate the statistic and so
    # type = 1 is 'longitudinal' in the main function type = 3 is 'independent'
    # type = 2 is a deprecated in between solution. It is completely working, 
    # but the performance ist not good enugh to be public visible
    
    # in this part the comparison matrices, statistics and so on are prepared. 
    # Afterwards the time series are evaluated with these statistics
    
    # 'longitudinal'
    if (type == 1) {
        # extended min function to avoid warnings if running min on empty
        # list
        emin <- function(x) {
            if (length(x) == 0) 
                return(Inf)
            return(min(x))
        }
        
        #prepare some additional lists
        ups <- c()
        extremas <- list()
        
        #go through all periods to measure
        for (period in periods) {
            maxtp <- ceiling(tp / period) * period
            
            # specify how the assymmetric scale peaks[] translates in this 
            # period setting
            perpeaks <- (floor(peaks[1] * period)):(ceiling(peaks[2] * 
                period))
            perpeaks <- perpeaks[perpeaks > 0 & perpeaks < period - 
                1]
            
            # the list of peakpositions for all phases and peak shapes under 
            # this period length is prepared
            peakl <- matrix(nrow = period * length(perpeaks), ncol = 1)
            peakl[, 1] <- rep(1:period, each = length(perpeaks))
            peakl <- lapply(peakl[, 1], seq, to = maxtp, by = period)
            
            # specify all peaks and troughs for all phases and peak shapes
            perpeakpos <- rep(1:period, each = length(perpeaks))
            peakpos <- c(peakpos, perpeakpos)
            perrelpeak <- rep(perpeaks, period)
            relpeak <- c(relpeak, perrelpeak)
            
            trough <- lapply(1:length(peakl), function(x) {
                peakl[[x]] + perrelpeak[[x]]
            })
            
            peakl <- lapply(peakl, function(x) {
                y = x %% maxtp
                return(y[y > 1 & y < tp])
            })
            
            trough <- lapply(trough, function(x) {
                y = x %% maxtp
                return(y[y > 1 & y < tp])
            })
            
            # check whether a peak or a trough comes first to evaluate if the
            # modell first rises or falls
            perups <- apply(rbind(sapply(peakl, emin), sapply(trough, 
                emin)), 2, function(v) {
                v[1] <= v[2]
            })
            ups <- c(perups, ups)
            
            # write it to the global exrtremas list
            periodlist <- c(periodlist, rep(period, length(perups)))
            extremas <- c(extremas, lapply(1:length(peakl), function(x) {
                sort(unique(c(peakl[[x]], trough[[x]])))
            }))
        }
        
        # filter out results without any peak or trough to avoid
        # meaningless results
        filter <- sapply(extremas, function(val) {
            length(val) != 0
        })
        ups <- ups[filter]
        periodlist <- periodlist[filter]
        extremas <- extremas[which(filter)]
        peakpos <- peakpos[filter]
        relpeak <- relpeak[filter]
        
        test.cases <- length(ups)
        
        trials <- measure.sequence
        
        coordtable <- rbind(c(1, cumsum(measure.sequence)[-tp] + 1), 
            cumsum(measure.sequence))
        
        # preparing the compare matrices for all different settings
        compmat <- lapply(1:test.cases, function(i) {
            manWilcoxMatrix(coords = lapply(c(1:tp), function(x) {
                if (coordtable[1, x] > coordtable[2, x]) 
                    return(c())
                coordtable[1, x]:coordtable[2, x]
            }), l = extremas[[i]], up = ups[i])
        })
        
        if (verbose) 
            message("\r\ncalculating distributions (", test.cases, 
                "):", appendLF = FALSE)
        distris <- numeric(test.cases)
        
        # calculate the statistics for all settings. To save computation time
        # the function generate code is used to check for similar previos 
        # calculated statistics
        for (i in 1:test.cases) {
            code <- generateCode(trials, extremas[[i]])
            
            ## cat(code,'\r\n')
            if (code %in% decoder$code) {
                if (verbose) message(".", appendLF = FALSE)
                pos <- which(decoder$code == code)
                distris[i] <- (decoder$harding[pos])
            } else {
                if (verbose) message("*", appendLF = FALSE)
                hardingpos <- length(hardings) + 1
                decoder <- rbind(decoder, data.frame(code = code, 
                    harding = hardingpos))
                hardings[[hardingpos]] <- harding(trials, extremas[[i]])$pval
                distris[i] <- (hardingpos)
            }
        }
        if (verbose) 
            message("\r\nDone")
    }
    
    # deprecated version 
    if (type == 2) {
        # new lists
        extremas <- list()
        sequences <- list()
        runind <- 1
        for (period in periods) {
            
            # prepare the coordtables and check which time points are repeats 
            # of others for a given period
            use.pers <- max(floor(tp/period), 1)
            perpeaks <- (floor(peaks[1] * period)):(ceiling(peaks[2] * 
                period))
            perpeaks <- perpeaks[perpeaks > 0 & perpeaks < period]
            coordtable <- rbind(c(1, cumsum(measure.sequence)[-tp] + 
                1), cumsum(measure.sequence))
            coords <- lapply(1:(period * use.pers), function(i) {
                if (i > tp) 
                    return(NULL)
                l <- do.call(c, lapply(seq(i, tp, period * use.pers), 
                    function(x) {
                        if (coordtable[1, x] > coordtable[2, x]) 
                            return(c())
                        return(coordtable[1, x]:coordtable[2, x])
                    }))
                return(l)
            })
            
            # generate the comparisonmatrices for a setting 
            group.size <- sapply(coords, length)
            for (phase in 1:period) {
                maxs <- seq(phase, period * use.pers, period)
                perextremas <- lapply(perpeaks, function(peak) {
                    mins <- (maxs + peak - 1)%%(period * use.pers) + 1
                    return(c(mins, maxs))
                })
                percompmat <- lapply(perextremas, function(extremlist) {
                    split <- extremlist[1]
                    extremlist <- sort((extremlist + 1 - split)%%(period * 
                        use.pers))
                    rotmat <- c(split:(period * use.pers), 1:(split - 1)
                        )[1:(period * use.pers)]
                    manWilcoxMatrix(coords = coords[rotmat], l = extremlist, 
                        up = TRUE)
                })
                
                peakpos <- c(peakpos, rep(phase, length(perpeaks)))
                relpeak <- c(relpeak, perpeaks)
                compmat <- c(compmat, percompmat)
                extremas <- c(extremas, lapply(perextremas, 
                    function(extremlist) {
                        split <- extremlist[1]
                        extremlist <- sort((extremlist + 1 - split) %% 
                            (period * use.pers))
                        
                    }
                ))
                
                periodlist <- c(periodlist, rep(period, length(percompmat)))
                sequences <- c(sequences, lapply(perextremas, 
                    function(extremlist) {
                        split <- extremlist[1]
                        rotmat <- c(split:(period * use.pers), 
                                    1:(split - 1))[1:(period * use.pers)]
                        return(group.size[rotmat])
                    }
                ))
            }
        }
        
        test.cases <- length(compmat)
        
        # calculate the statistics reuse similar statistics
        if (verbose) 
            message("\r\ncalculating distributions (", test.cases, 
                "):", appendLF = FALSE)
        distris <- numeric(test.cases)
        for (i in 1:test.cases) {
            code <- generateCode(sequences[[i]], extremas[[i]])
            if (code %in% decoder$code) {
                if (verbose) 
                    message(".", appendLF = FALSE)
                pos <- which(decoder$code == code)
                distris[i] <- (decoder$harding[pos])
            } else {
                if (verbose) 
                    message("*", appendLF = FALSE)
                hardingpos <- length(hardings) + 1
                decoder <- rbind(decoder, data.frame(code = code, 
                    harding = hardingpos))
                hardings[[hardingpos]] <- harding(sequences[[i]], 
                    extremas[[i]])$pval
                distris[i] <- (hardingpos)
            }
        }
        
        if (verbose) 
            message("\r\nDone")
        
    }
    
    # rain3 combine all datapoints in one circular Umbrella
    if (type == 3) {
        
        extremas <- list()
        sequences <- list()
        
        for (period in periods) {
            
            # for each period define the position of the peak
            perpeaks <- c((ceiling((1 - peaks[1]) * period)):(floor((1 - 
                peaks[2]) * period)))
            perpeaks <- perpeaks[perpeaks > 0 & perpeaks < period]
            
            # coordtable holds the beginning and ending index for each
            # timepoint
            coordtable <- rbind(c(1, cumsum(measure.sequence)[-tp] + 
                1), cumsum(measure.sequence))
            
            # for each phase setup the frames
            for (phase in 1:period) {
                
                # arrange the groups so that peak is at the first point and 
                # each group contains all measurement
                coordlists <- lapply(phase:(phase + period), function(i) {
                    i <- ((i - 1) %% period) + 1
                    if (i > tp) 
                        return(NULL)
                    l <- do.call(c, lapply(seq(i, tp, period), function(x) {
                        if (coordtable[1, x] > coordtable[2, x]) {
                            return(c())
                        }
                        return(coordtable[1, x]:coordtable[2, x])
                    }
                    ))
                    return(l)
                })
                
                # generate comparison matrices
                percompmat <- lapply(perpeaks, function(peak) {
                    manWilcoxMatrix(coords = coordlists, l = peak + 
                        1, up = FALSE)
                })
                
                # write all to global setting tables peakposition
                peakpos <- c(peakpos, rep(phase, length(perpeaks)))
                # Number of points from peak to trough
                relpeak <- c(relpeak, perpeaks)
                # comparison matrices
                compmat <- c(compmat, percompmat)
                # period lengths
                periodlist <- c(periodlist, rep(period, length(percompmat)))
                # number of points in the slopes
                sequences <- c(sequences, lapply(perpeaks, function(x) {
                    sapply(coordlists, length)
                }))
            }
        }
        test.cases <- length(compmat)
        
        if (verbose) 
            message("\r\ncalculating distributions (", test.cases, 
                "):", appendLF = FALSE)
        
        distris <- numeric(test.cases)
        for (i in 1:test.cases) {
            # generate codestring to avoid repeated calculation of
            # distributions
            code <- generateCode(sequences[[i]], c(1 + relpeak[i]), 
                cycl = TRUE)
            
            # if code ist still present enter the right harding distribution
            # into the decoder list else make new distribution
            if (code %in% decoder$code) {
                if (verbose)
                    message(".", appendLF = FALSE)
                pos <- which(decoder$code == code)
                distris[i] <- (decoder$harding[pos])
            } else {
                if (verbose)
                    message("*", appendLF = FALSE)
                hardingpos <- length(hardings) + 1
                decoder <- rbind(decoder, data.frame(code = code, 
                    harding = hardingpos))
                hardings[[hardingpos]] <- harding(sequences[[i]], 
                    c(1 + relpeak[i]), cycl = TRUE)$pval
                distris[i] <- (hardingpos)
            }
        }
        if (verbose) 
            message("\r\nDone")
    }
    # end of generating Statistical distributions
    
    # evaluating the Dataset
    if (verbose)
        message("\r\nEvaluating Datasets: ")
    
    # preparing resulttables
    pVal <- numeric(ncol(tSer))
    phase <- numeric(ncol(tSer))
    peak <- numeric(ncol(tSer))
    period <- numeric(ncol(tSer))
    
    # progressBar
    if (verbose) 
        pb <- txtProgressBar(min = 0, max = ncol(tSer), style = 3, 
            width = 30)
    
    # adjust.p.method correction
    col.adjp.method <- ifelse(adjp.method == "TSBH", "TSBH_0.05", 
        adjp.method)
    # for each single time series
    for (ind in 1:ncol(tSer)) {
        # reset progressbar
        if (verbose) 
            setTxtProgressBar(pb, ind)
        
        # pick column
        list <- tSer[, ind]
        
        # create comparison matrix of real Data
        testcomp <- sign(sapply(list, function(x) {
            list - x
        }))
        testcomp[is.na(testcomp)] <- 0
        
        # compare this compMatrix, with the ones calculated for different
        # settings
        scores <- sapply(seq_len(test.cases), function(i) {
            resmat <- compmat[[i]] * testcomp
            return(sum(resmat[resmat > 0])/2)
        })
        
        # for these scores find theaccording p-Values
        pvals <- sapply(seq_len(test.cases), function(i) {
            p <- hardings[[distris[i]]][scores[i]]
            if (length(p) == 0) 
                return(1)
            return(p)
        })
        # bad-value protection
        pvals[scores == 0] <- 1
        
        # if no 'real' pValues could be detected give direct output to
        # avoid problems
        if (min(pvals) == max(pvals)) {
            pVal[ind] <- 1
            best <- 1
        } else {
            # pvalue adjustments
            adjusted <- suppressWarnings(mt.rawp2adjp(pvals, 
                proc = adjp.method))
            best <- adjusted$index[1]
            
            if (raw.pVal){
              pVal[ind] <- adjusted$adjp[1, "rawp"]   ### Add option to get raw p-values
            } else {
              pVal[ind] <- adjusted$adjp[1, col.adjp.method]
            }
            if (is.na(pVal[ind])) 
              pVal[ind] <- 1
        }
        
        # write measurement parameters to output tabels
        phase[ind] <- peakpos[best]
        peak[ind] <- relpeak[best]
        period[ind] <- periodlist[best]
    }
    if (verbose) 
        message("\r\nDone\r\n")
    # return all
    return(data.frame(pVal = pVal, phase = phase, peak = peak, 
        period = period))
} 
