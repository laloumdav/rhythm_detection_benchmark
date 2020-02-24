#' fuction to calculate the man wilcox sum of two sets
#' 
#' compares all elements of x with all elements of y counts up the result score
#' for each comparison, where elem(x) < elem(y)
#' @param x set of numeric 
#' @param y set of numeric
#' @return the score
#' 
#' @author thaben
#' @noRd
lcomp <- function(x, y) {
    
    return(sum(sapply(x, function(z) {
        sum(ifelse(y > z, 1, 0))
    })))

}

#' calculate the jonkheere terpstra statistic for a series of sets of numbers
#' 
#' @param list a list of vectors of numeric
#' @return the resulting jonkheere terpstra statistic
#' 
#' @author thaben
#' @noRd
manWilcox <- function(list) {
    num.sets <- length(list)
    
    # compare each set with the sets at later position in the list
    # sums up a series of man wilcox comparisons 
    u <- sum(sapply(c(1:(num.sets - 1)), function(i) {
        sum(sapply((i + 1):num.sets, function(x) lcomp(list[[i]], 
            list[[x]])))
    }))

    return(u)
}

#' Calculate Umbrella statistic from a series of numeric vectors
#' 
#' @param x list of numeric vectors
#' @param l set of inflections
#' @param up logical whether the series between the start and the first 
#' inflection is rising or falling
#' @return the statistic
#' 
#' @author thaben
#' @noRd
calcUmbrellaM <- function(x, l, up = TRUE) {
    len <- length(x)
    direcs = rep(c(-1, 1), len = length(l) + 2)
    
    lext = c(1, l[l < len & l > 1], len)
    
    if (up) 
        direcs = direcs[-1]
    Al <- sum(sapply(1:(length(lext) - 1), function(i) {
        vals = x[lext[i]:lext[i + 1]]
        if (direcs[i] == -1) vals = rev(vals)
        return(manWilcox(vals))
    }))
    return(Al)
}


#' Calculation of a matrix, allowing fast man wilcox tests
#'
#' returns a comparisonMatrix of should be results 1 := row > col, 0 :=
#' Not compared or equal, -1 := row < col
#' @param coords list of vectors of numbers depicting the coordinates in the 
#' input vector of samples. Each element of the list represents a set of 
#' grouped samples
#' @param l numeric vector: series of inflection 
#' @param up boolean: if the first part of the series (1..l[1]) is expected to
#' rise
#' @return returns a matrix give the expexted > or < relations in a 
#' sample matrix
#'
#' @author thaben
#' @noRd
manWilcoxMatrix <- function(coords, l, up = TRUE) {
    # count the number of sample sets
    len = length(coords)
    
    # estimate the number of time points
    tp = max(do.call(c, coords))
    
    # add the initial and last time to the series of inflections and sort
    ls = l[l > 1 & l < len]
    ls = sort(c(1, ls, len))
    
    # prepare resulting compare matrix
    comparematrix = matrix(rep(0, tp^2), nrow = tp)
    
    # two markers indicating if for the current part the data should rise or 
    # fall
    rowbiggercol = 1
    colbiggerrow = -1
    
    if (up) {
        rowbiggercol = -1
        colbiggerrow = 1
    }
    
    # go through the coords list by subsets of pure rising or falling series
    for (i in seq_len(length(ls) - 1)) {
        for (x in ls[i]:(ls[i + 1] - 1)) {
            for (y in (x + 1):(ls[i + 1])) {
                comparematrix[coords[[x]], coords[[y]]] = 
                    comparematrix[coords[[x]], coords[[y]]] + rowbiggercol
                comparematrix[coords[[y]], coords[[x]]] = 
                    comparematrix[coords[[y]], coords[[x]]] + colbiggerrow
            }
        }
        rowbiggercol = -rowbiggercol
        colbiggerrow = -colbiggerrow
    }
    comparematrix = sign(comparematrix)
    return(comparematrix)
} 
