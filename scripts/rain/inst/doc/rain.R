### R code from vignette source 'rain.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: rain.Rnw:35-43
###################################################
set.seed(123)
times <- c(1: 24) * 2
sin <- 1 + 0.5 * sin(times / 24 * 2 * pi) + rnorm(24, 0, 0.3)
saw <- rep(13:24 / 18 , 2) + rnorm(24, 0, 0.3)
measure <- cbind(sin, saw)
require('lattice')
xyplot(t(measure)~rep(times, each=2) | c('sin', 'saw'),
    layout = c(1, 2), type = 'o', xlab = 'time', ylab = 'value', cex.lab = 0.6)


###################################################
### code chunk number 3: rain.Rnw:49-55
###################################################
require(rain)
rainresult <- rain(measure, period=24, 
                deltat=2, peak.border=c(0.1,0.9),
                verbose=FALSE
)
rainresult


###################################################
### code chunk number 4: rain.Rnw:99-102
###################################################
data(menetRNASeqMouseLiver)

colnames(menetRNASeqMouseLiver)


###################################################
### code chunk number 5: rain.Rnw:110-121
###################################################
results <- rain(t(menetRNASeqMouseLiver), deltat=4, period=24, nr.series=2,
                peak.border=c(0.3, 0.7), verbose=FALSE) 

best <- order(results$pVal)[1:10]

xyplot(as.matrix(menetRNASeqMouseLiver
            [best, (0:5 * 2 + rep(c(1, 2), each = 6))]) 
        ~rep(0:11 * 4 + 2, each = 10) |rownames(menetRNASeqMouseLiver)[best], 
        scales = list(y = list(relation = 'free')),
        layout = c(2, 5), type = 'b', pch = 16, xlab = 'time', 
        ylab = 'expression value', cex.lab = 1)


###################################################
### code chunk number 6: rain.Rnw:223-229
###################################################
plot(NULL, NULL, xlim = c(0, 1),ylim = c(0, 1), bty = 'n', xaxt = "n", 
    yaxt = "n", xlab = '', ylab = '', mar=c(0,0,1,0))
lines(c(0.2, 0.8), c(1.02, 1.02), lwd = 15, col = 'grey85', lend = 1)
lines(c(0, 0.3, 1), c(0, 1, 0), lwd = 2)
axis(3, c(0, 0.3, 1), labels = c('0', '', '1'), cex = 0.6)
text(0.65, 1.2, "c(0.2, 0.8)", col='grey75', xpd=TRUE)


###################################################
### code chunk number 7: rain.Rnw:244-245
###################################################
rainresult


