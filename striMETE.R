setwd('~/Dropbox/Research/maxentConcept')

library(meteR)
library(parallel)
library(socorro)

dataWD <- '~/Dropbox/Research/data/stri'

pdf('fig_striMETE.pdf', width = 4, height = 4)

striMETE <- sapply(list.files(dataWD), FUN = function(f) {
    x <- read.csv(file.path(dataWD, f))
    x <- x[x$year == max(x$year), ]
    esf <- meteESF(x$spp, x$count, x$dbh^2)
    thisSAD <- sad(esf)
    thisZ <- logLikZ(thisSAD, nrep = 499)$z
    thisSAR <- meteSAR(x$spp, x$count, x = x$x, y = x$y, 
                       Amin = diff(range(x$x)) * diff(range(x$y)) * 2^-8)
    
    par(mfrow = c(1, 3), mar = c(3, 3, 0, 0) + 0.5, mgp = c(2.5, 1, 0))
    
    plot(thisSAD, ptype = 'rad', log = 'y', add.legend = FALSE, yaxt = 'n')
    logAxis(2)
    legend('topright', legend = c(gsub('.csv', '', f), round(thisZ, 3)))
    
    plot(ipd(esf), ptype = 'rad', log = 'y', add.legend = FALSE, yaxt = 'n')
    logAxis(2)
    
    plot(thisSAR, log = 'xy', add.legend = FALSE, xaxt = 'n', yaxt = 'n')
    logAxis(1)
    logAxis(2)
    
    return(list(esf = esf, sar = thisSAR, sadZ2 = thisZ))
})

dev.off()
