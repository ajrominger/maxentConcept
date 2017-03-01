library(meteR)
library(parallel)

dataWD <- '~/Dropbox/Research/data/stri'

pdf('~/Dropbox/Research/fig_striSAD.pdf', width = 4, height = 4)

striZ2 <- sapply(list.files(dataWD), FUN = function(f) {
    x <- read.csv(file.path(dataWD, f))
    x <- x[x$year == max(x$year), ]
    esf <- meteESF(x$spp, x$count, x$dbh^2)
    thisSAD <- sad(esf)
    thisZ <- logLikZ(thisSAD, nrep = 499)$z
    
    plot(thisSAD, ptype = 'rad', log = 'y', add.legend = FALSE)
    legend('topright', legend = c(gsub('.csv', '', f), round(thisZ)))
    
    return(thisZ)
})

dev.off()