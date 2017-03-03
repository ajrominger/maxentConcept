library(MASS)
library(plyr)
library(socorro)

## source hidden function from meteR to deal with areas
source('~/Dropbox/Research/meteR/R/sar_helper_funs.R')

## wd storing the data
dataWD <- '~/Dropbox/Research/data/stri'

## helper function to make all calculations
ssadLL <- function(f, nrow, ncol) {
    dat <- read.csv(file.path(dataWD, f))
    dat <- dat[dat$year == max(dat$year), ]
    datArea <- .findAreas(dat$spp, dat$count, x = dat$x, y = dat$y, row = nrow, col = ncol)
    dat$cell <- paste(datArea$row, datArea$col, sep = ',')
    
    datMat <- tidy2mat(dat$cell, dat$spp, dat$count)
    datMat <- datMat[, colSums(datMat) > 0]
    
    llObs <- apply(datMat, 2, function(x) 
        unlist(fitdistr(x, 'negative binomial')[c('estimate', 'loglik')]))
    
    llThr <- apply(llObs, 2, function(x) 
        sum(dnbinom(rnbinom(100*nrow(datMat), size = x[1], mu = x[2]), 
                    size = x[1], mu = x[2], log = TRUE)) / 100)
    
    return(data.frame(t(llObs), llThr = llThr, N = colSums(datMat)))
}


## BCI ssad...all negative binomial
bciSSAD <- ssadLL('BCIS.csv', 5, 10)
plot(bciSSAD[, 3:4])
abline(0, 1, col = 'red')

## UCSC...all negative binomial
ucscSSAD <- ssadLL('UCSC.csv', 6, 4)
plot(ucscSSAD[, 3:4])
abline(0, 1, col = 'red')

## Cocoli...all negative binomial
cocoSSAD <- ssadLL('COCO.csv', 6, 2)
plot(cocoSSAD[, 3:4])
abline(0, 1, col = 'red')

## Paso...all negative binomial
pasoSSAD <- ssadLL('PASO.csv', 5, 10)
plot(pasoSSAD[, 3:4])
abline(0, 1, col = 'red')


plot(sort(pasoSSAD[, 1], TRUE), type = 'l', log = 'y', 
     ylim = range(pasoSSAD[, 1], bciSSAD[, 1]), col = 'red', lwd = 2)
points(sort(bciSSAD[, 1], TRUE), type = 'l', col = 'blue', lwd = 2)
points(sort(cocoSSAD[, 1], TRUE), type = 'l', col = 'skyblue', lwd = 2)
points(sort(ucscSSAD[, 1], TRUE), type = 'l', col = 'black', lwd = 2)

plot(sort(pasoSSAD[, 2], TRUE), type = 'l', log = 'y', 
     ylim = range(pasoSSAD[, 2], bciSSAD[, 2]), col = 'red', lwd = 2)
points(sort(bciSSAD[, 2], TRUE), type = 'l', col = 'blue', lwd = 2)
points(sort(cocoSSAD[, 2], TRUE), type = 'l', col = 'skyblue', lwd = 2)
points(sort(ucscSSAD[, 2], TRUE), type = 'l', col = 'black', lwd = 2)

plot(pasoSSAD$estimate.mu, pasoSSAD$N)
