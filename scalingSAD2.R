setwd('~/Dropbox/Research/maxentConcept')

library(meteR)
library(parallel)
library(socorro)
library(plyr)

## source hidden function from meteR to deal with areas
source('~/Dropbox/Research/meteR/R/sar_helper_funs.R')

## wd storing the data
dataWD <- '~/Dropbox/Research/data/stri'

## z-values for SADs across scales, for all plots

nscale <- 5

scaleZ <- mclapply(list.files(dataWD)[3:7], 
                   mc.cores = 6,
                   function(f) {
    x <- read.csv(file.path(dataWD, f), as.is = TRUE)
    ## only look at most recent census
    x <- x[x$year == max(x$year), ]
    
    ## logrithmic intervals for scaling back area
    xscale <- round(max(x$x)) * 2^(-(nscale-1):0)
    yscale <- round(max(x$y)) * 2^(-(nscale-1):0)
    
    ## loop over scales and calculate sad
    out <- sapply(1:length(xscale), function(i) {
        newx <- x[x$x <= xscale[i] & x$y <= yscale[i], ]
        thisSAD <- sad(meteESF(newx$spp, newx$count, newx$dbh^2))
        thisSim <- logLikZ(thisSAD, nrep = 499, return.sim = TRUE)
        
        ## z2 for spatial downscaling
        spatZ <- thisSim$z
        
        ## calculate z2 for each value in the simulated logLiks
        # I0 <- matrix(1, nrow = length(thisSim$sim), ncol = length(thisSim$sim))
        # diag(I0) <- 0
        # m <- as.vector(1/(length(thisSim$sim) - 1) * (I0 %*% thisSim$sim))
        # s <- sqrt(diag(outer(m, thisSim$sim, '-')^2 %*% I0) / (length(thisSim$sim) - 2))
        # trueZ <- ((thisSim$sim - m) / s)^2
        trueZ <- thisSim$sim
        trueZ <- c(mean(trueZ), quantile(trueZ, c(0.025, 0.975)))
        names(trueZ) <- c('mean', 'ciLo', 'ciHi')
        
        ## calculate z2 for random sample of individuals
        randX <- replicate(100, {
            x[sample(1:nrow(x), size = thisSAD$state.var['N0']), ]
        }, simplify = FALSE)
        
        randS <- unique(sapply(randX, function(d) length(unique(d$spp))))
        
        randZ <- lapply(randS, function(s) {
            theseX <- randX[randS == s]
            s <- logLikZ(sad(meteESF(theseX[[1]]$spp, theseX[[1]]$count)), 
                         nrep = 499, return.sim = TRUE)
            if(length(theseX) > 1) {
                zz <- sapply(2:length(theseX), function(i) {
                    obs <- logLik(sad(meteESF(theseX[[i]]$spp, theseX[[i]]$count)))
                    return(((obs - mean(s$sim)) / sd(s$sim))^2)
                })
                return(c(s$z, zz))
            } else {
                return(s$z)
            }
        })
        
        randZ <- unlist(randZ)
        randZ <- c(mean(randZ), quantile(randZ, c(0.025, 0.975)))
        names(randZ) <- c('mean', 'ciLo', 'ciHi')
            
        return(c(spatZ = spatZ, trueZ = trueZ, randZ = randZ))
    })
    
    areas <- xscale * yscale
    out <- data.frame(a = areas, t(out))
    
    return(out)
})

scaleZ <- do.call(rbind, scaleZ)
scaleZ$site <- rep(gsub('.csv', '', list.files(dataWD))[1:2], each = nscale)

write.csv(scaleZ, file = 'scaleZ2.csv', row.names = FALSE)
