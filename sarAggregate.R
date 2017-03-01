setwd('~/Dropbox/Research')

library(iNEXT)
library(meteR)
library(parallel)

data(BCI, package = 'vegan')
data(BCI.env, package = 'vegan')

bciTot <- colSums(BCI)

pdf('fig_bciSARaggregate.pdf', width = 6, height = 3)
out <- lapply(1:nrow(BCI), 
       # mc.cores = 6, 
       FUN = function(i) {
           bciOne <- BCI[i, ][BCI[i, ] > 0]
           bciOneMETE <- meteESF(as.character(1:length(bciOne)), bciOne)
           bciOneZ2 <- logLikZ(sad(meteESF(as.character(1:length(bciOne)), bciOne)), 
                               nrep = 499)$z
           bciOneSAR <- upscaleSAR(bciOneMETE, 1, 2^6)
           bciOneHill0 <- as.numeric(estimateD(bciOne, datatype = 'abundance',
                                               level = mean(rowSums(BCI)) * 
                                                   max(bciOneSAR$A))[1, 5:7])
           
           
           par(mfrow = c(1, 2), mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
           
           plot(sad(bciOneMETE), ptype = 'rad', log = 'y', add.legend = FALSE)
           legend('topright', legend = round(bciOneZ2, 3), bty = 'n')
           
           plot(bciOneSAR, log = 'xy', ylim = range(ncol(BCI), bciOneHill0, bciOneSAR$S))
           points(50, ncol(BCI), col = 'blue')
           arrows(x0 = 50, y0 = bciOneHill0[2], y1 = bciOneHill0[3], code = 3, 
                  length = 0.05, angle = 90,
                  col = 'red')
           points(50, bciOneHill0[1], col = 'red', bg = 'white', pch = 21)
           
           return(c(sadZ = bciOneZ2, SAR = max(bciOneSAR$S), hill = bciOneHill0))
})

dev.off()
