library(meteR)
library(socorro)
data(arth)
arth <- arth[-which(arth$mass == min(arth$mass)), ]

x <- meteESF(S0 = length(unique(arth$spp)), N0 = sum(arth$count), 
             E0 = sum(arth$mass^0.75 / min(arth$mass^0.75)))

PofE <- function(e, la1, la2, Z, N0) {
   1/Z * ((1 - exp(-(la1 + la2*e) * N0)) / (1 - exp(-(la1 + la2*e))) + 
              exp(-(la1 + la2*e) * N0) - 1) 
}

dat <- tapply(arth$mass^0.75 / min(arth$mass^0.75), arth$spp, mean)

datCDF <- simpECDF(dat)
thrCDF <- unlist(lapply(datCDF[, 1], function(epsilon) {
    integrate(function(e) PofE(e, x$La[1], x$La[2], x$Z, x$state.var['N0']), 
              lower = 1, upper = epsilon)$value
}))

plot(datCDF)
points(datCDF[, 1], thrCDF, col = hsv(alpha = 0.5))
