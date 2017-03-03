

## function to (properly) randomly sample an SAD
sample.sad <- function(x, n) {
    s <- sample(1:length(x), size = n, replace = TRUE, prob = x)
    out <- rowSums(outer(1:length(x), s, '=='))
    return(out)
}

P <- 0.25
plot.samp <- function(x, ...) {
    plot(x^P, ylim = c(0, max(x^P)), yaxt = 'n', ...)
    alab <- c(0, 10^seq(0, log(max(x)), by = 1))
    axis(2, at = alab^P, labels = alab)
}

points.samp <- function(x, ...) {
    points(x^P, ...)
}

poly.samp <- function(x, ...) {
    polygon(y = c(x, 0, 0)^P, x = c(1:length(x), length(x), 1), ...)
}


## common is common everywhere and rare is rare everywhere
x <- meteDist2Rank(sad(meteESF(S0 = 100, N0 = 10000)))
xsub1 <- sample.sad(x, round(sum(x)*0.5))
xsub2 <- sample.sad(x, round(sum(x)*0.5))
xcomb <- xsub1 + xsub2

plot.samp(x, type = 'n')
poly.samp(x, col = 'gray', border = NA)
points.samp(xsub1, type = 'l', col = hsv(0.05, 1, 0.8))
points.samp(xsub2, type = 'l', col = hsv(0.08, 0.8, 0.9))
points.samp(xcomb, type = 'l', col = hsv(0.65, 0.7, 1))


## common is common everywhere and rare turn over
x <- meteDist2Rank(sad(meteESF(S0 = 100, N0 = 10000)))
xsub1 <- sample.sad(x, round(sum(x)*0.5))
xsub2 <- sample.sad(x, round(sum(x)*0.5))

rareCutOff <- 0.005
rareCutOff <- ifelse(max(x)*rareCutOff < 0.5, 1, round(max(x)*rareCutOff))
xcomb <- ifelse(x > rareCutOff, xsub1 + xsub2, 0)
xcomb <- c(xcomb, xsub1[x <= rareCutOff], xsub2[x <= rareCutOff])
xcomb <- sort(xcomb[xcomb > 0], TRUE)

plot.samp(x, type = 'n', xlim = c(1, length(xcomb)))
poly.samp(x, col = 'gray', border = NA)
points.samp(xcomb)
points.samp(meteDist2Rank(sad(meteESF(1:length(xcomb), xcomb))), col = 'red', type = 'l')


## common is rare (=absent) somewhere and rare is rare everywhere
x <- meteDist2Rank(sad(meteESF(S0 = 100, N0 = 10000)))
xsub1 <- sample.sad(x, round(sum(x)*0.5))
xsub2 <- sample.sad(x, round(sum(x)*0.5))

rareCutOff <- 0.05
rareCutOff <- ifelse(max(x)*rareCutOff < 0.5, 1, round(max(x)*rareCutOff))
xcomb <- ifelse(x <= rareCutOff, xsub1 + xsub2, 0)
xcomb <- c(xcomb, xsub1[x > rareCutOff], xsub2[x > rareCutOff])
xcomb <- sort(xcomb[xcomb > 0], TRUE)

# plot.samp(x, type = 'n', xlim = c(1, length(xcomb)))
# poly.samp(x, col = 'gray', border = NA)
# points.samp(xcomb)
# points.samp(meteDist2Rank(sad(meteESF(1:length(xcomb), xcomb))), col = 'red', type = 'l')
plot(x, log = 'y', xlim = range(c(1, length(xcomb), length(x))), 
     type = 'l', col = 'blue', lwd = 3)
points(xcomb, col = gray(0, alpha = 0.5), pch = 16)
points(meteDist2Rank(sad(meteESF(1:length(xcomb), xcomb))), col = 'red', type = 'l')
