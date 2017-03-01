library(socorro)
library(meteR)
library(MASS)

randBisect <- function(n, b) {
	nn <- vector('list', b)
	
	nn[[1]] <- oneBisect(n)
	for(i in 2:b) {
		nn[[i]] <- unlist(lapply(nn[[i-1]], oneBisect))
	}
	
	return(nn)
}

oneBisect <- function(n) {
	if(n > 0) {
		n1 <- sample(n - 1, size = 1)
		n2 <- n - n1
		return(c(n1, n2))
	} else {
		return(rep(0, 2))
	}
	
}

n <- 500
b <- 8

foo <- randBisect(n, b)

par(mfrow = c(1, 2))

pred <- ssad(meteSSF(n0 = n, A = 2^0, A0 = 2^b))
plot(meteDist2Rank(pred), sort(foo[[b]], TRUE), xlab = 'McGlinn', ylab = 'METE')
abline(0, 1)


nb <- fitdistr(foo[[b]], 'negative binomial')
nb$estimate

plot(simpECDF(foo[[b]]), log = 'y', main = 'Almost neg binom, but k != 1')
curve(pnbinom(x, size = nb$estimate[1], mu = nb$estimate[2]), add = TRUE, col = 'red')
