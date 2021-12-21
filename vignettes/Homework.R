## -----------------------------------------------------------------------------
sigma = c(0.1, 1, 10, 50)
m = 100
X1 <- sqrt(-2*sigma[1]*log(1-runif(m)))
X2 <- sqrt(-2*sigma[2]*log(1-runif(m)))
X3 <- sqrt(-2*sigma[3]*log(1-runif(m)))
X4 <- sqrt(-2*sigma[4]*log(1-runif(m)))
par(mfrow=c(1,2))
hist(X1, main = "sigma = 0.1")
hist(X2, main = "sigma = 1")
par(mfrow=c(1,2))
hist(X3, main = "sigma = 10")
hist(X4, main = "sigma = 50")

