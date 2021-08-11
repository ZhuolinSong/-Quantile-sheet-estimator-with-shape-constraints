require(mgcv)
# first an un-penalized example - fit E(y)=a+bx subject to a>0
set.seed(0)

# Penalized example: monotonic penalized regression spline .....
n <- 100
# Generate data from a monotonic truth.
x <- runif(100) * 4 - 1
x <- sort(x)
f <- exp(4 * x) / (1 + exp(4 * x))
y <- f + rnorm(100) * 0.1
plot(x, y)
dat <- data.frame(x = x, y = y)
# Show regular spline fit (and save fitted object)
f.ug <- gam(y ~ s(x, k = 10, bs = "cr"))
lines(x, fitted(f.ug))
# Create Design matrix, constraints etc. for monotonic spline....
sm <- smoothCon(s(x, k = 10, bs = "cr"), dat, knots = NULL)[[1]]
F <- mono.con(sm$xp)
# get constraints
G <- list(X = sm$X, C = matrix(0, 0, 0), sp = f.ug$sp, p = sm$xp, y = y, w = y * 0 + 1)
G$Ain <- F$A
G$bin <- F$b
G$S <- sm$S
G$off <- 0

p <- pcls(G)
# fit spline (using s.p. from unconstrained fit)

fv <- Predict.matrix(sm, data.frame(x = x)) %*% p
lines(x, fv, col = 2)
