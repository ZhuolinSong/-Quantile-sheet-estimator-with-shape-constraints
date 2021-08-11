#generate data
x <- runif(400); y <- 4*x + 5;
noise <- rnorm(400, 0, 3*x + 0.5)
y <- y + noise

library(mgcv)
fitgam <- mgcv::gam(y~s(x))
#what scam::scam supposed to do
try(mgcv::gam(y~s(x), start = coef(fitgam)))

require(scam)
fitscam <- scam::scam(y ~ s(x))
# error occurs
try(scam::scam(y ~ s(x,bs="mpi"), start = coef(fitscam)))
# take weird starting values
scam(y ~ s(x,bs="mpi"), start = 1:400)