# mean
mean_g <- function(x, gmean = "sqrtsin") {
    switch(
        gmean,
        "linear" = 0.2 + 0.4 * x,
        "log" = log(x),
        "sin" = sin(2 * pi * x),
        "lsin" =  0.5 + 2 * x + sin(2 * pi * x - 0.5),
        "sqrtsin" = sqrt(x * (1 - x)) * sin(2 * pi * (1 + 2 ^ (-7 / 5)) / (x + 2 ^ (-7 / 5)))
    )
}
# error
err_sig <- function(x, xsig = "quad", derr = "norm") {
    v_sig <- switch(
        xsig,
        "const" = 0.2,
        "linear" = 0.2 + 0.2 * x,
        "quad" = 0.5 + 0.5 * (x - 1) ^ 2
    )
    switch(
        derr,
        "norm" = rnorm(x, 0),
        "chisq3" = rchisq(x, 3),
        "t3" = rt(x, 3),
        "t1" = rt(x, 1)
        #"laplace" = rla
    ) * v_sig
}
# error quantile
err_qt <- function(x, tau = 0.5, xsig = "quad", derr = "norm") {
    v_sig <- switch(
        xsig,
        "const" = 0.2,
        "linear" = 0.2 + 0.2 * x,
        "quad" = 0.5 + 0.5 * (x - 1) ^ 2
    )
    switch(
        derr,
        "norm" = qnorm(tau, 0),
        "chisq3" = qchisq(tau, 3),
        "t3" = qt(tau, 3),
        "t1" = qt(tau, 1)
        #"laplace" = rla
    ) * v_sig
}
obs_f <- function(x, gmean = "sqrtsin", xsig = "quad", derr = "norm") {
    mean_g(x, gmean) + err_sig(x, xsig, derr)
}
truth_f <- function(x, tau, gmean = "sqrtsin", xsig = "quad", derr = "norm") {
    mean_g(x, gmean) + err_qt(x, tau, xsig, derr)
}

# x <- sort(x)
# plot(x, obs_f(x, derr = "chisq3"))
# for (tau in seq(0.001, 0.999, length.out = 10)) {
#     lines(comp[[3]][[1]]$x, truth_f(comp[[3]][[1]]$x, tau, gmean = "linear", xsig = "linear", derr = "t3"))
# }
