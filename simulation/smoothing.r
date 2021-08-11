library(devtools)
library(parallel)
devtools::load_all()
RNGkind("L'Ecuyer-CMRG")
set.seed(999983)

taus <- seq(0.05, 0.95, length.out = 10)
k1 <- 11; ktau <- 2; m1 <- 4; mtau <- 3
v_k <- c(ktau, k1); v_m <- c(mtau, m1); v_dim <- v_k + v_m
n <- 400
ncores <- 16
m <- 20
v_sptau <- 10 ^ (-5 : 0); v_spx <- 10^ (-5:4)


# linear
f1 <- function(x) {0.2 + 0.4 * x}
# logarithm
f2 <- function(x) {log(x)}
# sinusoidal
f3 <- function(x) {sin(2 * pi * x)}
# square root sinusoidal
f4 <- function(x) {
    sqrt(x * (1 - x)) * sin(2 * pi * (1 + 2 ^ (-7 / 5)) / (x + 2 ^ (-7 / 5)))
}
obs_f <- function(sim, x) {
    if (sim == 1) {
        f1(x) + rnorm(length(x), mean = rep(0, length(x)), sd = 0.4 * (0.5 + x))
    } else if (sim == 2) {
        f2(x) + rnorm(length(x), mean = rep(0, length(x)), sd = 0.4 * (0.5 + x))
    } else if (sim == 3) {
        f3(x) + rnorm(length(x), mean = rep(0, length(x)), sd = 0.4 * (0.5 + x))
    } else if (sim == 4) {
        f4(x) + rnorm(length(x), mean = rep(0, length(x)), sd = 0.4 * (0.5 + x))
    }
}
truth_f <- function(sim, x, tau) {
    if (sim == 1) {
        f1(x) + qnorm(tau, sd = 0.4 * (0.5 + x))
    } else if (sim == 2) {
        f2(x) + qnorm(tau, sd = 0.4 * (0.5 + x))
    } else if (sim == 3) {
        f3(x) + qnorm(tau, sd = 0.4 * (0.5 + x))
    } else if (sim == 4) {
        f4(x) + qnorm(tau, sd = 0.4 * (0.5 + x))
    }
}

smooth_code2 <- function(sp, m = 20, sim = 2) {
    # L2 loss
    L2loss <- matrix(rep(0, 1 * length(taus)), nrow = 1)
    colnames(L2loss) <- taus
    rownames(L2loss) <- c("quantile sheets") #, "smooth_quan", "coco_quan")
    for (j in seq(m)) {
        #train
        x <- runif(n); ord <- order(x)
        x <- x[ord]
        y <- obs_f(sim, x)
        # test
        test_x <- seq(min(x), max(x), length.out = 1e4)
        #quantile sheets
        qs_fit <- qs_scam(x, y, maxit = 10, ntaus = 10,
                    arg_bs = c("tedmi", "ps"),
                    opt = "efs", dims = v_dim,
                    sp = sp[1:2])
        #smooth quan
#        smooth_fit <- f_smoothdescent(y, x, v_k, v_m,
#                        alpha = 0.1, bet = 0.8,
#                        v_smooth = sp,
#                        h = max((((log(length(y))+1)/length(y))^0.4), 0.05),
#                        s_ntau = 100,
#                        maxit = 500, eps = 1e-5, eta = 1e-3)
#        #coco_quan
#        coco_fit <- f_descent(y, x, v_k, v_m,
#                    alpha = 0.1, bet = 0.8, v_smooth = sp,
#                    maxit = 500, eps = 1e-5)
#
        for (i in seq_along(taus)) {
            tau <- taus[i]
            truth <- truth_f(sim, test_x, tau)
            L2loss[1, i] <- L2loss[1, i] + mean((c(predict(qs_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
#            L2loss[2, i] <- L2loss[2, i] + mean((c(cq_prediction(tau, smooth_fit, newdata = test_x)) - truth)^2)
#            L2loss[3, i] <- L2loss[3, i] + mean((c(cq_prediction(tau, coco_fit, newdata = test_x)) - truth)^2)
       }
    }
    list(loss = L2loss / m,
        model = qs_fit)
}


sim_smooth <- mclapply(v_sptau, outer_loop <- function(sptau) {
    mclapply(v_spx, inner_loop <- function(spx) {
        smooth_code2(sp = c(sptau, spx, spx), m)
    }, mc.cores = ncores)
}, mc.cores = ncores)

save(sim_smooth, file = "sim_smooth.RData")



lapply(sim_smooth, inner <- function(gri) {
        lapply(gri, fn <- function(gr) {
            if (!is.character(gr$loss)) {
                rowMeans(gr$loss)
            } else {
                0
            }
        })
})
