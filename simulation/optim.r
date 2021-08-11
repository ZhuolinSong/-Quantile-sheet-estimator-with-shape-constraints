devtools::load_all()
library(parallel)
library(quantregGrowth)
library(quantreg)
RNGkind("L'Ecuyer-CMRG")
set.seed(999983)

ncores <- 4
n <- 400
m <- 20
lengthout <- 1024
v_taus <- seq(0.001, 0.999, length.out = lengthout)
k1 <- 11; ktau <- 2; m1 <- 4; mtau <- 3
v_k <- c(ktau, k1); v_m <- c(mtau, m1); v_dim <- v_k + v_m


# compare functions(increasing in x & no constraint)
optim_code <- function(m = 20, taus = v_taus, ...) {
    # L2 loss
    L2loss <- matrix(rep(0, length(taus) * 5), nrow = 5)
    colnames(L2loss) <- taus
    rownames(L2loss) <- c("bfgs", "optim",
                        "nlm", "nlm.fd",  "efs")
    # time
    m_time <- matrix(rep(0, 3 * 5), nrow = 5)
    rownames(m_time) <- rownames(L2loss)
    for (j in seq(m)) {
        #train
        x <- runif(n); ord <- order(x)
        x <- x[ord]
        y <- obs_f(x, ...)
        # test
        test_x <- seq(min(x), max(x), length.out = 1e4)
        # bfgs
        m_time[1, ] <- m_time[1, ] + system.time(bfgs_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = 10,
                    arg_bs = c("tesmi1", "cr"),
                    opt = "bfgs", dims = v_dim, tol = 1e-2))[1:3]
        # optim
        m_time[2, ] <- m_time[2, ] + system.time(optim_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = 10,
                    arg_bs = c("tesmi1", "cr"),
                    opt = "optim", dims = v_dim, tol = 1e-2))[1:3]
        # nlm
        m_time[3, ] <- m_time[3, ] + system.time(nlm_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = 10,
                    arg_bs = c("tesmi1", "cr"),
                    opt = "nlm", dims = v_dim, tol = 1e-2))[1:3]
        # nlm.fd
        m_time[4, ] <- m_time[4, ] + system.time(nlmfd_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = 10,
                    arg_bs = c("tesmi1", "cr"),
                    opt = "nlm.fd", dims = v_dim, tol = 1e-2))[1:3]
        # efs
        m_time[5, ] <- m_time[5, ] + system.time(efs_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = 10,
                    arg_bs = c("tesmi1", "cr"),
                    opt = "efs", dims = v_dim, tol = 1e-2))[1:3]
        for (i in seq_along(taus)) {
            tau <- taus[i]
            truth <- truth_f(test_x, tau, ...)
            L2loss[1, i] <- L2loss[1, i] + mean((c(predict(bfgs_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[2, i] <- L2loss[2, i] + mean((c(predict(optim_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[3, i] <- L2loss[3, i] + mean((c(predict(nlm_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[4, i] <- L2loss[4, i] + mean((c(predict(nlmfd_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[5, i] <- L2loss[5, i] + mean((c(predict(efs_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
        }
    }

    list(L2loss = L2loss / m,
        time = m_time / m,
        model = list(bfgs_fit, optim_fit,
        nlm_fit, nlmfd_fit, efs_fit))
}
# comp <- optim_code(1, gmean = "linear", xsig = "linear", derr = "t3")
# comp <- compare_no(1, gmean = "sqrtsin", xsig = "quad", derr = "chisq3")

v_gmean <- c("linear", "log", "sin", "lsin", "sqrtsin")
v_xsig <- c("const", "linear", "quad")
v_derr <- c("norm", "chisq3", "t3", "t1")

opt <- mclapply(v_gmean, gmean_loop <- function(gmean) {
    mclapply(v_xsig, xsig_loop <- function(xsig) {
        mclapply(v_derr, derr_loop <- function(derr) {
            optim_code(m, gmean = gmean, xsig = xsig, derr = derr)
        }, mc.cores = ncores)
    }, mc.cores = ncores)
}, mc.cores = ncores)


save(opt, file = "sim_opt.RData")