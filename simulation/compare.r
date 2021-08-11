devtools::load_all()
library(parallel)
library(quantregGrowth)
library(quantreg)
RNGkind("L'Ecuyer-CMRG")
set.seed(999983)

ncores <- 16
n <- 512; m <- 10; lengthtau <- 1024; lengthtest <- 1e4;
v_taus <- seq(0.001, 0.999, length.out = lengthtau)
k1 <- 11; ktau <- 2; m1 <- 4; mtau <- 3
v_k <- c(ktau, k1); v_m <- c(mtau, m1); v_dim <- v_k + v_m


quantile_regession <- function(x, y, constr, taus = v_taus) {
    fit_quanreg <- list()
    for (i in seq_along(taus)) {
        tau <- taus[i]
        fit_quanreg[[i]] <- rqss(y ~ qss(x, constraint = constr), tau = tau, data = data.frame(x, y))
    }
    return(fit_quanreg)
}

# compare functions(increasing in x & no constraint)
compare_incr <- function(m = 20, taus = v_taus, ...) {
    # L2 loss
    L2loss <- matrix(rep(0, length(taus) * 5), nrow = 5)
    colnames(L2loss) <- taus
    rownames(L2loss) <- c("quantile_sheet(gam)", "quantile sheets(scam)",
                        "cq_reg", "quantreg",  "quantregGrowth")
    # time
    m_time <- matrix(rep(0, 3 * 5), nrow = 5)
    rownames(m_time) <- rownames(L2loss)

    for (j in seq(m)) {
        #train
        x <- runif(n); ord <- order(x)
        x <- x[ord]
        y <- obs_f(x, ...)
        # test
        test_x <- seq(min(x), max(x), length.out = lengthtest)
        #quantile sheets
        m_time[1, ] <- m_time[1, ] + system.time(gam_fit <- qs_gam(x, y,
                maxit = 50, ntaus = 10,
                dims = v_dim, ords = NA,
                sp = NULL, tol = 1e-3))[1:3]
        m_time[2, ] <- m_time[2, ] + system.time(scam_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = 10,
                    arg_bs = c("tedmi", "cr"),
                    opt = "nlm.fd", dims = v_dim))[1:3]
        # cq.reg
        m_time[3, ] <- m_time[3, ] + system.time(cq_fit <- cq.reg(x, y,
                    arg_bs = c("tedmi", "cr"), opt = "bfgs",
                    dims = v_dim, span = 0.05))[1:3]
        # quantreg
        m_time[4, ] <- m_time[4, ] + system.time(quan_fit <- quantile_regession(x, y, "I"))[1:3]
        # quantregGrowth non-crossing+ monotone
        m_time[5, ] <- m_time[5, ] + system.time(fit_wp <-  gcrq(y ~ ps(x, monotone = 1),
                tau = taus, data = data.frame(x, y)))[1:3]
        test_fit <- predict(fit_wp, data.frame(x = test_x))

        for (i in seq_along(taus)) {
            tau <- taus[i]
            truth <- truth_f(test_x, tau, ...)
            L2loss[1, i] <- L2loss[1, i] + mean((c(predict(gam_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[2, i] <- L2loss[2, i] + mean((c(predict(scam_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[3, i] <- L2loss[3, i] + mean((predict(cq_fit, data.frame(x1 = test_x, tau = tau)) - truth)^2)
            L2loss[4, i] <- L2loss[4, i] + mean((predict(quan_fit[[i]], data.frame(x = test_x)) - truth)^2)
            L2loss[5, i] <- L2loss[5, i] + mean((test_fit[, i] - truth)^2)
        }
    }

    list(L2loss = L2loss / m,
        time = m_time / m,
        model = list(gam_fit, scam_fit,
        cq_fit, quan_fit, fit_wp))
}
compare_no <- function(m = 20, taus = v_taus, ...) {
    # L2 loss
    L2loss <- matrix(rep(0, length(taus) * 5), nrow = 5)
    colnames(L2loss) <- taus
    rownames(L2loss) <- c("quantile_sheet(gam)", "quantile sheets(scam)",
                        "cq_reg", "quantreg",  "quantregGrowth")
    # time
    m_time <- matrix(rep(0, 3 * 5), nrow = 5)
    rownames(m_time) <- rownames(L2loss)
    for (j in seq(m)) {
        #data set
        x <- sort(runif(n))
        y <- obs_f(x, ...)
        # test
        test_x <- seq(min(x), max(x), length.out = lengthtest)
        #quantile sheets
        m_time[1, ] <- m_time[1, ] + system.time(gam_fit <- qs_gam(x, y,
                maxit = 50, ntaus = 10,
                dims = v_dim, ords = NA,
                sp = NULL, tol = 1e-3))[1:3]
        m_time[2, ] <- m_time[2, ] + system.time(scam_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = 10,
                    arg_bs = c("tesmi1", "cr"),
                    opt = "nlm.fd", dims = v_dim))[1:3]
        # cq.reg
        m_time[3, ] <- m_time[3, ] + system.time(cq_fit <- cq.reg(x, y,
                    arg_bs = c("tesmi1", "cr"), opt = "bfgs",
                    dims = v_dim, span = 0.05))[1:3]
        # quantreg
        m_time[4, ] <- m_time[4, ] + system.time(
                    quan_fit <- quantile_regession(x, y, "N"))[1:3]
        # quantregGrowth non-crossing+ monotone
        m_time[5, ] <- m_time[5, ] + system.time(fit_wp <-  gcrq(y ~ ps(x),
                tau = taus, data = data.frame(x, y)))[1:3]
        test_fit <- predict(fit_wp, data.frame(x = test_x))

        for (i in seq_along(taus)) {
            tau <- taus[i]
            truth <- truth_f(test_x, tau, ...)
            L2loss[1, i] <- L2loss[1, i] + mean((c(predict(gam_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[2, i] <- L2loss[2, i] + mean((c(predict(scam_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[3, i] <- L2loss[3, i] + mean((c(predict(cq_fit, data.frame(x1 = test_x, tau = tau))) - truth)^2)
            L2loss[4, i] <- L2loss[4, i] + mean((predict(quan_fit[[i]], data.frame(x = test_x)) - truth)^2)
            L2loss[5, i] <- L2loss[5, i] + mean((test_fit[, i] - truth)^2)
        }
    }

    list(L2loss = L2loss / m,
        time = m_time / m,
        model = list(gam_fit, scam_fit,
        cq_fit, quan_fit, fit_wp))
}
# comp <- compare_no(1, gmean = "linear", xsig = "linear", derr = "t3")
# comp <- compare_no(1, gmean = "sqrtsin", xsig = "quad", derr = "chisq3")

v_gmean <- c("linear", "log", "sin", "lsin", "sqrtsin")
v_xsig <- c("const", "linear", "quad")
v_derr <- c("norm", "chisq3", "t3", "t1")

comp <- mclapply(v_gmean, gmean_loop <- function(gmean) {
    mclapply(v_xsig, xsig_loop <- function(xsig) {
        mclapply(v_derr, derr_loop <- function(derr) {
            compare_no(m, gmean = gmean, xsig = xsig, derr = derr)
        }, mc.cores = ncores)
    }, mc.cores = ncores)
}, mc.cores = ncores)


save(comp, file = "sim_compare.RData")