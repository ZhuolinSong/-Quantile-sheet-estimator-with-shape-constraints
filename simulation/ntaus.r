library(devtools)
library(parallel)
devtools::load_all()
RNGkind("L'Ecuyer-CMRG")
set.seed(999983)

ncores <- 8
n <- 512; m <- 10; lengthtau <- 1024; lengthtest <- 1e4;
v_taus <- seq(0.001, 0.999, length.out = lengthtau)
k1 <- 11; ktau <- 2; m1 <- 4; mtau <- 3
v_k <- c(ktau, k1); v_m <- c(mtau, m1); v_dim <- v_k + v_m


ntaus_code <- function(ntaus, m = 20, taus = v_taus, ...)  {
    # L2 loss
    L2loss <- matrix(rep(0, length(taus) * 2), nrow = 2)
    colnames(L2loss) <- taus
    rownames(L2loss) <- c("quantile_sheet(gam)", "quantile sheets(scam)")
    # time
    m_time <- matrix(rep(0, 3 * 2), nrow = 2)
    rownames(m_time) <- rownames(L2loss)
    for (j in seq(m)) {
        #data set
        x <- sort(runif(n))
        y <- obs_f(x, ...)
        # test
        test_x <- seq(min(x), max(x), length.out = lengthtest)
        #quantile sheets
        #quantile sheets
        m_time[1, ] <- m_time[1, ] + system.time(gam_fit <- qs_gam(x, y,
                maxit = 50, ntaus = ntaus,
                dims = v_dim, ords = NA,
                sp = NULL, tol = 1e-3))[1:3]
        m_time[2, ] <- m_time[2, ] + system.time(scam_fit <- qs_scam(x, y,
                    maxit = 20, ntaus = ntaus,
                    arg_bs = c("tesmi1", "cr"),
                    opt = "nlm.fd", dims = v_dim))[1:3]

        for (i in seq_along(taus)) {
            tau <- taus[i]
            truth <- truth_f(test_x, tau, ...)
            L2loss[1, i] <- L2loss[1, i] + mean((c(predict(gam_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
            L2loss[2, i] <- L2loss[2, i] + mean((c(predict(scam_fit, data.frame(tau = tau, x1 = test_x))) - truth)^2)
       }
    }
    list(L2loss = L2loss / m,
        time = m_time / m,
        model = list(gam_fit, scam_fit))
}

v_ntau <- 10^ (1:3)
v_gmean <- c("linear", "log", "sin", "lsin", "sqrtsin")
v_xsig <- c("const", "linear", "quad")
v_derr <- c("norm", "chisq3", "t3", "t1")

sim_ntaus <- mclapply(v_gmean, gmean_loop <- function(gmean) {
    mclapply(v_xsig, xsig_loop <- function(xsig) {
        mclapply(v_derr, derr_loop <- function(derr) {
            mclapply(v_ntau, ntau_loop <- function(ntaus) {
                ntaus_code(ntaus, m, gmean = gmean, xsig = xsig, derr = derr)
            }, mc.cores = 1)
        }, mc.cores = 1)
    }, mc.cores = ncores)
}, mc.cores = ncores)

save(sim_ntaus, file = "sim_ntaus.RData")