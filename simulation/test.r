devtools::load_all()

#examples for tau_star.r
k1 <- 4; ktau <- 4; m1 <- 4; mtau <- 4; y <- rnorm(200)
n1 <- f_n1(1:200, k1, m1); beta <- rnorm(prod(ktau + mtau, k1 + m1))
til_beta <- tilde_beta(beta, k1 + m1)
m_sigma <- f_sigma(ktau + mtau, k1 + m1)
f_qscop(seq(0, 1, 0.1), n1, til_beta, m_sigma, ktau, mtau)
f_qscop(0, n1, til_beta, m_sigma, ktau, mtau) < f_qscop(1, n1, til_beta, m_sigma, ktau, mtau)

tilde_beta(beta, k1 + m1)

f_qscop(.5, t(n1[1, ]), til_beta, m_sigma, ktau, mtau)

f_qzeros(0.5, 2:3, y, n1, til_beta, m_sigma, ktau, mtau)

v_taus <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
sapply(1:200, function(x) {
        f_qscop(v_taus[x], t(n1[x, ]), til_beta, m_sigma, ktau, mtau)
}) - y


#examples for derivative.r
v_taustar <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
m_s <- f_mats(ktau + mtau, k1 + m1, rep(1, 3))

m_g1 <- f_g1(ktau, mtau)
m_g2 <- f_g2(ktau, mtau)
m_sigmatau <- f_sigmatau(ktau + mtau)

m_h1 <- f_h1(n1, m_g1, m_g2, m_sigmatau)
m_h2 <- f_h2(n1, m_ntaustar, m_g1, m_sigmatau)

v_lossderiv <- f_lossderiv(til_beta, c(ktau + mtau, k1 + m1), m_sigma, m_h1, m_h2)

f_grad(beta, v_lossderiv, m_s)
m_hess <- f_hess(v_lossderiv, c(ktau + mtau, k1 + m1), m_s)
MASS::ginv(m_hess)
## validate the integration
v_taus <- seq(0, 1, 0.001)
cbind(colMeans(diag(v_taus) %*% spline_basis(v_taus, ktau, mtau, c(1, 0))))

tcrossprod(m_g1, m_sigmatau) %*% (n_m2 - tcrossprod(m_g2, m_sigmatau) %*% n_m2)

#examples for criteria.r
## f_loss
f_loss(til_beta, v_taustar, y, m_sigma, m_h1, m_h2)
s_criteria <- f_criteria(beta, til_beta, v_taustar, y, m_sigma, m_h1, m_h2, m_s)

#examples for descent.r
devtools::load_all()
set.seed(202104)
k1 <- 11; ktau <- 2; m1 <- 4; mtau <- 3; n <- 400
v_dim <- c(ktau + mtau, k1 + m1)
v_x <- runif(n); v_y <- log(v_x); sigma <- 0.4 * (0.5 + v_x)
v_y <- v_y + rnorm(n, mean = rep(0, n), sd = sigma)
v_taus <- seq(0.05, 0.95, 0.15); til_beta <- NULL

#v_beta <- rnorm(prod(v_dim))
#til_beta <- tilde_beta(v_beta, v_dim[2])
#init <- cq.reg(v_x, v_y, c("tesmi1", "ps"), "bfgs", c(5, 15))
#til_beta <- sum(init$coefficients[-1])
dev.new()
fit <- f_descent(v_y, v_x, c(ktau, k1), c(mtau, m1), til_beta,
                alpha = 0.1, bet = 0.1, v_smooth = c(1e-2, 1e-15, 1e-15),
                maxit = 100, eps = 1e-3, trace = seq(0.05, 0.95, 0.15))


#examples for bb-gd.r
fit <- f_bbgd(v_y, v_x, c(k1, k2), c(m1, m2), til_beta,
                v_smooth = c(10, 10, 10),
                maxit = 1, eps = 1e-4, trace = seq(0.05, 0.95, 0.15))

## f_newton
fit <- f_newton(v_y, v_x, c(k1, k2), c(m1, m2), NULL,
        alpha = 0.1, bet = 0.1, v_smooth = c(10, 10, 2),
        maxit = 100, thresh = 1e-3, trace = seq(0.05, 0.95, 0.15))


v_y <- sin(2 * pi * v_x)
v_y <- v_y + rnorm(n, mean = rep(0, n), sd = sigma)
fit <- f_descent(v_y, v_x, c(k1, k2), c(m1, m2), NULL,
                alpha = 0.2, bet = 0.1, v_smooth = c(10, 10, 2),
                maxit = 5, thresh = 1e-3, trace = v_taus)

ptm <- proc.time()# Start the clock!
fit <- f_descent(v_y, v_x, c(k1, k2), c(m1, m2), NULL,
                alpha = 0.2, bet = 0.1, v_smooth = c(10, 10, 2),
                maxit = 100, thresh = 1e-3, trace = NULL)
print(proc.time() - ptm)# Stop the clock


#examples for initial.r
devtools::load_all()
set.seed(202104)
k1 <- 4; k2 <- 4; m1 <- 4; m2 <- 4; n <- 400
v_x <- runif(n); v_y <- log(v_x); sigma <- 0.4 * (0.5 + v_x)
v_y <- v_y + rnorm(n, mean = rep(0, n), sd = sigma)
n1 <- f_n1(v_x, k1, m1)
m_sigma <- f_sigma(k1 + m1, k2 + m2)
m_s <- f_mats(k1 + m1, k2 + m2, rep(10, 3))
til_beta <- f_init(v_y, v_x, m_sigma, n1, c(k1, k2), c(m1, m2), m_s)

og_beta(til_beta, k1 + m1)

f_qscop(seq(0.05, 0.95, 0.15), n1, til_beta, m_sigma, k2, m2)

#examples for visualization.r
cq_prediction(0.1, fit, newdata = runif(10))
cq_plot(c(0.05, 0.25, 0.5, 0.75, 0.95), fit)

#examples for smoothQR.r
devtools::load_all()
set.seed(202104)
k1 <- 11; ktau <- 6; m1 <- 4; mtau <- 4; n <- 400
v_x <- runif(n); v_y <- log(v_x); sigma <- 0.4 * (0.5 + v_x)
v_y <- v_y + rnorm(n, mean = rep(0, n), sd = sigma)

#setup
v_dim <- c(ktau + mtau, k1 + m1)
n1 <- f_n1(v_x, k1, m1); beta <- rnorm(prod(v_dim[1], v_dim[2]))
m_sigmatau <- f_sigmatau(v_dim[1])
m_g1 <- f_g1(ktau, mtau)
m_g2 <- f_g2(ktau, mtau)
#h1
m_h1 <- f_smoothh1(n1, m_g1, m_g2, m_sigmatau)
m_h0 <- m_h1$h0
m_h1 <- m_h1$h1

#htau #try sparse matrix
s_ntau <- 100; v_taus <- seq(0, 1, length.out = s_ntau)
m_ndesign <- kronecker(f_ntau(v_taus, ktau, mtau), n1)
m_sigma <- f_sigma(ktau + mtau, k1 + m1)
til_beta <- tilde_beta(beta, k1 + m1)
h <- max((((log(n) + 1) / n) ^ 0.4), 0.05)

ptm <- proc.time()# Start the clock!

v_res <- f_smoothres(m_ndesign, v_y, m_sigma, til_beta, h)

print(proc.time() - ptm)# Stop the clock
#Using the Gaussian kernel
v_pres <- pnorm(v_res)
m_htau <- f_smoothhtau(m_ndesign, v_pres, s_ntau)
#deriv
loss_deriv <- f_smoothderiv(til_beta, v_dim, m_sigma, m_htau, m_h1)
#grad
m_s <- f_mats(v_dim[1], v_dim[2], rep(1, 3))
f_smoothgrad(beta, loss_deriv / n, m_s)
#loss
v_h0sig <- c(crossprod(m_h0, m_sigma))
f_smoothloss(v_res, v_pres, v_h0sig, til_beta, s_ntau, h)
f_smoothcrit(v_res, v_pres, v_h0sig, til_beta, beta, s_ntau, h, m_s)

#smoothdescent
devtools::load_all()
set.seed(202104)
k1 <- 11; ktau <- 2; m1 <- 4; mtau <- 3; n <- 400
v_x <- runif(n); v_y <- log(v_x); sigma <- 0.4 * (0.5 + v_x)
v_y <- v_y + rnorm(n, mean = rep(0, n), sd = sigma)

fit <- f_smoothbbgd(v_y, v_x, c(ktau, k1), c(mtau, m1), init = NULL,
                v_smooth = c(0.01, 0.01, 0.01),
                h = max((((log(length(v_y))+1)/length(v_y))^0.4), 0.05),
                s_ntau = 50,
                maxit = 50, eps = 1e-5, eta = 1e-8,
                trace = seq(0.05, 0.95, 0.15))

dev.new()

fit <- f_smoothdescent(v_y, v_x, c(ktau, k1), c(mtau, m1), init = NULL,
                alpha = 0.1, bet = 0.8,
                v_smooth = c(1e-3, rep(1e-2, 2)),
                h = max((((log(length(v_y))+1)/length(v_y))^0.4), 0.05),
                s_ntau = 100,
                maxit = 10, eps = 1e-5, eta = 1e-3,
                trace = seq(0.05, 0.95, 0.15))

dev.new()

ptm <- proc.time()# Start the clock!

        regfit <- cq.reg(v_x, v_y,
                        c("tedmi", "cr"),
                        "bfgs", c(10, 15),
                        span = 0.05, sp = c(2.593664e-7, 5.28e-10))

print(proc.time() - ptm)# Stop the clock





#qs.r
devtools::load_all()
set.seed(202106)
n <- 400
x <- sort(runif(n))
y <- obs_f(x, gmean = "log", xsig = "linear", derr = "norm")


regfit <- cq.reg(x, y,
                c("tedmi", "cr"),
                "bfgs", c(10, 15),
                span = 0.05, sp = c(2.593664e-7, 5.28e-10))
cq_regplot(regfit)

# quantile sheet
dev.new()
system.time(scam_fit <- qs_scam(x, y, maxit = 10, ntaus = 10,
                        arg_bs = c("tesmi1", "cr"),
                        opt = "nlm.fd", dims = c(5, 15), ords = NA,
                        analysis = F, trace = T, sp = c(1e-4, 1e-3)))
qs_loss(x, y, scam_fit)
cq_regplot(fit)



dev.new()
v_time <- system.time(gam_fit <- qs_gam(x, y, maxit = 50, ntaus = 10,
                dims = c(5, 15), ords = NA,
                analysis=F, trace = T, sp = c(0, 0), tol = 1e-3))

## count crossings
ptm <- proc.time()# Start the clock!

        ntest <- 1e3; lengthout <- 1e3;
        test_x <- seq(0, 1, length.out = ntest)
        taus <- seq(0, 1, length.out = lengthout)
        tau <- taus[cut(seq_len(ntest*lengthout), breaks = lengthout, labels = FALSE)]
        dat <- data.frame(tau = tau, x1 = rep(test_x, lengthout))
        v_pred <- predict(scam_fit, dat)
        m_pred <- matrix(c(v_pred), nrow = lengthout, byrow = T)
        # cq_regplot(gam_fit)
        # for (i in seq_along(taus)) {
        #         lines(test_x, m_pred[i, ])
        # }
        sum(apply(diff(m_pred) < 0, 1, fun <- function(x) {
                                        if (any(x)) 1
                                        else 0}))

print(proc.time() - ptm)# Stop the clock

any(diff(m_pred) < 0)

# expectile sheet
dev.new()
es_fit <- es_scam(x, y, maxit = 10, ntaus = 10,
                arg_bs = c("tedmi", "ps"),
                opt = "efs", dims = c(5, 15), ords = NA,
                analysis = F, trace = T, sp = c(10, 1e-5))

dev.new()
es_fit <- es_gam(x, y, maxit = 50, ntaus = 30,
                dims = c(5, 15), ords = NA,
                analysis = F, trace = T, sp = NULL, tol = 1e-5)

# examples in cross_valid
devtools::load_all()
set.seed(202106)
n <- 400
x <- sort(runif(n))
y <- obs_f(x, gmean = "log", xsig = "linear", derr = "norm")

v_spx <- 10^ (-8:0); v_taus <- 10^ (-1)
system.time(cv_qs <- cv_qsscam(x, y, 5,
        v_spx, v_taus, maxit = 10, ntaus = 10,
        arg_bs = c("tesmi1", "cr"),
        opt = "nlm.fd", dims = c(5, 15), tol = 1e-2))

scam_fit <- qs_scam(x, y, maxit = 10, ntaus = 10,
                        arg_bs = c("tesmi1", "cr"),
                        opt = "nlm.fd", dims = c(5, 15), ords = NA,
                        analysis = F, trace = T, sp = c(.1, 1), tol = 1e-2)
