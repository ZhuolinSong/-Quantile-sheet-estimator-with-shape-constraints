#' f_newton
#'
#' The newton algorithm to minimize the criteria
#'
#' @param v_y   vector of y
#' @param v_x   vector of x
#' @param v_k0  vector of internal knots for covariates
#' @param v_m   vector of orders for covariates
#'
#' @param init initialization for $\beta$
#' @param alpha threshold $\alpha$, within $(0, 0.5)$.
#' @param bet   shrinkage rate $\bet$, within $(0, 1)$.
#' @param v_smooth  vector of smoothing parameters
#' @param maxit maximum number of iterations
#' @param thresh    stopping critria threshhold
#' @param trace indicator whether to plot trace
#'
#' @import MASS
#'
#' @return the estimated $\beta$
#'
f_newton <- function(v_y, v_x, v_k0, v_m, init,
                        alpha = 0.2, bet = 0.1, v_smooth = rep(1, 3),
                        maxit = 1e4, thresh = 1e-3, trace = NULL) {

#-------------------------------------------------------
##### set up #####
#-------------------------------------------------------
    function_call <- match.call()
    x_ord <- order(v_x)
    v_x <- v_x[x_ord]
    v_y <- v_y[x_ord]
    iter <- 0L
    #Define dimensions
    k1 <- v_k0[1]; k2 <- v_k0[2]; m1 <- v_m[1]; m2 <- v_m[2]
    #Require
    m_n1 <- f_n1(v_x, k1, m1)
    m_g1 <- f_g1(k2, m2)
    m_g2 <- f_g2(k2, m2)
    m_sigmatau <- f_sigmatau(k2 + m2)
    m_sigma <- f_sigma(k1 + m1, k2 + m2)
    m_s <- f_mats(k1 + m1, k2 + m2, v_smooth)
    n_m2 <- rep(0, k2 + m2); n_m2[k2 + m2] <- 1

    #Precalculate
    m_h1 <- f_h1(m_n1, m_g1, m_g2, m_sigmatau, n_m2)

    #Initialization
    if (is.null(init)) {
        init <- f_init(v_y, v_x, m_sigma, m_n1, v_k0, v_m, m_s)
    }
    til_beta <- init
    v_beta <- og_beta(til_beta, k1 + m1, k2 + m2)

    #Calculate H2
    v_taustar <- f_tausearch(v_y, m_n1, til_beta, m_sigma, k2, m2)
    m_ntau <- f_ntau(v_taustar, k2, m2)
    m_h2 <- f_h2(m_n1, m_ntau, m_g1, m_sigmatau)

    #Loss
    s_loss <- f_loss(v_beta, til_beta, v_taustar, v_y, m_sigma, m_h1, m_h2, m_s)
    ##trace
    if (!is.null(trace)) {
        plot(v_x, v_y,
            type = "p", pch = 16, cex = 0.8, col = "grey",
            xlim = grDevices::extendrange(v_x, f = .5),
            ylim = grDevices::extendrange(v_y, f = .5),
            xlab = "", ylab = "",
            sub = paste("Iter:", iter, ", Error=", round(Inf, 3))
        )
        for (i in seq_along(trace)) {
            v_med <- f_qscop(trace[i], m_n1, til_beta, m_sigma, k2, m2)
            lines(v_x, v_med, col = i + 1)
        }
    }


    #Repeat
    while (maxit > iter) {
        iter <- iter + 1

        #Calculate gradient
        v_deriv <- f_lossderiv(til_beta, v_k0 + v_m, m_sigma, m_h1, m_h2)
        v_grad <- f_grad(v_beta, v_deriv, m_s)
        #Calculate Hessian
        m_hess <- f_hess(v_deriv, v_k0 + v_m, m_s)

        ##Loss
        t <- 1
        #Line search
        step_size <-  c(MASS::ginv(m_hess) %*% v_grad)
        grad_size <- alpha * crossprod(v_grad, step_size)
        repeat{
            new_beta <- v_beta - t * step_size
            new_tilbeta <- tilde_beta(new_beta, k1 + m1, k2 + m2)
            new_loss <- f_loss(new_beta, new_tilbeta, v_taustar, v_y, m_sigma, m_h1, m_h2, m_s)
            t <- bet * t
            if (new_loss <= (s_loss -  t * grad_size)) {
                break
            }
        }

        #Update v_beta
        v_beta <- new_beta
        til_beta <- new_tilbeta

        #Update H2 (most time consuming)
#ptm <- proc.time()# Start the clock!
        v_taustar <- f_tausearch(v_y, m_n1, til_beta, m_sigma, k2, m2)
#print(proc.time() - ptm)# Stop the clock
        m_ntau <- f_ntau(v_taustar, k2, m2)
        m_h2 <- f_h2(m_n1, m_ntau, m_g1, m_sigmatau)

        ##trace
        if (!is.null(trace)) {
            plot(v_x, v_y,
                type = "p", pch = 16, cex = 0.8, col = "grey",
                xlim = grDevices::extendrange(v_x, f = .5),
                ylim = grDevices::extendrange(v_y, f = .5),
                xlab = "", ylab = "",
                sub = paste("Iter:", iter, "Error=", round(s_loss, 3))
            )
            for (i in seq_along(trace)) {
                v_med <- f_qscop(trace[i], m_n1, til_beta, m_sigma, k2, m2)
                lines(v_x, v_med, col = i + 1)
            }
        }
        #Check stop criterion
        converge <- (abs(new_loss - s_loss) <= thresh * s_loss)
        s_loss <- new_loss
        if (converge) break

    }

#-------------------------------------------------------
##### Output #####
#-------------------------------------------------------
    out <- list(
        beta = v_beta,
        til_beta = til_beta,
        converge = converge,
        iter = iter,
        loss = s_loss,
        tau = v_taustar,
        x = v_x,
        y = v_y,
        n1 = m_n1,
        m_sigma = m_sigma,
        mat_s = m_s,
        call = function_call)

    return(out)
}