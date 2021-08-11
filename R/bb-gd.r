#' f_bbgd
#'
#' The descent algorithm with Barzilar-borwein stepsize
#'
#' @param v_y   vector of y
#' @param v_x   vector of x
#' @param v_k0  vector of internal knots for covariates
#' @param v_m   vector of orders for covariates
#'
#' @param init initialization for $\beta$
#' @param v_smooth  vector of smoothing parameters
#' @param maxit maximum number of iterations
#' @param eps   stopping critria threshhold
#' @param eta   stopping criteria for gradient size
#' @param trace indicator whether to plot trace
#'
#'
#' @return the estimated $\beta$
#'
f_bbgd <- function(v_y, v_x, v_k0, v_m, init = NULL,
                        v_smooth = rep(1, 3),
                        maxit = 1e4, eps = 1e-3, eta = NULL, trace = NULL) {

#-------------------------------------------------------
##### set up #####
#-------------------------------------------------------
    function_call <- match.call()
    x_ord <- order(v_x)
    v_x <- v_x[x_ord]
    v_y <- v_y[x_ord]
    s_center <- mean(v_x)
    s_sd <- sd(v_x)
    v_x <- c(scale(v_x))
    iter <- 0L
    converge <- F
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
    til_beta0 <- init
    v_beta0 <- og_beta(til_beta0, k1 + m1, k2 + m2)
    v_grad0 <- 0

    #Calculate H2 (most time consuming)
    v_taustar <- f_tausearch(v_y, m_n1, til_beta0, m_sigma, k2, m2)
    m_ntau <- f_ntau(v_taustar, k2, m2)
    m_h2 <- f_h2(m_n1, m_ntau, m_g1, m_sigmatau)

    ##Loss
    s_loss <- f_loss(v_beta0, til_beta0, v_taustar, v_y, m_sigma, m_h1, m_h2, m_s)

    ##trace
    if (!is.null(trace)) {
        con_quan <- list(
                    x = v_x,
                    y = v_y,
                    center = s_center,
                    sd = s_sd,
                    v_k0 = v_k0,
                    v_m = v_m,
                    m_sigma = m_sigma,
                    n1 = m_n1,
                    til_beta = til_beta0,
                    iter = iter,
                    loss = s_loss,
                    alpha = 0,
                    bet = 0)
        cq_plot(trace, con_quan)
    }

    #Standard gradient updates
    v_deriv <- f_lossderiv(til_beta0, v_k0 + v_m, m_sigma, m_h1, m_h2)
    v_grad1 <- f_grad(v_beta0, v_deriv, m_s)
    #Calculate new coefficients
    s_step <- 0.25
    v_beta1 <- v_beta0 - s_step * v_grad1
    til_beta1 <- tilde_beta(v_beta1, k1 + m1, k2 + m2)
    new_loss <- f_loss(v_beta1, til_beta1, v_taustar, v_y, m_sigma, m_h1, m_h2, m_s)

    #Repeat
    while (maxit > iter) {
        iter <- iter + 1

        #Calculate difference
        v_delta <- v_beta1 - v_beta0
        v_gdiff <- v_grad1 - v_grad0

        #bb stepsize
        s_delgrad <- crossprod(v_delta, v_gdiff)
        step1 <- crossprod(v_delta) / s_delgrad
        step2 <- s_delgrad / crossprod(v_gdiff)
        s_step <- min(step1, step2, 100)
        if (s_step < 0) s_step <- 0.25

        #update v_beta0, til_beta0, v_grad0, s_loss
        v_beta0 <- v_beta1
        til_beta0 <- til_beta1
        v_grad0 <- v_grad1
        s_loss <- new_loss

        #Update H2 (most time consuming)
#ptm <- proc.time()# Start the clock!
        v_taustar <- f_tausearch(v_y, m_n1, til_beta1, m_sigma, k2, m2)
#print(proc.time() - ptm)# Stop the clock
        m_ntau <- f_ntau(v_taustar, k2, m2)
        m_h2 <- f_h2(m_n1, m_ntau, m_g1, m_sigmatau)

        #Calculate gradient
        v_deriv <- f_lossderiv(til_beta1, v_k0 + v_m, m_sigma, m_h1, m_h2)
        v_grad1 <- f_grad(v_beta1, v_deriv, m_s)

        #Update v_beta1, til_beta1, new_loss
        v_beta1 <- v_beta0 - s_step * v_grad1
        til_beta1 <- tilde_beta(v_beta1, k1 + m1, k2 + m2)
        new_loss <- f_loss(v_beta1, til_beta1, v_taustar, v_y, m_sigma, m_h1, m_h2, m_s)

        #Check stop criterion
        if (is.null(eta)) {
            converge <- (abs(new_loss - s_loss) <= eps * abs(s_loss))
        } else {
            converge <- (norm_2(v_grad1) / norm_2(v_beta0) <= eta)
        }

        ##trace
        if (!is.null(trace)) {
            con_quan$til_beta <- til_beta1
            con_quan$loss <- new_loss
            con_quan$iter <- iter
            cq_plot(trace, con_quan)
        }

        #If gradient size criteria
        if (converge) break

    }

#-------------------------------------------------------
##### Output #####
#-------------------------------------------------------

    out <- list(
        beta = v_beta1,
        til_beta = til_beta1,
        gradient = v_grad1,
        converge = converge,
        stepsize = s_step,
        grad_size = norm_2(v_grad1),
        iter = iter,
        loss = new_loss,
        tau = v_taustar,
        x = v_x,
        y = v_y,
        center = s_center,
        sd = s_sd,
        v_k0 = v_k0,
        v_m = v_m,
        n1 = m_n1,
        m_sigma = m_sigma,
        mat_s = m_s,
        alpha = 0,
        bet = 0,
        call = function_call)

    return(out)
}