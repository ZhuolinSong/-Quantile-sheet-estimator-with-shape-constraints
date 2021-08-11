#' f_descent
#'
#' The descent algorithm using back tracking line search stepsize
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
#' @param eps   stopping critria threshhold
#' @param eta   stopping criteria for gradient size
#' @param trace indicator whether to plot trace
#'
#'
#' @return the estimated $\beta$
#'
f_descent <- function(v_y, v_x, v_k0, v_m, init = NULL,
                        alpha = 0.3, bet = 0.5, v_smooth = rep(1, 3),
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
    s_n <- length(v_y)
    #Define dimensions
    k1 <- v_k0[2]; ktau <- v_k0[1]; m1 <- v_m[2]; mtau <- v_m[1]
    v_dim <- v_k0 + v_m
    #Require
    m_n1 <- f_n1(v_x, k1, m1)
    m_g1 <- f_g1(ktau, mtau)
    m_g2 <- f_g2(ktau, mtau)
    m_sigmatau <- f_sigmatau(v_dim[1])
    m_sigma <- f_sigma(v_dim[1], v_dim[2])
    m_s <- f_mats(v_dim[1], v_dim[2], v_smooth)

    #Precalculate h1
    m_h1 <- f_h1(m_n1, m_g1, m_g2, m_sigmatau)

    #Initialization
    if (is.null(init)) {
        init <- f_init(v_y, v_x, m_sigma, m_n1, v_k0, v_m, m_s)
    }
    til_beta <- init
    v_beta <- og_beta(til_beta, v_dim[1], v_dim[2])
    v_grad <- 0

    #Calculate H2 (most time consuming)
    v_taustar <- f_tausearch(v_y, m_n1, til_beta, m_sigma, ktau, mtau)
    m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
    m_h2 <- f_h2(m_n1, m_ntaustar, m_g1, m_sigmatau)

    ##Loss
    s_loss <- f_criteria(v_beta, til_beta, v_taustar, v_y, m_sigma, m_h1, m_h2, m_s)

    ##trace
    ##trace
    if (!is.null(trace)) {
        v_taus <- trace
        con_quan <- list(x = v_x, y = v_y, center = s_center, sd = s_sd,
                    v_k0 = v_k0, v_m = v_m, m_sigma = m_sigma, n1 = m_n1,
                    til_beta = til_beta)
        sub <- paste("Iter:", iter, ", Error=", round(s_loss, 3))
        main <- paste("alpha=", alpha, "beta=", bet,
                "Smoothing:", v_smooth[1:2], "dim:", v_dim)
        cq_plot(trace, con_quan, sub = sub, main = main)
    }

    #Repeat
    while (maxit > iter) {
        iter <- iter + 1

        #Calculate gradient
        v_deriv <- f_lossderiv(til_beta, v_k0 + v_m, m_sigma, m_h1, m_h2) / s_n
        v_grad <- f_grad(v_beta, v_deriv, m_s)
        grad_size <-  sum(v_grad ^ 2)

        #backtracking line search
        t <- 1
        repeat{
            new_beta <- v_beta - t * v_grad
            new_tilbeta <- tilde_beta(new_beta, v_dim[2])
            new_loss <- f_criteria(new_beta, new_tilbeta, v_taustar, v_y, m_sigma, m_h1, m_h2, m_s)
            t <- bet * t
            if (new_loss <= (s_loss - alpha * t * grad_size) || t < 1e-4) {
                break
            }
        }

        #Update v_beta
        v_beta <- new_beta
        til_beta <- new_tilbeta

        #Update H2 (most time consuming)
#ptm <- proc.time()# Start the clock!
        v_taustar <- f_tausearch(v_y, m_n1, til_beta, m_sigma, ktau, mtau)
#print(proc.time() - ptm)# Stop the clock
        m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
        m_h2 <- f_h2(m_n1, m_ntaustar, m_g1, m_sigmatau)

        ##trace
        if (!is.null(trace)) {
            con_quan$til_beta <- til_beta
            sub <- paste("Iter:", iter, ", Error=", round(new_loss, 3))
            cq_plot(trace, con_quan, sub = sub, main = main)
        }

        #Check stop criterion
        if (is.null(eta)) {
            converge <- (abs(new_loss - s_loss) <= eps * abs(s_loss))
        } else {
            converge <- (sqrt(grad_size) / norm_2(v_beta) <= eta)
        }

        s_loss <- new_loss
        #If gradient size criteria
        if (is.na(converge)) {
            converge <- F
            break
        }
        if (converge) break
    }

#-------------------------------------------------------
##### Output #####
#-------------------------------------------------------
    out <- list(
        beta = v_beta,
        til_beta = til_beta,
        gradient = v_grad,
        converge = converge,
        grad_size = sum(v_grad ^ 2),
        iter = iter,
        loss = s_loss,
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
        alpha = alpha,
        bet = bet,
        call = function_call)

    return(out)
}