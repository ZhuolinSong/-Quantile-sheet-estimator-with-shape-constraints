#' f_smoothh1
#'
#' Construct the $K_\tau K_1$ vector $h_1$
#'
#' @param n1 matrix n1 output from f_n1
#' @param g1 $G_1$ diagonal matrix
#' @param g2 $G_2$ diagonal matrix
#' @param sigmatau $\Sigma_tau$ diagonal matrix
#'
#'
#' @return $K_\tau K_1$ vector $H_1$
#'
#' @examples
#' \dontrun{
#' devtools::load_all()
#'k1 <- 4; ktau <- 4; m1 <- 4; mtau <- 4; n <- 400
#' v_dim <- c(ktau + mtau, k1 + m1)
#' v_x <- runif(n); v_y <- log(v_x); sigma <- 0.4 * (0.5 + v_x)
#' n1 <- f_n1(v_x, k1, m1); beta <- rnorm(prod(v_dim[1], v_dim[2]))
#' m_sigmatau <- f_sigmatau(v_dim[1])
#' m_g1 <- f_g1(ktau, mtau)
#' m_g2 <- f_g2(ktau, mtau)
#' #h1
#' m_h1 <- f_smoothh1(n1, m_g1, m_g2, m_sigmatau)
#' m_h1$h1
#' m_h1$h0
#' }
f_smoothh1 <- function(n1, g1, g2, sigmatau) {
    #form N^{m_tau+1}(1) and N^{m_tau+2}(\tau)|_0^1
    nm1 <- rep(0, nrow(sigmatau)); nm1[nrow(sigmatau)] <- 1
    h_1 <- nm1 - tcrossprod(g2, sigmatau) %*% nm1
    h_1 <- tcrossprod(g1, sigmatau) %*% h_1
    h_1 <- kronecker(h_1, t(n1))
    h_1 <- rowSums(h_1)

    h_0 <- tcrossprod(g1, sigmatau) %*% nm1
    h_0 <- kronecker(h_0, t(n1))
    h_0 <- 0.5 * rowSums(h_0) - h_1
    return(list(h1 = h_1,
                h0 = h_0))
}



#' f_smoothres
#'
#' Construct the $n n_tau$ residual vector
#'
#' @param ndesign   design matrix $n n_\tau \times K_\tau K_1$ $N(\tau, x)$
#' @param v_y       vector of observations y_i
#' @param sigma  $\Sigma$ diagonal matrix
#' @param til_beta  tilde_beta coefficients
#'
#' @param h         bandwidth value
#' @return $n n_tau$ residual vector
#'
#'
#' @examples
#' \dontrun{
#' s_ntau <- 500; v_taus <- seq(0, 1, length.out = s_ntau)
#' m_ndesign <- kronecker(f_ntau(v_taus, ktau, mtau), n1)
#' m_sigma <- f_sigma(ktau + mtau, k1 + m1)
#' til_beta <- tilde_beta(beta, k1 + m1)
#' h <- max((((log(n) + 1) / n) ^ 0.4), 0.05)
#' v_res <- f_smoothres(m_ndesign, v_y, m_sigma, til_beta)
#' }
f_smoothres <- function(ndesign, v_y, sigma, til_beta, h) {
    #Calculate the residual
    c(ndesign %*% sigma %*% til_beta - v_y) / h
}


#' f_smoothhtau
#'
#' Construct the $K_\tau K_1$ vector $h_\tau$
#'
#' @param ndesign   design matrix $n n_\tau \times K_\tau K_1$ $N(\tau, x)$
#' @param pres  error function of the $K_\tau n_tau$ residual vector devided by h
#' @param s_ntau    number of taus values
#'
#'
#' @return $K_\tau K_1$ vector $h_\tau$
#'
#'
#' @examples
#' \dontrun{
#' s_ntau <- 500; v_taus <- seq(0, 1, length.out = s_ntau)
#' m_ndesign <- kronecker(f_ntau(v_taus, ktau, mtau), n1)
#' m_sigma <- f_sigma(ktau + mtau, k1 + m1)
#' til_beta <- tilde_beta(beta, k1 + m1)
#' h <- max((((log(n) + 1) / n) ^ 0.4), 0.05)
#'
#' ptm <- proc.time()# Start the clock!
#' v_res <- f_smoothres(m_ndesign, v_y, m_sigma, til_beta)
#' print(proc.time() - ptm)# Stop the clock
#' #Using the Gaussian kernel
#' v_pres <- pnorm(v_res)
#' m_htau <- f_smoothhtau(m_ndesign, v_pres, h, s_ntau)
#' }
f_smoothhtau <- function(ndesign, pres, s_ntau) {
    c(crossprod(pres, ndesign) / s_ntau)
}



#' f_smoothderiv
#'
#' Construct the $K_\tau K_1$ gradient vector for the loss
#'
#' @param til_beta the current estimation fro $\tilde{\beta}$
#' @param v_dim vector contain the dimension of covariates
#'             $(K_\tau, K_1, K_2, K_3, \cdots)$
#'
#' @param m_sigma $K_\tau K_1K_2 \times K_1K_2 \Sigma$ diagonal block matix.
#' @param h1 $K_\tau K_1K_2$ vector
#' @param h2 $K_\tau K_1K_2$ vector
#'
#'
#' @return the $K_\tau K_1$ gradient vector for the loss
#'
#' @examples
#' \dontrun{
#' f_smoothderiv(til_beta, v_dim, m_sigma, m_htau, m_h1)
#' }
f_smoothderiv <- function(til_beta, v_dim, m_sigma, htau, h1) {
    m_c <- f_matc(v_dim, til_beta)
    lossderiv <- m_c %*% crossprod(m_sigma, (htau - h1))
    return(c(lossderiv))
}

#' f_smoothgrad
#'
#' Construct the $K_\tau K_1K_2$ gradient vector for $\beta$
#'
#' @param beta the current estimation fro $\beta$
#' @param loss_deriv the derivative of the loss
#' @param s $K_\tau K_1K_2 \times K_1K_2$ $S$ matrix
#'
#'
#' @return the $K_\tau K_1K_2$ gradient vector for $\beta$
#'
#' @examples
#' \dontrun{
#' f_grad(beta, v_lossderiv, m_s)
#' }
f_smoothgrad <- function(beta, loss_deriv, s) {
    loss_deriv + c(s %*% beta)
}

#' f_smoothloss
#'
#' Return the quantile loss criteria
#'
#' @param res   $K_\tau n_tau$ residual vector
#' @param pres  error function of the $K_\tau n_tau$ residual vector devided by h
#' @param h0sig    $K_\tau K_1$ vector $H_0$ times m_sigma
#' @param til_beta  estimation for $\tilde{\beta}$
#'
#' @param s_ntau    number of taus values
#' @param h bandwidth value
#' @param s_n number of observation
#'
#' @return the smooth quantile loss
#'
f_smoothloss <- function(res, pres, h0sig, til_beta,
                        s_ntau, h, s_n) {
    loss <- crossprod(h0sig, til_beta)
    resloss <- h / 2 * sqrt(2 / pi) * sum(exp(- (- res / 2)^2))
    resloss <- resloss - crossprod(res, (1 - 2 * pres))
    loss <- loss + resloss / s_ntau
    return(c(loss) / s_n)
}


#' f_smoothcrit
#'
#' Return the smooth QR criteria
#'
#' @param res   $n n_tau$ residual vector
#' @param pres  error function of the $K_\tau n_tau$ residual vector devided by h
#' @param h0sig    $K_\tau K_1$ vector $H_0$ times m_sigma
#' @param til_beta  estimation for $\tilde{\beta}$
#' @param beta the current estimation for $\beta$
#'
#' @param s_ntau    number of taus values
#' @param h bandwidth value
#' @param s_n number of observation
#' @param s $K_1K_2 \times K_1K_2$ $S$ matrix
#'
#' @return the smooth quantile criteria
#'
f_smoothcrit <- function(res, pres, h0sig, til_beta, beta,
                        s_ntau, h, s_n, s) {
    loss <- f_smoothloss(res, pres, h0sig, til_beta, s_ntau, h, s_n)
    criteria <- loss + crossprod(beta, s %*% beta)
    return(c(criteria))
}

#' f_smoothbbgd
#'
#' The descent algorithm using bb stepsize
#'
#' @param v_y   vector of y
#' @param v_x   vector of x
#' @param v_k0  vector of internal knots for covariates
#' @param v_m   vector of orders for covariates
#'
#' @param init initialization for $\beta$
#' @param h     bandwidth value
#' @param s_ntau    number of taus for numerical integration
#' @param v_smooth  vector of smoothing parameters
#' @param maxit maximum number of iterations
#' @param eps   stopping critria threshhold
#' @param eta   stopping criteria for gradient size
#' @param trace indicator whether to plot trace
#'
#'
#' @return the estimated $\beta$
#'
f_smoothbbgd <- function(v_y, v_x, v_k0, v_m, init = NULL,
                        v_smooth = rep(1, 3),
                        h=max((((log(length(v_y))+1)/length(v_y))^0.4), 0.05),
                        s_ntau = 200,
                        maxit = 1e4, eps = 1e-3, eta = 1e-3,
                        trace = NULL, echostep = F) {

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
    # Set up design matrix
    v_taus <- seq(0, 1, length.out = s_ntau)
    m_ndesign <- kronecker(f_ntau(v_taus, ktau, mtau), m_n1)

    #Precalculate h1, h0, h0 * Sigma
    m_h1 <- f_smoothh1(m_n1, m_g1, m_g2, m_sigmatau)
    m_h0 <- m_h1$h1
    m_h1 <- m_h1$h1
    v_h0sig <- c(crossprod(m_h0, m_sigma))

    #Initialization
    if (is.null(init)) {
        init <- f_init(v_y, v_x, m_sigma, m_n1, v_k0, v_m, m_s)
    }
    til_beta0 <- init
    v_beta0 <- og_beta(til_beta0, v_dim[1], v_dim[2])
    v_grad0 <- 0

    #Calculate v_res (most time consuming)
    v_res <- f_smoothres(m_ndesign, v_y, m_sigma, til_beta0, h)
    #Using the Gaussian kernel for h_\tau
    v_pres <- pnorm(v_res)
    m_htau <- f_smoothhtau(m_ndesign, v_pres, s_ntau)

    ##Loss
    s_loss <- f_smoothcrit(v_res, v_pres, v_h0sig, til_beta0, v_beta0,
                            s_ntau, h, s_n, m_s)

    ##trace
    if (!is.null(trace)) {
        v_taus <- trace
        con_quan <- list(x = v_x, y = v_y, center = s_center, sd = s_sd,
                    v_k0 = v_k0, v_m = v_m, m_sigma = m_sigma, n1 = m_n1,
                    til_beta = til_beta0)
        sub <- paste("Iter:", iter, ", Error=", round(s_loss, 3))
        main <- paste("Smoothing:", v_smooth[1:2], "dim:", v_dim)
        cq_plot(trace, con_quan, sub = sub, main = main)
    }

    #Standard gradient updates
    v_deriv <- f_smoothderiv(til_beta0, v_dim, m_sigma, m_htau, m_h1)
    v_grad1 <- f_smoothgrad(v_beta0, v_deriv / s_n, m_s)

    #Calculate new coefficients
    v_beta1 <- v_beta0 - v_grad1
    til_beta1 <- tilde_beta(v_beta1, v_dim[2])
    new_loss <- f_smoothcrit(v_res, v_pres, v_h0sig, til_beta1, v_beta1,
                            s_ntau, h, s_n, m_s)

    #Repeat
    while (maxit > iter) {
        iter <- iter + 1

        #Calculate difference
        v_delta <- v_beta1 - v_beta0
        v_gdiff <- v_grad1 - v_grad0

        #bb stepsize
        s_delgrad <- crossprod(v_delta, v_gdiff)
        if (s_delgrad < 0) {
            if (echostep) print(s_delgrad)
            s_step <- 1
        } else  {
            step1 <- crossprod(v_delta) / s_delgrad
            step2 <- s_delgrad / crossprod(v_gdiff)
            s_step <- min(step1, step2, 10)
            if (echostep) print(s_step)
        }

        #update v_beta0, til_beta0, v_grad0, s_loss
        v_beta0 <- v_beta1
        til_beta0 <- til_beta1
        v_grad0 <- v_grad1
        s_loss <- new_loss

        #Calculate v_res (most time consuming)
        v_res <- f_smoothres(m_ndesign, v_y, m_sigma, til_beta0, h)
        #Using the Gaussian kernel for h_\tau
        v_pres <- pnorm(v_res)
        m_htau <- f_smoothhtau(m_ndesign, v_pres, s_ntau)

        #Calculate gradient
        v_deriv <- f_smoothderiv(til_beta0, v_dim, m_sigma, m_htau, m_h1)
        v_grad1 <- f_smoothgrad(v_beta0, v_deriv / s_n, m_s)

        #Update v_beta1, til_beta1
        v_beta1 <- v_beta0 - s_step * v_grad1
        til_beta1 <- tilde_beta(v_beta1, v_dim[2])
        new_loss <- f_smoothcrit(v_res, v_pres, v_h0sig, til_beta1, v_beta1,
                            s_ntau, h, s_n, m_s)

        ##trace
        if (!is.null(trace)) {
            con_quan$til_beta <- til_beta1
            sub <- paste("Iter:", iter, ", Error=", round(new_loss, 3))
            cq_plot(trace, con_quan, sub = sub, main = main)
        }

        #Check stop criterion
        if (is.null(eta)) {
            converge <- (abs(new_loss - s_loss) <= eps * abs(s_loss))
        } else {
            converge <- (norm_2(v_grad1) / norm_2(v_beta0) <= eta)
        }

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
        beta = v_beta1,
        til_beta = til_beta1,
        gradient = v_grad1,
        converge = converge,
        grad_size = sum(v_grad1 ^ 2),
        iter = iter,
        loss = s_loss,
        x = v_x,
        y = v_y,
        center = s_center,
        sd = s_sd,
        v_k0 = v_k0,
        v_m = v_m,
        n1 = m_n1,
        m_sigma = m_sigma,
        mat_s = m_s,
        call = function_call)

    return(out)
}





#' f_smoothdescent
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
f_smoothdescent <- function(v_y, v_x, v_k0, v_m, init = NULL,
                        alpha = 0.3, bet = 0.8, v_smooth = rep(1, 3),
                        h=max((((log(length(v_y))+1)/length(v_y))^0.4), 0.05),
                        s_ntau = 200,
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
    # Set up design matrix
    v_taus <- seq(0, 1, length.out = s_ntau)
    m_ndesign <- kronecker(f_ntau(v_taus, ktau, mtau), m_n1)

    #Precalculate h1, h0, h0 * Sigma
    m_h1 <- f_smoothh1(m_n1, m_g1, m_g2, m_sigmatau)
    m_h0 <- m_h1$h1
    m_h1 <- m_h1$h1
    v_h0sig <- c(crossprod(m_h0, m_sigma))

    #Initialization
    if (is.null(init)) {
        init <- f_init(v_y, v_x, m_sigma, m_n1, v_k0, v_m, m_s)
    }
    til_beta <- init
    v_beta <- og_beta(til_beta, v_dim[1], v_dim[2])
    v_grad <- 0

    #Calculate v_res (most time consuming)
    v_res <- f_smoothres(m_ndesign, v_y, m_sigma, til_beta, h)
    #Using the Gaussian kernel for h_\tau
    v_pres <- pnorm(v_res)
    m_htau <- f_smoothhtau(m_ndesign, v_pres, s_ntau)

    ##Loss
    s_loss <- f_smoothcrit(v_res, v_pres, v_h0sig, til_beta, v_beta,
                            s_ntau, h, s_n, m_s)

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
        v_deriv <- f_smoothderiv(til_beta, v_dim, m_sigma, m_htau, m_h1)
        v_grad <- f_smoothgrad(v_beta, v_deriv / s_n, m_s)
        grad_size <-  sum(v_grad ^ 2)
        #backtracking line search
        t <- 1
        repeat{
            new_beta <- v_beta - t * v_grad
            new_tilbeta <- tilde_beta(new_beta, v_dim[2])
            new_loss <- f_smoothcrit(v_res, v_pres, v_h0sig,
                                new_tilbeta, new_beta,
                                s_ntau, h, s_n, m_s)
            t <- bet * t
            if (new_loss <= (s_loss - alpha * t * grad_size) || t < 1e-4) {
                break
            }
        }

        #Update v_beta
        v_beta <- new_beta
        til_beta <- new_tilbeta

        #Calculate v_res (most time consuming)
        v_res <- f_smoothres(m_ndesign, v_y, m_sigma, til_beta, h)
        #Using the Gaussian kernel for h_\tau
        v_pres <- pnorm(v_res)
        m_htau <- f_smoothhtau(m_ndesign, v_pres, s_ntau)

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