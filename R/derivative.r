#' f_h1
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
#' v_taustar <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
#' m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
#' m_s <- f_mats(ktau + mtau, k1 + m1, rep(1, 3))
#' m_g1 <- f_g1(ktau, mtau)
#' m_g2 <- f_g2(ktau, mtau)
#' m_sigmatau <- f_sigmatau(ktau + mtau)
#' m_h1 <- f_h1(n1, m_g1, m_g2, m_sigmatau)
#' m_h2 <- f_h2(n1, m_ntaustar, m_g1, m_sigmatau)
#' }
f_h1 <- function(n1, g1, g2, sigmatau) {
    #form N^{m_tau+2}(\tau)|_0^1
    nm2 <- rep(0, nrow(sigmatau)); nm2[nrow(sigmatau)] <- 1
    #Calculate $H_1$
    h_1 <- tcrossprod(g1, sigmatau) %*% tcrossprod(g2, sigmatau)
    h_1 <- -h_1 %*% nm2
    h_1 <- kronecker(h_1, t(n1))
    h_1 <- rowSums(h_1)
    return(h_1)
}


#' f_h2
#'
#' Construct the $K_\tauK_1$ vector $H_2$
#'
#' @param n1 matrix n1 output from f_n1
#' @param ntaustar matrix ntau output from f_ntaustar
#' @param g1 $G_1$ diagonal matrix
#' @param sigmatau $\Sigma_tau$ diagonal matrix
#'
#'
#' @return $K_\tau K_1$ vector $H_2$
#'
#' @import mgcv
#'
#' @examples
#' \dontrun{
#' v_taustar <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
#' m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
#' m_s <- f_mats(ktau + mtau, k1 + m1, rep(1, 3))
#' m_g1 <- f_g1(ktau, mtau)
#' m_g2 <- f_g2(ktau, mtau)
#' m_sigmatau <- f_sigmatau(ktau + mtau)
#' m_h1 <- f_h1(n1, m_g1, m_g2, m_sigmatau)
#' m_h2 <- f_h2(n1, m_ntaustar, m_g1, m_sigmatau)
#' }
f_h2 <- function(n1, ntaustar, g1, sigmatau) {
    h_2 <- tcrossprod(g1, sigmatau) %*% t(ntaustar)
#ptm <- proc.time()# Start the clock!
    h_2 <- mgcv::tensor.prod.model.matrix(list(t(h_2), n1))
    h_2 <- colSums(h_2)
#print(proc.time() - ptm)# Stop the clock

    return(h_2)
}

#' f_lossderiv
#'
#' Construct the $K_\tau K_1$ gradient vector for the loss
#'
#' @param til_beta the current estimation fro $\tilde{\beta}$
#' @param v_dim vector contain the dimension of covariates
#'             $(K_\tau, K_1, K_2, K_3, \cdots)$
#'
#' @param m_sigma $K_\tau K_1K_2 \times K_1K_2 \Sigma$ diagonal block matix.
#' @param h1 $K_\tau K_1K_2$ spline vector for the first covariate
#' @param h2 $K_\tau K_1K_2$ spline vector for the tau covariate
#'
#'
#' @return the $K_\tau K_1$ gradient vector for the loss
#'
#' @examples
#' \dontrun{
#' v_taustar <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
#' m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
#' m_s <- f_mats(ktau + mtau, k1 + m1, rep(1, 3))
#' m_g1 <- f_g1(ktau, mtau)
#' m_g2 <- f_g2(ktau, mtau)
#' m_sigmatau <- f_sigmatau(ktau + mtau)
#' m_h1 <- f_h1(n1, m_g1, m_g2, m_sigmatau)
#' m_h2 <- f_h2(n1, m_ntaustar, m_g1, m_sigmatau)
#' v_lossderiv <- f_lossderiv(til_beta, c(ktau + mtau, k1 + m1), m_sigma, m_h1, m_h2)
#' }
f_lossderiv <- function(til_beta, v_dim, m_sigma, h1, h2) {
    m_c <- f_matc(v_dim, til_beta)
    lossderiv <- - m_c %*% crossprod(m_sigma, (h1 + h2))
    return(c(lossderiv))
}

#' f_grad
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
#' v_taustar <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
#' m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
#' m_s <- f_mats(ktau + mtau, k1 + m1, rep(1, 3))
#' m_g1 <- f_g1(ktau, mtau)
#' m_g2 <- f_g2(ktau, mtau)
#' m_sigmatau <- f_sigmatau(ktau + mtau)
#' m_h1 <- f_h1(n1, m_g1, m_g2, m_sigmatau)
#' m_h2 <- f_h2(n1, m_ntaustar, m_g1, m_sigmatau)
#' v_lossderiv <- f_lossderiv(til_beta, c(ktau + mtau, k1 + m1), m_sigma, m_h1, m_h2)
#' f_grad(beta, v_lossderiv, m_s)
#' }
f_grad <- function(beta, loss_deriv, s) {
    loss_deriv + c(s %*% beta)
}

#' f_hess
#'
#' Construct the $K_\tau K_1K_2 \times K_\tau K_1K_2$ Hessian matrix for $\beta$
#'
#' @param loss_deriv the derivative of the loss
#' @param v_dim vector contain the dimension of covariates
#'             $(K_\tau, K_1, K_2, K_3, \cdots)$
#' @param s $K_\tau K_1K_2 \times K_\tau K_1K_2$ $S$ matrix
#'
#'
#' @return the $K_\tau K_1K_2$ gradient vector for $\beta$
#'
#' @examples
#' \dontrun{
#' v_taustar <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
#' m_ntaustar <- f_ntaustar(v_taustar, ktau, mtau)
#' m_s <- f_mats(ktau + mtau, k1 + m1, rep(1, 3))
#' m_g1 <- f_g1(ktau, mtau)
#' m_g2 <- f_g2(ktau, mtau)
#' m_sigmatau <- f_sigmatau(ktau + mtau)
#' n_m2 <- rep(0, ktau + mtau)
#' n_m2[ktau + mtau] <- 1
#' m_h1 <- f_h1(n1, m_g1, m_g2, m_sigmatau, n_m2)
#' m_h2 <- f_h2(n1, m_ntaustar, m_g1, m_sigmatau)
#' v_lossderiv <- f_lossderiv(til_beta, c(ktau + mtau, k1 + m1), m_sigma, m_h1, m_h2)
#' f_grad(beta, v_lossderiv, m_s)
#' m_hess <- f_hess(v_lossderiv, c(ktau + mtau, k1 + m1), m_s)
#' }
f_hess <- function(loss_deriv, v_dim, s) {
    v_idx <- seq_len(v_dim[2])
    loss_deriv[v_idx] <- 0
    return(diag(loss_deriv) + s)
}