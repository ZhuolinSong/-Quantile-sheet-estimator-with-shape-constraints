
#' f_qscop
#'
#' Estimate the quantile given specific $(\tau, x)$
#'
#' @param tau conditional probability
#' @param n1 the first covariates
#' @param til_beta the current estimation fro $\tilde{\beta}$
#' @param m_sigma $K_\tau K_1 \times K_\tau K_1 \Sigma$ diagonal block matix.
#' @param ktau internal knots of tau
#' @param mtau order for tau spline
#'
#'
#' @return An estimated the quantile given specific $(\tau, x)$
#'
#' @examples
#' \dontrun{
#' k1 <- 4; ktau <- 4; m1 <- 4; mtau <- 4; y <- rnorm(200)
#' n1 <- f_n1(1:200, k1, m1); beta <- rnorm(prod(ktau + mtau, k1 + m1))
#' til_beta <- tilde_beta(beta, k1 + m1)
#' m_sigma <- f_sigma(ktau + mtau, k1 + m1)
#' f_qscop(seq(0, 1, 0.1), n1, til_beta, m_sigma, ktau, mtau)
#' f_qscop(0, n1, til_beta, m_sigma, ktau, mtau) < f_qscop(1, n1, til_beta, m_sigma, ktau, mtau)

#' }

f_qscop <- function(tau, n1, til_beta, m_sigma, ktau, mtau) {
    m_n <- kronecker(spline_basis(tau, ktau, mtau, c(1, 0)), n1)
    out <- m_n %*% m_sigma %*% til_beta
    return(c(out))
}



#' f_qzeros
#'
#' zeros function to find the root $\tau_i$ for each observs
#'
#' @param tau conditional probability
#' @param idx index of the observations
#' @param v_y vector of y
#' @param n1 spline basis for the first covariates
#' @param ... other 4 args for f_qscop in order (til_beta, m_sigma, ktau, mtau)
#'
#'
#' @return Difference bw response and estimation $Q(\tau, x_i) - y_i$
#'
#' @examples
#' \dontrun{
#' k1 <- 4; ktau <- 4; m1 <- 4; mtau <- 4; y <- rnorm(200)
#' n1 <- f_n1(1:200, k1, m1); beta <- rnorm(prod(ktau + mtau, k1 + m1))
#' til_beta <- tilde_beta(beta, k1 + m1)
#' m_sigma <- f_sigma(ktau + mtau, k1 + m1)
#' f_qzeros(0.5, 2:3, y, n1, til_beta, m_sigma, ktau, mtau)
#' }
f_qzeros <- function(tau, idx, v_y, n1, ...) {
    n1 <- n1[idx, ]
    if (length(idx) == 1) {
        n1 <- t(n1)
    }
    s_dif <- f_qscop(tau, n1, ...) - v_y[idx]
    if ((tau == 0 && s_dif > 0) || (tau == 1 && s_dif < 0)) {
        return(0)
    }

    return(s_dif)
}


#' f_tausearch
#'
#' search the $\tau^*$ using the unitroot function
#'
#' @param v_y vector of y
#' @param n1 spline basis for the first covariates
#' @param ... other 4 args for f_qscop in order (til_beta, m_sigma, k2, m2)
#'
#'
#' @return $\tau^*$
#'
#' \dontrun{
#' k1 <- 4; ktau <- 4; m1 <- 4; mtau <- 4; y <- rnorm(200)
#' n1 <- f_n1(1:200, k1, m1); beta <- rnorm(prod(ktau + mtau, k1 + m1))
#' til_beta <- tilde_beta(beta, k1 + m1)
#' m_sigma <- f_sigma(ktau + mtau, k1 + m1)
#' f_qzeros(0.5, 2:3, y, n1, til_beta, m_sigma, ktau, mtau)
#' v_taus <- f_tausearch(y, n1, til_beta, m_sigma, ktau, mtau)
#' sapply(1:200, function(x) {
#'      f_qscop(v_taus[x], t(n1[x, ]), til_beta, m_sigma, ktau, mtau)
#' }) - y
#' }
f_tausearch <- function(v_y, n1, ...) {
    taustar <- sapply(seq_along(v_y), function(x) {
                    uniroot(f_qzeros,
                            idx = x,
                            v_y = v_y,
                            n1 = n1,
                            ...,
                            interval = c(0, 1))$root})
    return(taustar)
}