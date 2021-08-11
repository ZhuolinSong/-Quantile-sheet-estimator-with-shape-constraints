#' f_sigmatau
#'
#' Construct the $\Sigma_tau$ diagonal matrix
#' given the dimension $K$
#'
#' @param k dimension of the diagonal matrix
#'
#' @return $K\times K \Sigma_tau$ diagonal matix,
#'          with $\M \Sigma_\tau_{ij} = 1$ if $i\geq j$ and $0$ otherwise.
#'
f_sigmatau <- function(k) {
    sigmatau <- matrix(0L, k, k)
    sigmatau[lower.tri(sigmatau, diag = T)] <- 1
    return(sigmatau)
}


#' f_sigma
#'
#' Construct the $\Sigma$ matrix
#' given the dimension $k_\tau, K_1$
#'
#' @param ktau dimension of the taus
#' @param k1 dimension of the first covariate
#'
#' @return $k_\tauK_1 \times k_\tauK_1 \Sigma$ diagonal block matix.
#'
f_sigma <- function(ktau, k1) {
    sigmatau <- f_sigmatau(ktau)
    return(kronecker(sigmatau, diag(k1)))
}

#' f_mate
#'
#' Construct the $E_{K_\tau}$ diagonal matrix
#' given the dimension $K_\tau$
#'
#' @param ktau dimension of the diagonal matrix
#'
#' @return $K_\tau \times K_\tau E_{K_\tau}$ diagonal matix,
#'          with $\M E_{11} = 1$ and $0$ anywhere else.
#'
f_mate <- function(ktau) {
    m_e <- matrix(0L, ktau, ktau)
    m_e[1, 1] <- 1
    return(m_e)
}

#' f_matf
#'
#' Construct the $(K_\tau-2) \times K_\tau$ $F_{K_\tau}$ matrix
#' given the dimension $K$
#'
#' @param ktau dimension of the diagonal matrix
#'
#' @return $(K_\tau-2)\times K_\tau F_{K_\tau}$,
#'          aka the first order differencing matrix,
#'          without the first row.
#'
f_matf <- function(ktau) {
    m_dif <- diff_matrix(ktau, 1)[-1, ]
    if (is.vector(m_dif)) {
        m_dif <- t(m_dif)
    }
    return(m_dif)
}

#' f_d11
#'
#' Construct the $K_{\tau}(K_1-2) \times K_{\tau}K_1$ $D_{11}$ matrix
#' given the dimension $K$
#'
#' @param ktau dimension of the taus
#' @param k1 dimension of the first covariate
#'
#' @return $K_{\tau}(K_1-2) \times K_{\tau}K_1$ $D_{11}$ matrix
#'
f_d11 <- function(ktau, k1) {
    kronecker(f_mate(ktau), diff_matrix(k1, 2))
}

#' f_d12
#'
#' Construct the $K_{\tau}(K_1-1) \times K_{\tau}K_1$ $D_{12}$ matrix
#' given the dimension $K$
#'
#' @param ktau dimension of the taus
#' @param k1 dimension of the first covariate
#'
#' @return $K_{\tau}(K_1-1) \times K_{\tau}K_1$ $D_{12}$ matrix
#'
f_d12 <- function(ktau, k1) {
    kronecker((diag(ktau) - f_mate(ktau)), diff_matrix(k1, 1))
}

#' f_dtau
#'
#' Construct the $K_{\tau - 2}K_1 \times K_{\tau}K_1$ $D_{\tau}$ matrix
#' given the dimension $K$
#'
#' @param ktau dimension of the taus
#' @param k1 dimension of the first covariate
#'
#' @return $K_{\tau - 2}K_1 \times K_{\tau}K_1$ $D_{\tau}$ matrix
#'
f_dtau <- function(ktau, k1) {
    kronecker(f_matf(ktau), diag(k1))
}

#' f_mats
#'
#' Construct the $K_1K_2 \times K_1K_2$ $S$ matrix
#' given the dimension $K_1, K_2$
#'
#' @param ktau dimension of the taus
#' @param k1 dimension of the first covariate
#' @param lambdas vector of lambdas(penalty parameters)
#'
#' @return $K_{K_\tau}K_1 \timesK_{K_\tau}K_1$ $S$ matrix
#'
f_mats <- function(ktau, k1, lambdas) {
    s <- lambdas[1] * crossprod(f_dtau(ktau, k1))
    s <- s + lambdas[2] * crossprod(f_d11(ktau, k1))
    s <- s + lambdas[3] * crossprod(f_d12(ktau, k1))
    return(s)
}

#' f_matc
#'
#' Construct the $K_{K_\tau}K_1 \times $K_{K_\tau}K_1$ diagonla matrix $C$
#' given the dimension \code{v_dim}
#'
#' @param v_dim vector contain the dimension of covariates
#'             $(K_\tau, K_1, K_2, \cdots)$
#'
#' @param til_beta $K_\tauK_1\cdots$ vector of the coefficient $\tilde{\beta}$
#'
#' @return K_{K_\tau}K_1\cdots \times K_{K_\tau}K_1\cdots$ diagonal matrix $C$,
#'         with diagonal equals derivative of $\tilde{\beta}$.
#'          $\M C_{jj} = \begin{cases}
#'                      1, &\text{if } \tilde{\beta}_j=\beta_j\\
#'                       \exp(\beta_j), &\text{otherwise}.
#'                      \end{cases}$
#'
f_matc <- function(v_dim, til_beta) {
    if (prod(v_dim) != length(til_beta)) {
        print("v_dim:"); print(v_dim); print("til_beta:"); print(til_beta)
        stop("f_matc: input dimension invalid \n")
    }
    v_idx <- seq_len(v_dim[2])
    til_beta[v_idx] <- 1L
    m_c <- diag(til_beta)
    return(m_c)
}


#' f_g1
#'
#' Construct the $(K_0+m) \times (K_0+m)$ or $K \times K$ $G_1$ diagonal matrix
#' given the dimension $K = K_0 + m$
#'
#' @param k0 number of internal knots
#' @param m order of the basis splines
#'
#' @return $K\times K G_1$ diagonal matix,
#'          with $G_1_{ii} = (t_{i+m} - t_i)/m$.
#'
f_g1 <- function(k0, m) {
    h <- 1 / (1 + k0)
    v_knots <- seq_len(k0) * h
    v_knots <- c(rep(0, m), v_knots, rep(1, m))
    out <- diag(diff(x = v_knots, lag = m) / m)
    return(out)
}


#' f_g2
#'
#' Construct the $(K_0+m) \times (K_0+m)$ or $K \times K$ $G_2$ diagonal matrix
#' given the dimension $K = K_0 + m$
#'
#' @param k0 number of internal knots
#' @param m order of the basis splines
#'
#' @return $K\times K G_2$ diagonal matix,
#'          with $G_2_{ii} = (t_{i+m+1} - t_i)/m$.
#'
f_g2 <- function(k0, m) {
    out <- f_g1(k0, m + 1)
    out <- out[-1, -1]
    return(out)
}



#' f_n1
#'
#' Construct the $n \times (K_1)$ spline matrix for the first covariate
#' given the dimension $K_1 = k1 + m1$
#'
#' @param x observations for the first covariate
#' @param k1 number of internal knots
#' @param m1 order of the basis splines
#'
#'
#' @return $n \times K_1$ spline matrix for the first covariate,
#'
f_n1 <- function(x, k1, m1) {
    out <- spline_basis(x, k1, m1)
    return(out)
}


#' f_ntau
#'
#' Construct the $n \times (K_1)$ spline designfor the taus
#' given the dimension $K_\tau = k\tau + m\tau$
#'
#' @param v_tau observations for the tau
#' @param ktau number of internal knots
#' @param mtau order of the basis splines
#'
#'
#' @return $n \times K_\tau$ spline design for the taus,
#'
f_ntau <- function(v_tau, ktau, mtau) {
    spline_basis(v_tau, ktau, mtau, x_range = c(1, 0))
}


#' f_ntau
#'
#' Construct the $n \times (K_\tau)$ spline matrix $N_{\tau^*}$
#' given the dimension $K_\tau = k\tau + m\tau$
#'
#' @param taustar estimated tau for each observations
#' @param k\tau number of internal knots
#' @param m\tau order of the basis splines
#'
#' @import splines
#'
#' @return $n \times K_\tau$ intergrated spline matrix $N_{\tau^*}$
#'
f_ntaustar <- function(taustar, ktau, mtau) {
    out <- spline_basis(taustar, ktau, mtau + 1, x_range = c(1, 0))
    out <- out[, -1]
    return(out)
}
