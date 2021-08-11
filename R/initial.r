#' f_init
#'
#' Construct an initial for $\beta$ by minimizing least square criteria
#'
#' @param v_y   vector of y
#' @param v_x   vector of x
#' @param m_sigma   $K_1K_2 \times K_1K_2 \Sigma$ diagonal block matix.
#' @param n1    $n \times K_1$ spline matrix for the first covariate,
#' @param v_k0  vector of internal knots for covariates
#' @param v_m   vector of orders for covariates
#' @param s     $K_1K_2 \times K_1K_2$ $S$ matrix
#'
#' @import mgcv
#'
#' @return an initial for $\beta$
#'
#' @examples
#' \dontrun{
#' devtools::load_all()
#' set.seed(202104)
#' k1 <- 4; k2 <- 4; m1 <- 4; m2 <- 4; n <- 400
#' v_x <- runif(n); v_y <- log(v_x); sigma <- 0.4 * (0.5 + v_x)
#' v_y <- v_y + rnorm(n, mean = rep(0, n), sd = sigma)
#' n1 <- f_n1(v_x, k1, m1)
#' m_sigma <- f_sigma(k1 + m1, k2 + m2
#' m_s <- f_mats(k1 + m1, k2 + m2, rep(10, 3))
#' f_init(v_y, v_x, m_sigma, n1, c(k1, k2), c(m1, m2), m_s)
#' }
#'
f_init <- function(v_y, v_x, m_sigma, n1, v_k0, v_m, s) {
    #Define dimensions
    k1 <- v_k0[2]; ktau <- v_k0[1]; m1 <- v_m[2]; mtau <- v_m[1]

    #1. Estimate the conditional probability
    v_taus <- cond_pr(v_x, v_y)
    m_ntau <- f_ntau(v_taus, ktau, mtau)
    #2. Construct Design matrix
    m_design <- mgcv::tensor.prod.model.matrix(list(m_ntau, n1)) %*% m_sigma

    #3. Construct Constraints
    v_idx <- seq_len(k1 + m1)
    a_in <- diag(prod(v_k0 + v_m)) * length(v_y)
    a_in <- a_in[-v_idx, ]
    b_in <- rep(0, nrow(a_in))

    l_pcls <- list(X     = m_design,
                    C    = matrix(0, 0, 0),
                    sp   = 1,
                    p    = rep(1, ncol(a_in)),
                    off  = 0,
                    S    = list(s),
                    Ain  = a_in,
                    bin  = b_in,
                    y    = v_y,
                    w    = v_y * 0 + 1)
    mgcv::pcls(l_pcls)
}
