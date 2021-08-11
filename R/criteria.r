#' f_loss
#'
#' Return the quantile loss criteria
#'
#' @param til_beta the current estimation fro $\tilde{\beta}$
#' @param taustar vector of estimated tau for each observations
#' @param v_y vector of y
#'
#' @param m_sigma $K_1K_2 \times K_1K_2 \Sigma$ diagonal block matix.
#' @param h1 $K_1K_2$ vector $H_1$
#' @param h2 $K_1K_2$ vector $H_2$
#'
#' @return the quantile loss
#'
f_loss <- function(til_beta, taustar, v_y, m_sigma, h1, h2) {
    loss <- crossprod((taustar - 0.5), v_y)
    loss <- loss - crossprod((h1 + h2), m_sigma %*% til_beta)
    loss <- loss / length(v_y)
    return(c(loss))
}


#' f_criteria
#'
#' Return the QR criteria
#'
#' @param beta the current estimation for $\beta$
#' @param til_beta the current estimation fro $\tilde{\beta}$
#' @param taustar vector of estimated tau for each observations
#' @param v_y vector of y
#'
#' @param m_sigma $K_1K_2 \times K_1K_2 \Sigma$ diagonal block matix.
#' @param h1 $K_1K_2$ vector $H_1$
#' @param h2 $K_1K_2$ vector $H_2$
#' @param s $K_1K_2 \times K_1K_2$ $S$ matrix
#'
#' @return the quantile loss
#'
f_criteria <- function(beta, til_beta, taustar, v_y, m_sigma, h1, h2, s) {
    loss <- f_loss(til_beta, taustar, v_y, m_sigma, h1, h2)
    criteria <- loss + crossprod(beta, s %*% beta)
    return(c(criteria))
}