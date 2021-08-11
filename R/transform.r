#' tilde_beta
#'
#' transformed coefficient $\tilde{\beta}$
#'
#' @param v_beta  coefficients for scop
#' @param k1 dimension of the first covariate
#'
#'
#' @return transformed coefficient $\tilde{\beta}$
#'
#' @examples
#' \dontrun{
#' k1 <- 4; ktau <- 4; m1 <- 4; mtau <- 4; y <- rnorm(200)
#' n1 <- f_n1(1:200, k1, m1); beta <- rnorm(prod(ktau + mtau, k1 + m1))
#' til_beta <- tilde_beta(beta, ktau + mtau, k1 + m1)
#' }

tilde_beta <- function(v_beta, k1) {
    v_idx <- seq_len(k1)
    v_beta[-v_idx] <- exp(v_beta[-v_idx])
    #cap the maximum values
    v_beta[which(v_beta == Inf)] <- 1e8
    return(v_beta)
}

#' og_beta
#'
#' transformed coefficient from $\tilde{\beta}$ to $\beta$
#'
#' @param til_beta  coefficients for scop
#' @param ktau dimension of the taus
#' @param k1 dimension of the first covariate
#'
#'
#' @return transformed coefficient $\beta$
#'
#' @examples
#' \dontrun{
#' k1 <- 4; ktau <- 4; m1 <- 4; mtau <- 4; y <- rnorm(200)
#' n1 <- f_n1(1:200, k1, m1); beta <- rnorm(prod(ktau + mtau, k1 + m1))
#' til_beta <- tilde_beta(beta, k1 + m1)
#' og_beta(til_beta, ktau + mtau, k1 + m1)
#' }

og_beta <- function(til_beta, ktau, k1) {
    #Index of til_betas
    v_idx <- seq_len(k1)
    v_idx <- seq_len(k1 * ktau)[-v_idx]
    #deal with NAs
    v_naidx <- intersect(v_idx, which(til_beta < 0))
    til_beta[v_naidx] <- 0
    #Transform
    til_beta[v_idx] <- log(til_beta[v_idx])
    # Deal with 0s
    til_beta[which(til_beta == -Inf)] <- -1e5
    return(til_beta)
}