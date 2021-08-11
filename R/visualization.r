#' cq_plot
#' Plot continuous quantile curves
#'
#' @param taus  vectors of specified conditional probability
#' @param con_quan an continuous quantile object includes
#'        \item{y}: vector of y
#'        \item{x}: vector of x
#'        \item{v_k0}:  vector of internal knots for covariates
#'        \item{v_m}:  vector of orders for covariates
#'        \item{m_sigma}:   matrix $\Sigma$
#'        \item{n1}:  design matrix of first covariate
#'        \item{til_beta}:  coefficient $\tilde{\beta}$
#'        \item{iter}   number of iteration
#'        \item{loss} Error or loss estimated
#'        \item{alpha}  accepted decrease in loss
#'        \item{bet}    shrinkage rate $\bet$, within $(0, 1)$.
#'
#' @keywords quantiles
#' @export
cq_plot <- function(taus, con_quan, ...) {
    #Extract values
    y <- con_quan$y; x <- con_quan$x; s_center <- con_quan$center
    s_sd <- con_quan$sd
    #transfer back
    x <- x * s_sd + s_center
    # Draw data
    plot(x, y, ...,
        type = "p", pch = 16, cex = 0.8, col = "grey",
        xlim = grDevices::extendrange(x, f = .1),
        ylim = grDevices::extendrange(y, f = .1),
        xlab = "", ylab = ""
    )
    # Calculate predictions
    v_pred <- cq_prediction(taus, con_quan)
    for (i in seq_along(taus)) {
        lines(x, v_pred[, i], col = i + 1)
    }
}

#' cq_prediction
#' prediction quantiles for orginal dataset given any conditional probability
#'
#' @param taus  vectors of specified conditional probability
#' @param con_quan an continuous quantile object includes
#'        \item{y}: vector of y
#'        \item{x}: vector of x
#'        \item{v_k0}:  vector of internal knots for covariates
#'        \item{v_m}:  vector of orders for covariates
#'        \item{m_sigma}:   matrix $\Sigma$
#'        \item{n1}:  design matrix of first covariate
#'        \item{til_beta}:  coefficient $\tilde{\beta}$
#'        \item{iter}   number of iteration
#'        \item{loss} Error or loss estimated
#'        \item{alpha}  accepted decrease in loss
#'        \item{bet}    shrinkage rate $\bet$, within $(0, 1)$.
#'
#' @param newdata a new dataset of first covariate
#'
#' @keywords quantiles
#' @export
cq_prediction <- function(taus, con_quan, newdata = NULL) {
    #Define dimensions
    k1 <- con_quan$v_k0[2]; ktau <- con_quan$v_k0[1]
    m1 <- con_quan$v_m[2]; mtau <- con_quan$v_m[1]
    #Extract values
    y <- con_quan$y; x <- con_quan$x
    m_sigma <- con_quan$m_sigma; til_beta <- con_quan$til_beta;
    s_center <- con_quan$center; s_sd <- con_quan$sd

    #Whether use newdata
    if (is.null(newdata)) {
        n1 <- con_quan$n1
        s_n <- length(y)
    } else {
        if (max(newdata) > max(x) || min(newdata) < min(x)) {
            stop("Can only prodict in original range.")
        }
        newdata <- (newdata - s_center) / s_sd
        n1 <- f_n1(scale(newdata), k1, m1)
        s_n <- length(newdata)
    }

    #Form design matrix
    m_n <- kronecker(spline_basis(taus, ktau, mtau, c(1, 0)), n1)

    # Calculate predictions
    v_pred <- c(m_n %*% m_sigma %*% til_beta)

    # Rearrange the orders
    m_pred <- matrix(v_pred, nrow = length(y), ncol = length(taus))

    return(m_pred)
}