#' cq.reg
#'
#' Estimate conditional quantile in two steps:
#' 1. estimate the conditional probability for each observation: tau=p(y|x)
#'          Default to use local method over knn method,
#'          if both arguement supply or unspecify,
#'          only the local method will use.
#' 2. apply scam package to fix x to y with one more predictor tau,
#'    constrain tau to be monotonely nondecreasing.
#'
#' @param x vector or matrix whose rows contains predictors
#'          of each observations.
#' @param y vector of the variable we want to estimate its condtional cdf
#' @param span  fraction of total distance to be the bandwidth, between 0 and 1,
#'             default to be 0.5 if bandwidth not supply.
#' @param bandwidth the bandwidth of x to evaluate the cdf,
#'                 overwrite the span input.
#' @param rnn   fraction of total number of observations to be knn,
#'              between 0 and 1. Set span to be NULL to use knn,
#'              default to be 0.5, if not supply
#' @param knn   the number of nearest neighbor to estimate the cdf,
#'              should be less than the number of observation and greater than 0
#'              (\code{knn < length(y)}.
#'              setting this value overwrite the rnn
#' @param ...   additional arguement pass into scam (eg. sp, gamma, weights)
#'
#' @import scam
#'
#' @return a scam object
#'
#' @examples
#' devtools::load_all()
#' x <- soy$V5
#' y <- soy$V6
#' fit <- cq.reg(x, y)
#' cq_plot(fit)

cq.reg <- function(x, y, arg_bs = c("tesmi1", "ps"),
                opt = "efs", dims = c(5, 15), ords = NA,
                span = 0.05, bandwidth=NULL, rnn=NULL, knn=NULL,
                analysis=F, sp = NULL, ...) {
    #1. Estimate the conditional probability
    tau <- cond_pr(x, y, span, bandwidth, rnn, knn)

    if (is.vector(x)) {
        s_k <- 1
    } else if (is.matrix(x)) {
        s_k <- ncol(x)
    } else {
        stop("input x type not support")
    }

    xnam <- c("tau", paste0("x", seq_len(s_k)))
    newx <- cbind(tau, x)

    dat <- data.frame(y, newx)
    colnames(dat) <- c("y", xnam)

    #2. apply scam package
    fit <- scam(y ~ s(tau, x1, k = dims, m = ords, bs = arg_bs),
            data = dat, optimizer = opt, sp = sp, ...)

    if (analysis) {
        par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
        plot(fit, se = TRUE)
        plot(fit, pers = TRUE, theta = 30, phi = 40)
        plot(y, fit$fitted.values,
            xlab = "Real Data",
            ylab = "Fitted data", pch = ".", cex = 3)

        plot(x, y)
        sort_x <- sort(x)
        for (i in 1:10) {
            newx <- data.frame(i / 10, sort_x)
            colnames(newx) <- xnam
            lines(sort_x, predict(fit, newx), col = i + 1)
        }
    }

    fit$x <- x
    fit$xnam <- xnam
    return(fit)
}


#' @rdname cq.reg
#' @export
cq_regplot <- function(fit, v_taus = seq(0.05, 0.95, 0.15), ...) {
    x <- fit$x
    y <- fit$y
    xnam <- fit$xnam
    sort_x <- sort(x)

    xlim <- grDevices::extendrange(x, f = .01)
    ylim <- grDevices::extendrange(y, f = .01)
    xlab <- "x"
    ylab <- "y"

    plot(x, y, ...,
        type = "p", pch = 16, cex = 0.8, col = "grey",
        xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

    for (i in seq_along(v_taus)) {
        newx <- data.frame(v_taus[i], sort_x)
        colnames(newx) <- xnam
        lines(sort_x, predict(fit, newx), col = i + 1)
}

#    lgd <- c(
 #           expression(paste(mu, "=0.1")),
#            expression(paste(mu, "=0.2")),
#            expression(paste(mu, "=0.3")),
#            expression(paste(mu, "=0.4")),
#            expression(paste(mu, "=0.5")),
#            expression(paste(mu, "=0.6")),
#            expression(paste(mu, "=0.7")),
#            expression(paste(mu, "=0.8")),
#            expression(paste(mu, "=0.9")),
#            expression(paste(mu, "=1")))

#    legend("topleft",
#            lty = 1, legend = rev(lgd), col = 11:2,
#            lwd = 1, cex = 1, bty = "n"
#        )
}