#' Quantile Sheets with pya's constraint
#'
#' Estimate conditional quantile sheet Schlossmacher's
#'  Alternating b/w weighted regression and recomputing weights
#' @param x vector or matrix whose rows contains predictors
#'          of each observations.
#' @param y vector of the variable we want to estimate its condtional cdf
#' @param arg_bs  argument for b-spline
#' @param dims knots number
#' @param ords  order of splien
#' @param opt   optimization methods
#' @param maxit maximum iteration for Schlossmacherl's method
#' @param ntaus number of taus to be considered
#' @param analysis  indicator whther plot analysis
#' @param trace
#' @param sp    smoothing coefficients
#' @param tol   tolerance for Schlossmacher's
#' @param ...   additional arguement pass into scam (eg. gamma, weights)
#'
#' @import scam
#'
#' @return a scam object
#'

qs_scam <- function(x, y, maxit = 10, ntaus = 10,
                arg_bs = c("tesmi1", "ps"),
                opt = "bfgs", dims = c(5, 15), ords = NA,
                analysis=F, trace=F, sp = NULL, tol = 5e-2,
                bnd = 1e-10, ...) {

    if (is.vector(x)) {
        s_k <- 1
    } else if (is.matrix(x)) {
        s_k <- ncol(x)
    } else {
        stop("input x type not support")
    }

    #initialization
    fit <- cq.reg(x, y, arg_bs, opt, dims, ords, sp = sp, ...)
    v_coef <- coef(fit)
    #fit <- f_descent(y, x, dims - ords, ords, maxit = 0)
    #v_coef <- fit$beta

    #repeat Y ntaus times
    ystar <- rep(y, ntaus)

    #data
    v_taus <- seq(bnd, 1 - bnd, length.out = ntaus)
    v_wog <- v_taus[cut(seq_along(ystar), breaks = ntaus, labels = FALSE)]
    xnam <- c("tau", paste0("x", seq_len(s_k)))
    newx <- cbind(v_wog, rep(x, ntaus))
    dat <- data.frame(ystar, newx)
    colnames(dat) <- c("y", xnam)

    #get info
    v_pred <- predict(fit, dat)
    #v_pred <- c(cq_prediction(v_taus, fit))
    v_resid <- ystar - v_pred

    #Schlossmacher's method
    for (iter in seq(maxit)) {
        #weights
        v_idx <- ystar <= v_pred
        v_weight <- v_wog
        v_weight[v_idx] <- 1 - v_weight[v_idx]
        v_weight <- v_weight / sqrt(v_resid^2 + 1e-8*max(abs(v_resid))^2)
        v_weight <- v_weight / mean(v_weight)
        #fit the model
        fit <- scam(y ~ s(tau, x1, k = dims, m = ords, bs = arg_bs),
            data = dat, optimizer = opt, sp = c(sp), weights = c(v_weight),
            start = NULL, ...)

        #get info
        v_pred <- fit$fitted.values
        v_resid <- fit$residual
        s_diff <- norm_2(v_coef - c(coef(fit))) / norm_2(v_coef)
        v_coef <- c(coef(fit))
        #plot the trace
        if (trace) {
            fit$x <- x; fit$y <- y; fit$xnam <- xnam
            sub <- paste("Iter:", iter, "diff:", round(s_diff, 3))
            main <- paste("Smoothing:", fit$sp, "dim:", dims)
            cq_regplot(fit, sub = sub, main = main)
        }
        #terminate
        if (s_diff <= tol) break
    }

    fit$x <- x; fit$y <- y; fit$xnam <- xnam
    if (analysis) {
        dev.new()
        par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
        plot(fit, se = TRUE)
        plot(fit, pers = TRUE, theta = 30, phi = 40)
        plot(ystar, fit$fitted.values,
            xlab = "Real Data",
            ylab = "Fitted data", pch = ".", cex = 3)
        cq_regplot(fit, sub = sub, main = main)
    }
    fit$ntaus <- ntaus
    fit$bnd <- bnd
    fit$dat <- dat
    fit$wog <- v_wog
    fit$idx <- v_idx
    return(fit)
}




#' Quantile Sheets(original)
#'
#' Estimate conditional quantile sheet Schlossmacher's
#'  Alternating b/w weighted regression and recomputing weights
#' @param x vector or matrix whose rows contains predictors
#'          of each observations.
#' @param y vector of the variable we want to estimate its condtional cdf
#' @param ords  order of splien
#' @param opt   optimization methods
#' @param maxit maximum iteration for Schlossmacherl's method
#' @param ntaus number of taus to be considered
#' @param analysis  indicator whther plot analysis
#' @param trace
#' @param sp    smoothing coefficients
#' @param tol   tolerance for Schlossmacher's
#' @param ...   additional arguement pass into scam (eg. gamma, weights)
#'
#' @import mgcv
#'
#' @return a gam object
#'

qs_gam <- function(x, y, maxit = 10, ntaus = 10,
                dims = c(10, 10), ords = NA,
                analysis=F, trace=F, sp = NULL, tol = 1e-3, ...) {

    if (is.vector(x)) {
        s_k <- 1
    } else if (is.matrix(x)) {
        s_k <- ncol(x)
    } else {
        stop("input x type not support")
    }
    #repeat Y ntaus times
    ystar <- rep(y, ntaus)

    #data
    v_taus <- seq(0, 1, length.out = ntaus)
    v_wog <- v_taus[cut(seq_along(ystar), breaks = ntaus, labels = FALSE)]
    xnam <- c("tau", paste0("x", seq_len(s_k)))
    newx <- cbind(v_wog, rep(x, ntaus))
    dat <- data.frame(ystar, newx)
    colnames(dat) <- c("y", xnam)

    #initialization
    fit <- mgcv::gam(y ~ te(tau, x1, k = dims, m = ords),
            data = dat, sp = c(sp), ...)
    v_coef <- coef(fit)

    #get info
    v_pred <- predict(fit, dat)
    v_resid <- ystar - v_pred

    #Schlossmacher's method
    for (iter in seq(maxit)) {
        #weights
        v_idx <- ystar <= v_pred
        v_weight <- v_wog
        v_weight[v_idx] <- 1 - v_weight[v_idx]
        v_weight <- v_weight / sqrt(v_resid^2 + 1e-8*max(abs(v_resid))^2)
        v_weight <- v_weight / mean(v_weight)
        #fit the model
        fit <- mgcv::gam(y ~ te(tau, x1, k = dims, m = ords),
            data = dat, sp = c(sp), weights = c(v_weight),
            start = v_coef, ...)

        #terminate
        s_diff <- norm_2(v_coef - coef(fit)) / norm_2(v_coef)
        v_coef <- coef(fit)
        #plot the trace
        if (trace) {
            fit$x <- x; fit$y <- y; fit$xnam <- xnam
            sub <- paste("Iter:", iter, "diff:", round(s_diff, 3))
            main <- paste("Smoothing:", fit$sp, "dim:", dims)
            cq_regplot(fit, sub = sub, main = main)
        }
        if (s_diff <= tol) break
        #get info
        v_pred <- fit$fitted.values
        v_resid <- ystar - v_pred
    }

    fit$x <- x; fit$y <- y; fit$xnam <- xnam
    if (analysis) {
        dev.new()
        par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
        plot(fit, se = TRUE)
        plot(fit, pers = TRUE, theta = 30, phi = 40)
        plot(ystar, fit$fitted.values,
            xlab = "Real Data",
            ylab = "Fitted data", pch = ".", cex = 3)
        cq_regplot(fit, sub = sub, main = main)
    }
    fit$ntaus <- ntaus
    fit$dat <- dat
    fit$wog <- v_wog
    fit$idx <- v_idx
    return(fit)
}


#' Expectile Sheets with pya's constraint
#'
#' Estimate conditional quantile sheet iterative evalution
#' @param x vector or matrix whose rows contains predictors
#'          of each observations.
#' @param y vector of the variable we want to estimate its condtional cdf
#' @param arg_bs  argument for b-spline
#' @param dims knots number
#' @param ords  order of splien
#' @param opt   optimization methods
#' @param maxit maximum iteration for Schlossmacherl's method
#' @param ntaus number of taus to be considered
#' @param analysis  indicator whther plot analysis
#' @param trace
#' @param sp    smoothing coefficients
#' @param tol   tolerance for Schlossmacher's
#' @param ...   additional arguement pass into scam (eg. gamma, weights)
#'
#' @import scam
#'
#' @return a scam object
#'
es_scam <- function(x, y, maxit = 10, ntaus = 10,
                arg_bs = c("tesmi1", "ps"),
                opt = "nlm.fd", dims = c(5, 10), ords = NA,
                analysis=F, trace=F, sp = NULL, tol = 1e-2, ...) {

    if (is.vector(x)) {
        s_k <- 1
    } else if (is.matrix(x)) {
        s_k <- ncol(x)
    } else {
        stop("input x type not support")
    }

    #initialization
    fit <- cq.reg(x, y, arg_bs, opt, dims, ords, sp = sp, ...)
    v_coef <- coef(fit)

    #repeat Y ntaus times
    ystar <- rep(y, ntaus)

    #data
    v_taus <- seq(0, 1, length.out = ntaus)
    v_wog <- v_taus[cut(seq_along(ystar), breaks = ntaus, labels = FALSE)]
    xnam <- c("tau", paste0("x", seq_len(s_k)))
    newx <- cbind(v_wog, rep(x, ntaus))
    dat <- data.frame(ystar, newx)
    colnames(dat) <- c("y", xnam)

    #get info
    v_pred <- predict(fit, dat)

    #Schlossmacher's method
    for (iter in seq(maxit)) {
        #weights
        v_idx <- ystar <= v_pred
        v_weight <- v_wog
        v_weight[v_idx] <- 1 - v_weight[v_idx]
        v_weight <- v_weight / mean(v_weight)
        #fit the model
        fit <- scam(y ~ s(tau, x1, k = dims, m = ords, bs = arg_bs),
            data = dat, optimizer = opt, sp = c(sp), weights = c(v_weight),
            start = NULL, ...)

        #terminate
        s_diff <- norm_2(v_coef - coef(fit)) / norm_2(v_coef)
        v_coef <- coef(fit)
        #plot the trace
        if (trace) {
            fit$x <- x; fit$y <- y; fit$xnam <- xnam
            sub <- paste("Iter:", iter, "diff:", round(s_diff, 3))
            main <- paste("Smoothing:", fit$sp, "dim:", dims)
            cq_regplot(fit, sub = sub, main = main)
        }
        if (s_diff <= tol) break
        #get info
        v_pred <- fit$fitted.values
    }

    fit$x <- x; fit$y <- y; fit$xnam <- xnam
    if (analysis) {
        dev.new()
        par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
        plot(fit, se = TRUE)
        plot(fit, pers = TRUE, theta = 30, phi = 40)
        plot(ystar, fit$fitted.values,
            xlab = "Real Data",
            ylab = "Fitted data", pch = ".", cex = 3)
        cq_regplot(fit, sub = sub, main = main)
    }
    fit$ntaus <- ntaus
    fit$dat <- dat
    fit$wog <- v_wog
    fit$idx <- v_idx
    return(fit)
}



#' Expectile Sheets(original)
#'
#' Expectile sheet 
#'  
#' @param x vector or matrix whose rows contains predictors
#'          of each observations.
#' @param y vector of the variable we want to estimate its condtional cdf
#' @param ords  order of splien
#' @param opt   optimization methods
#' @param maxit maximum iteration for Schlossmacherl's method
#' @param ntaus number of taus to be considered
#' @param analysis  indicator whther plot analysis
#' @param trace
#' @param sp    smoothing coefficients
#' @param tol   tolerance for Schlossmacher's
#' @param ...   additional arguement pass into scam (eg. gamma, weights)
#'
#' @import mgcv
#'
#' @return a gam object
#'

es_gam <- function(x, y, maxit = 10, ntaus = 10,
                dims = c(10, 10), ords = NA,
                analysis=F, trace=F, sp = NULL, tol = 1e-3, ...) {

    if (is.vector(x)) {
        s_k <- 1
    } else if (is.matrix(x)) {
        s_k <- ncol(x)
    } else {
        stop("input x type not support")
    }
    #repeat Y ntaus times
    ystar <- rep(y, ntaus)

    #data
    v_taus <- seq(0, 1, length.out = ntaus)
    v_wog <- v_taus[cut(seq_along(ystar), breaks = ntaus, labels = FALSE)]
    xnam <- c("tau", paste0("x", seq_len(s_k)))
    newx <- cbind(v_wog, rep(x, ntaus))
    dat <- data.frame(ystar, newx)
    colnames(dat) <- c("y", xnam)

    #initialization
    fit <- mgcv::gam(y ~ te(tau, x1, k = dims, m = ords),
            data = dat, sp = c(sp), ...)
    v_coef <- coef(fit)

    #get info
    v_pred <- predict(fit, dat)

    #Schlossmacher's method
    for (iter in seq(maxit)) {
        #weights
        v_idx <- ystar <= v_pred
        v_weight <- v_wog
        v_weight[v_idx] <- 1 - v_weight[v_idx]
        v_weight <- v_weight / mean(v_weight)
        #fit the model
        fit <- mgcv::gam(y ~ te(tau, x1, k = dims, m = ords),
            data = dat, sp = c(sp), weights = c(v_weight),
            start = v_coef, ...)

        #terminate
        s_diff <- norm_2(v_coef - coef(fit)) / norm_2(v_coef)
        v_coef <- coef(fit)
        #plot the trace
        if (trace) {
            fit$x <- x; fit$y <- y; fit$xnam <- xnam
            sub <- paste("Iter:", iter, "diff:", round(s_diff, 3))
            main <- paste("Smoothing:", fit$sp, "dim:", dims)
            cq_regplot(fit, sub = sub, main = main)
        }
        if (s_diff <= tol) break
        #get info
        v_pred <- fit$fitted.values
    }

    fit$x <- x; fit$y <- y; fit$xnam <- xnam
    if (analysis) {
        dev.new()
        par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
        plot(fit, se = TRUE)
        plot(fit, pers = TRUE, theta = 30, phi = 40)
        plot(ystar, fit$fitted.values,
            xlab = "Real Data",
            ylab = "Fitted data", pch = ".", cex = 3)
        cq_regplot(fit, sub = sub, main = main)
    }
    fit$ntaus <- ntaus
    fit$dat <- dat
    fit$wog <- v_wog
    fit$idx <- v_idx
    return(fit)
}