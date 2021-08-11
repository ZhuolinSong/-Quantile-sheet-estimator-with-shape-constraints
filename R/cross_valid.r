#' Cross Validationg
#' for Quantile Sheets with pya's constraint
#'
#'
#' @param x vector or matrix whose rows contains predictors
#'          of each observations.
#' @param y vector of the variable we want to estimate its condtional cdf
#' @param k k-fold cross-validation
#' @param xsp spar range(from 0 to 1):
#'              either a vector of tau smoothing parameters
#'              or an integer indicating the range of the selection range
#'              : (1 - 0.3) / range * (0:range) + 0.3
#' @param tausp spar range(from 0 to 1):
#'              either a vector of tau smoothing parameters
#'              or an integer indicating the range of the selection range
#'              : (1 - 0.3) / range * (0:range) + 0.3
#' 
#' @param ...   additional arguement pass into qs_scam
#' (eg. maxit = 10, ntaus = 10, arg_bs = c("tesmi1", "ps"),
#' opt = "bfgs", dims = c(5, 10), ords = NA,
#' analysis = F, trace=F, sp = NULL, tol = 1e-2,)
#'
#' @import scam
#'
#' @return a scam object
#'

cv_qsscam <- function(x, y, k = 5, xsp=10^(-5:0), tausp=1, ...) {
    s_n <- length(y) # Define the number of observation
    if (s_n < k) {
        stop(paste("folds:", k, "is greater than n:", s_n))
    }
    v_tausp <- tausp
    
    if (length(xsp) == 1) {
        v_xsp <- (1 - 0.3) / xsp * (0:xsp) + 0.3
    } else {
        v_xsp <- xsp
    }

    loss <- matrix(0, nrow = length(v_tausp), ncol = length(v_xsp))
    # Randomly shuffle the data
    y <- y[sample(s_n)]
    v_folds <- cut(seq(s_n), breaks = k, labels = FALSE) # Create k folds

    for (fold in 1:k) { # Perform k fold cross validation
        test_idx <- which(v_folds == fold, arr.ind = TRUE) # Segement data by fold
        test_y <- y[test_idx]; test_x <- x[test_idx]
        train_y <- y[-test_idx]; train_x <- x[-test_idx];
        # Loop the range
        for (i in seq_along(v_tausp)) {
            tausp <- v_tausp[i]
            for (j in seq_along(v_xsp)) {
                xsp <- v_xsp[j]
                # Build the model
                fit <- qs_scam(train_x, train_y,
                        sp = c(tausp, xsp), ...)
                # Store the loss
                loss[i, j] <- loss[i, j] + qs_loss(test_x, test_y, fit) / k
            }
        }
    }
    cv_idx <- arrayInd(which.min(loss), dim(loss))
    tausp <- v_tausp[cv_idx[1]]
    xsp <- v_xsp[cv_idx[2]]
    return(list(sel = c(tausp, xsp),
            loss = loss))
}


#' qs_loss
#'
#' Return the quantile loss criteria
#'
#' @param x vector or matrix whose rows contains predictors
#'          of each observations.
#' @param y vector of the variable we want to estimate its condtional cdf
#' @param fit vector of y produced by qs_scam, qs_gam, ex_scam or ex_gam
#'
#'
#' @return the quantile sheet loss
#'
qs_loss <- function(x, y, fit) {
    if (is.vector(x)) {
        s_k <- 1
    } else if (is.matrix(x)) {
        s_k <- ncol(x)
    } else {
        stop("input x type not support")
    }
    ntaus <- fit$ntaus
    bnd <- fit$bnd
    # form data
    #repeat Y ntaus times
    ystar <- rep(y, ntaus)
    #data
    v_taus <- seq(bnd, 1 - bnd, length.out = ntaus)
    v_wog <- v_taus[cut(seq_along(ystar), breaks = ntaus, labels = FALSE)]
    xnam <- c("tau", paste0("x", seq_len(s_k)))
    newx <- cbind(v_wog, rep(x, ntaus))
    dat <- data.frame(ystar, newx)
    colnames(dat) <- c("y", xnam)
    #reweight
    v_pred <- c(predict(fit, dat))
    v_idx <- ystar <= v_pred
    v_weight <- v_wog
    v_weight[v_idx] <- 1 - v_weight[v_idx]

    #calculate loss
    sum(abs(ystar - v_pred) * v_weight)
}