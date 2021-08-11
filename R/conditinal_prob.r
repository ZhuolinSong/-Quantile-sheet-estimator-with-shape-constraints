#' cond.pr
#'
#' Estimate conditional probability for each observation
#' Default to use local method over knn method,
#' if both arguement supply or unspecify, only the local method will use.
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
#'
#' @return vector of conditional probability \code{p(Y < y | x)}
#'
#' @examples
#' cond_pr(1:5, 1:5) #[1] 0.33 0.5 0.6 0.75 1
#' cond_pr(1:5, 5:1, span=0) #[1] 1 1 1 1 1
#' cond_pr(1:5, 5:1, span=1) #[1] 1.0 0.8 0.6 0.4 0.2
#' cond_pr(1:5, 5:1, span=NULL) [1] 1. 0.66 0.66 0.66 0.33
#' cond_pr(1:5, 5:1, span=NULL, rnn=1) [1] 1.0 0.8 0.6 0.4 0.2
#' cond_pr(1:5, 5:1, span=NULL, rnn=0) [1] 1 1 1 1 1
#' cond_pr(1:5, 5:1, span=NULL, rnn=0, knn=1) [1] 1. 0.6 0.66 0.66 0.5
#' cond_pr(1:5, 5:1, span=NULL, rnn=0, knn=3) [1] 1.00 0.75 0.60 0.50 0.25
#' cond_pr(1:5, 5:1, bandwidth=1) [1] 1 0.6 0.6 0.6 0.5
#' cond_pr(matrix(1:12, nrow=6), 6:1) [1] 1. 0.75 0.6 0.6 0.5 0.33
#' cond_pr(matrix(1:12, nrow=6), 6:1, span = 0.6) [1] 1.0.8 0.67 0.5 0.4 0.25
#' cond_pr(matrix(1:12, nrow=6), 6:1, bandwidth = 5) 1. 0.8 0.6 0.5 0.4 0.25
#' cond_pr(matrix(1:12, nrow=6), 6:1, span = NULL) [1] 1. 0.75 0.6 0.6 0.5 0.25
#' cond_pr(matrix(1:12, nrow=6), 6:1, span = NULL, rnn = 0)
#' cond_pr(matrix(1:12, nrow=6), 6:1, span = NULL, rnn = 1)
#' cond_pr(matrix(1:12, nrow=6), 6:1, span = NULL, rnn = 1, knn=1)

cond_pr <- function(x, y,
                    span = 0.1, bandwidth=NULL,
                    rnn = NULL, knn=NULL) {
    #--------------------------------------
    ### Checking---------------------------
    #--------------------------------------
    if (is.vector(x)) {
        if ((s_n <- length(x)) != length(y)) {
            stop("input error: length(x)!=length(y)!")
        }
        is_vec <- T
    } else {
        if ((s_n <- nrow(x)) != length(y)) {
            stop("input error: nrow(x)!=length(y)!")
        }
        is_vec <- F
    }

    use_local <- F
    if (!is.null(span) || !is.null(bandwidth)) {
        if ((span < 0 || span > 1) & missing(bandwidth)) {
            stop("input error: span out of bound, no bandwidth supply!")
        }
        use_local <- T
    } else {
        if (missing(knn) || is.null(knn)) {
            if (is.null(rnn)) {
                rnn <- 0.5
            } else if (rnn < 0 || rnn > 1) {
                stop("input error: rnn out of bound, no knn supply!")
            }
            knn <- floor(rnn * s_n)
        }

        if (knn >= s_n) {
            print("warning: knn > s_n, set knn = s_n.")
            return(order(y) / s_n)
        }
    }

    #---------------------------------------------------
    ### Calculate probability---------------------------
    #---------------------------------------------------
    if (use_local) {
        output <- local_method(s_n, x, y, span, bandwidth, is_vec)
    } else {
        output <- knn_method(s_n, x, y, knn, is_vec)
    }

    return(output)
}



local_method <- function(s_n, x, y, span, bandwidth, is_vec) {
    # preallocate vector for condition probability
    c_prob <- c()
    if (is_vec) {
        s_range <- max(x) - min(x)
        if (missing(bandwidth) || is.null(bandwidth)) {
            bandwidth <- span * s_range
        }
        if (bandwidth >= s_range) {
            print("warning: bandwidth > maxwidth, bandwidth = maxwidth.")
            return(order(y) / s_n)
        }
        # Estimate conditional probability
        for (i in seq_len(s_n)) {
            idx <- which(abs(x - x[i]) <= bandwidth)
            c_prob[i] <- mean(y[idx] <= y[i])
        }
    } else {
        # Preallocate matrix for distance
        m_dist <- matrix(0, nrow = s_n, ncol = s_n)
        # Estimate distance between each x_i
        for (i in seq_len(s_n)) {
            diff <- x[- (1:i), ] - tcrossprod(rep(1, s_n - i), x[i, ])
            m_dist[i, - (1:i)] <- rowSums(diff^2)
        }
        m_dist <- sqrt(m_dist)
        s_range <- max(m_dist) - min(m_dist)
        if (missing(bandwidth) || is.null(bandwidth)) {
            bandwidth <- span * s_range
        }
        if (bandwidth >= s_range) {
            print("warning: bandwidth > maxwidth, bandwidth = maxwidth.")
            return(order(y) / s_n)
        }
        m_dist <- m_dist + t(m_dist)
        for (i in seq_len(s_n)) {
            idx <- which(m_dist[i, ] <= bandwidth)
            c_prob[i] <- mean(y[idx] <= y[i])
        }
    }

    return(c_prob)
}


knn_method <- function(s_n, x, y, knn, is_vec) {
    # preallocate vector for condition probability
    c_prob <- c()
    if (is_vec) {
        ord <- order(x)
        # Estimate conditional probability
        for (i in seq_len(s_n)) {
            ord_i <- which(ord == i)
            l_idx <- max(1, ord_i - knn)
            u_idx <- min(s_n, ord_i + knn)
            set <- ord[l_idx:u_idx]
            v_dist <- abs(x[set] - x[i])
            # find the closest k elements
            kth <- sort(v_dist, partial = knn + 1)[knn + 1]
            idx <- which(v_dist <= kth)
            c_prob[i] <- mean(y[set[idx]] <= y[i])
        }
    } else {
        # Estimate distance between each x_i
        for (i in seq_len(s_n)) {
            diff <- abs(x - tcrossprod(rep(1, s_n), x[i, ]))
            v_dist <- rowSums(diff^2)

            kth <- sort(v_dist, partial = knn + 1)[knn + 1]
            idx <- which(v_dist <= kth)
            c_prob[i] <- mean(y[idx] <= y[i])
        }
    }

    return(c_prob)
}