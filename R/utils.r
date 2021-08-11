view_list <- function(list) {

  for (item in 1:length(list)) {

    print(head(list[[item]]))

  }
}


#' row-wise centering
#'
#' center the col of matrix x
#' @param x a matrix
#'
#' @keywords norm internal
center_colmeans <- function(x) {
  xcenter <- colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

#' The 2norm or Frobenius norm function
#'
#'
#' @param x vector or matrix to calculate its 2norm or Frobenius norm
#' @return 2norm or Frobenius norm
#'
#' @keywords norm internal
norm_2 <- function(x) sqrt(sum(x^2))



#' Differencing matrix
#'
#' @param k number of coeficient to be diffed
#' @param m order of differencing
#'
#' @keywords
#'
diff_matrix <- function(k, m) {
  m_diff <- diag(k)
  while (m > 0) {
    m_diff <- diff(m_diff)
    m <- m - 1
  }
  if (is.vector(m_diff)) {
    m_diff <- t(m_diff)
  }
  return(m_diff)
}



#' P-spline basis
#'
#' Construct the $n \times K$ spline matrix using quiv-distant knots
#' $K = k_0 + m$
#'
#' @param x observations for the first covariate
#' @param k0 number of internal knots
#' @param m order of the basis splines
#' @param x_range range of the spline basis,
#'                default to be NULL then use the c(max, min)
#'
#' @import splines
#'
#' @return $n \times K$ spline matrix for the first covariate,
#'
spline_basis <- function(x, k0, m, x_range=NULL) {
  if (is.null(x_range)) {
    x_range <- c(max(x), min(x))
  }

  big <- x_range[1]; small <- x_range[2]
  h <- (big - small) / (k0 + 1)
  v_knots <- seq_len(k0) * h + small
  v_knots <- c(rep(small, m), v_knots, rep(big, m))
  out <- splines::splineDesign(v_knots, x, ord = m)
  return(out)
}
